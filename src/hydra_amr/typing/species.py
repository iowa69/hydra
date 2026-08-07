"""Species identification and mapping onto AMRFinderPlus organism groups.

Two independent signals are combined. The primary one is the MLST result: seven
housekeeping loci matching a scheme's alleles at full length is species-level
evidence, and Hydra ships a curated scheme -> organism table because the map
distributed with some MLST installations is stale. When Mash sketches are
installed they provide a second, assembly-wide signal that can confirm or
override a weak MLST call.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config import BUNDLED_DATA
from ..db.manager import DatabaseStore
from ..records import MlstCall, SpeciesCall
from ..utils import LOG, have, run

#: Mash distance below which a sketch match is considered same-species.
MASH_SPECIES_DISTANCE = 0.05
MASH_STRONG_DISTANCE = 0.02


@dataclass(frozen=True)
class SchemeOrganism:
    genus: str
    species: str
    organism: str  # AMRFinderPlus taxgroup, may be empty

    @property
    def name(self) -> str:
        if self.genus and self.species:
            return f"{self.genus} {self.species}"
        if self.genus:
            return f"{self.genus} sp."
        return "unknown"


def load_scheme_table(store: DatabaseStore | None = None) -> dict[str, SchemeOrganism]:
    """Curated scheme -> organism table, extended by the installed MLST map."""
    table: dict[str, SchemeOrganism] = {}
    curated = BUNDLED_DATA / "scheme_species.tsv"
    if curated.exists():
        with open(curated) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                parts += [""] * (4 - len(parts))
                scheme = parts[0].strip()
                if not scheme:
                    continue
                table[scheme] = SchemeOrganism(parts[1].strip(), parts[2].strip(), parts[3].strip())
    if store is not None:
        try:
            entry = store.require_installed("pubmlst")
        except Exception:
            return table
        installed_map = store.root / entry["path"] / "scheme_species_map.tab"
        if installed_map.exists():
            with open(installed_map) as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    parts += [""] * (3 - len(parts))
                    scheme = parts[0].strip()
                    # Curated entries win; the installed map only fills gaps.
                    if scheme and scheme not in table:
                        table[scheme] = SchemeOrganism(parts[1].strip(), parts[2].strip(), "")
    return _inherit_organism_by_species(table)


def _inherit_organism_by_species(table: dict[str, SchemeOrganism]) -> dict[str, SchemeOrganism]:
    """Give a scheme the AMRFinderPlus organism curated for its species.

    The curated table is keyed by scheme name, but scheme names are a naming
    convention rather than data: a store built straight from PubMLST can hold a
    scheme the curated table has never seen, or the same scheme under a different
    name. What identifies the organism is the genus and species, so any entry that
    knows those but not its taxgroup borrows it from a curated entry that does.
    Without this, an unrecognised scheme name silently disables point mutations
    for that isolate.
    """
    by_species: dict[tuple[str, str], str] = {}
    for entry in table.values():
        if entry.organism and entry.genus:
            by_species.setdefault((entry.genus, entry.species), entry.organism)
    # A genus-only fallback covers schemes named for a species complex.
    by_genus: dict[str, str] = {}
    for (genus, _species), organism in by_species.items():
        if by_genus.setdefault(genus, organism) != organism:
            by_genus[genus] = ""  # ambiguous within the genus, so do not guess

    for scheme, entry in table.items():
        if entry.organism or not entry.genus:
            continue
        organism = by_species.get((entry.genus, entry.species)) or by_genus.get(entry.genus, "")
        if organism:
            table[scheme] = SchemeOrganism(entry.genus, entry.species, organism)
    return table


class SpeciesIdentifier:
    def __init__(self, store: DatabaseStore):
        self.store = store
        self.scheme_table = load_scheme_table(store)
        self.sketches = self._find_sketches()

    def _find_sketches(self) -> list[Path]:
        try:
            entry = self.store.require_installed("species")
        except Exception:
            return []
        path = self.store.root / entry["path"]
        return sorted(path.glob("*.msh"))

    def organism_for_scheme(self, scheme: str) -> SchemeOrganism | None:
        return self.scheme_table.get(scheme)

    def sketch_species(self, assembly: Path, threads: int = 1) -> SpeciesCall | None:
        """Species from Mash sketches alone, independent of any MLST result.

        Run before MLST so the genus can steer scheme selection: several PubMLST
        schemes now carry alleles from neighbouring genera, and a scheme chosen
        on allele counts alone can land on the wrong organism.
        """
        if not self.sketches or not have("mash"):
            return None
        return self._mash(assembly, threads)

    def identify(self, sample: str, mlst: MlstCall, assembly: Path | None = None,
                 threads: int = 1, sketch: SpeciesCall | None = None) -> SpeciesCall:
        """Combine the MLST call with an optional Mash screen."""
        call = SpeciesCall()
        entry = self.scheme_table.get(mlst.scheme) if mlst.scheme != "-" else None
        if entry is not None and entry.genus:
            fraction = (mlst.loci_found / mlst.loci_total) if mlst.loci_total else 0.0
            if mlst.sequence_type not in ("-", "") or fraction >= 0.85:
                confidence = "strong"
            elif fraction >= 0.5:
                confidence = "good"
            else:
                confidence = "weak"
            call = SpeciesCall(
                name=entry.name, genus=entry.genus, species=entry.species,
                confidence=confidence,
                evidence=f"MLST scheme {mlst.scheme} ({mlst.loci_found}/{mlst.loci_total} exact loci)",
                organism=entry.organism or None,
            )
        if sketch is None and assembly is not None and self.sketches and have("mash"):
            sketch = self._mash(assembly, threads)
        if sketch is not None:
            # A whole-genome sketch resolves species inside a complex that MLST
            # cannot: the Klebsiella pneumoniae scheme also types K. variicola
            # and K. quasipneumoniae, and reporting all three as K. pneumoniae
            # loses the distinction the scheme was never able to make.
            disagrees = sketch.name.lower() != call.name.lower()
            confident = sketch.distance is not None and sketch.distance < MASH_STRONG_DISTANCE
            # A disagreement about *genus* is settled by the sketch, not the scheme,
            # and at the looser species distance rather than the override distance.
            #
            # The two are not equal evidence. A sketch compares the whole genome; a
            # scheme compares seven or eight housekeeping loci, and several PubMLST
            # schemes have accumulated alleles matching neighbouring genera -- the
            # EnteroBase E. coli scheme types some Klebsiella at 8/8 exact loci. When
            # the two name different genera, the locus count is the thing that should
            # give way. Inside one genus the old rule stands: resolving K. variicola
            # from K. pneumoniae is what the strict threshold is for.
            genus_conflict = bool(
                sketch.genus and call.genus
                and sketch.genus.lower() != call.genus.lower()
                and sketch.distance is not None
                and sketch.distance < MASH_SPECIES_DISTANCE)
            if call.confidence in ("none", "weak") or (confident and disagrees) or genus_conflict:
                if call.confidence != "none":
                    sketch.evidence += f"; MLST scheme {mlst.scheme} suggested {call.name}"
                if genus_conflict and not confident:
                    # Won on genus, but by a sketch that is not itself strong. Say so
                    # rather than replacing one overconfident answer with another.
                    sketch.confidence = "good" if sketch.confidence == "strong" else sketch.confidence
                    sketch.evidence += (" (genus taken from the sketch: a scheme match"
                                        " cannot outvote a whole-genome distance)")
                # Never inherit the rejected call's organism: running the
                # AMRFinderPlus rules for the wrong taxgroup is worse than
                # running none at all.
                return sketch
            distance = f" (d={sketch.distance:.4f})" if sketch.distance is not None else ""
            call.evidence += f"; Mash {sketch.name}{distance}"
            # Two independent methods naming different genera is not strong evidence
            # for either of them. The scheme still wins the call -- a sketch that
            # missed MASH_STRONG_DISTANCE has not earned an override -- but the
            # confidence must not claim more than the evidence supports.
            #
            # This is not hypothetical. One isolate here matched the EnteroBase
            # E. coli scheme at 8/8 exact loci while its sketch said Klebsiella
            # pneumoniae at d=0.0210, a thousandth outside the override threshold.
            # It was reported as Escherichia coli with "strong" confidence, and the
            # E. coli mutation catalogue was applied to a Klebsiella genome.
            same_genus = (sketch.genus or "").lower() == (call.genus or "").lower()
            if disagrees and not same_genus and sketch.genus and call.genus:
                call.confidence = "weak"
                call.evidence += (f"; genus disagreement, sketch says {sketch.genus}"
                                  f" and the scheme says {call.genus}")
        return call

    def identify_reads(self, reads: Path, threads: int = 1) -> SpeciesCall | None:
        """Species from raw reads, so a FASTQ-only sample still gets an organism.

        Without this, read-only input cannot reach the organism-specific mutation
        catalogues at all and the user has to supply ``--organism`` by hand.
        """
        if not self.sketches or not have("mash"):
            return None
        return self._mash(reads, threads, from_reads=True)

    def _mash(self, query: Path, threads: int, from_reads: bool = False) -> SpeciesCall | None:
        best: tuple[float, str] | None = None
        for sketch in self.sketches:
            command = ["mash", "dist", "-p", str(max(1, threads))]
            if from_reads:
                # Treat the input as reads and ignore k-mers seen once, which are
                # almost all sequencing error.
                command += ["-r", "-m", "2"]
            command += [str(sketch), str(query)]
            try:
                proc = run(command, check=False)
            except Exception as exc:  # noqa: BLE001 - mash is optional
                LOG.debug("mash failed on %s: %s", sketch, exc)
                return None
            if proc.returncode != 0:
                LOG.debug("mash dist returned %d for %s", proc.returncode, sketch)
                continue
            for line in (proc.stdout or "").splitlines():
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                try:
                    distance = float(parts[2])
                except ValueError:
                    continue
                if best is None or distance < best[0]:
                    best = (distance, parts[0])
        if best is None or best[0] > MASH_SPECIES_DISTANCE:
            return None
        distance, reference = best
        name = _clean_sketch_name(reference)
        if not name:
            LOG.debug("mash matched %s but it carries no taxon name; ignoring", reference)
            return None
        genus, _, species = name.partition(" ")
        return SpeciesCall(
            name=name, genus=genus, species=species,
            confidence="strong" if distance < MASH_STRONG_DISTANCE else "good",
            evidence=f"Mash sketch {Path(reference).name}", distance=distance,
            organism=self.organism_for_name(name),
        )

    def organism_for_name(self, name: str) -> str | None:
        """Map a species name onto the AMRFinderPlus taxgroup that covers it."""
        for entry in self.scheme_table.values():
            if entry.organism and entry.name.lower() == name.lower():
                return entry.organism
        genus = name.split()[0] if name else ""
        # Taxgroups that cover a whole genus, plus the complexes whose scheme
        # entry is not a plain binomial and so never matches by name.
        genus_only = {
            "Escherichia": "Escherichia", "Shigella": "Escherichia",
            "Salmonella": "Salmonella", "Campylobacter": "Campylobacter",
            "Burkholderia": None, "Klebsiella": None,
        }
        if genus == "Burkholderia":
            species = name.split()[1].lower() if len(name.split()) > 1 else ""
            return {"cepacia": "Burkholderia_cepacia", "mallei": "Burkholderia_mallei",
                    "pseudomallei": "Burkholderia_pseudomallei"}.get(species)
        if genus == "Klebsiella":
            species = name.split()[1].lower() if len(name.split()) > 1 else ""
            # The Kp-complex members share the K. pneumoniae curation.
            if species in ("pneumoniae", "variicola", "quasipneumoniae",
                           "quasivariicola", "africana"):
                return "Klebsiella_pneumoniae"
            if species in ("oxytoca", "michiganensis", "grimontii", "pasteurii"):
                return "Klebsiella_oxytoca"
            return None
        return genus_only.get(genus)

    def known_organisms(self) -> list[str]:
        return sorted({e.organism for e in self.scheme_table.values() if e.organism})


def _clean_sketch_name(reference: str) -> str:
    """Turn a sketch reference id into a readable ``Genus species`` label.

    Sketches are named ``Genus_species/GCF_000000000.fna.gz``: the taxon lives in
    the directory component, and the file itself is only an assembly accession.
    Returns an empty string when no taxon name can be recovered, so the caller
    can fall back to another line of evidence rather than report an accession.
    """
    parts = [p for p in str(reference).replace("\\", "/").split("/") if p]
    taxon = parts[-2] if len(parts) >= 2 else parts[-1] if parts else ""
    for suffix in (".fna.gz", ".fasta.gz", ".fna", ".fasta", ".fa", ".gz"):
        if taxon.endswith(suffix):
            taxon = taxon[: -len(suffix)]
            break
    tokens = [t for t in taxon.replace("_", " ").split() if t]
    if not tokens:
        return ""
    genus = tokens[0]
    # Assembly accessions are not taxon names.
    if genus.upper().startswith(("GCF", "GCA")) or not genus[0].isalpha():
        return ""
    species = tokens[1].lower() if len(tokens) >= 2 else ""
    if species in ("unknown", "sp", "spp"):
        species = "sp."
    return f"{genus.capitalize()} {species}".strip()
