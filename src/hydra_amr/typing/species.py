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
            if call.confidence in ("none", "weak") or (confident and disagrees):
                if call.confidence != "none":
                    sketch.evidence += f"; MLST scheme {mlst.scheme} suggested {call.name}"
                # Never inherit the rejected call's organism: running the
                # AMRFinderPlus rules for the wrong taxgroup is worse than
                # running none at all.
                return sketch
            distance = f" (d={sketch.distance:.4f})" if sketch.distance is not None else ""
            call.evidence += f"; Mash {sketch.name}{distance}"
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
