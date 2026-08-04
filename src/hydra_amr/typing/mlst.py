"""Automatic multi-locus sequence typing.

Hydra does not ask which scheme to use: every sample in a batch is searched
against a single BLAST database holding the alleles of all installed PubMLST
schemes, and the scheme that recovers the most loci at full length wins. The
sequence type is then read from that scheme's profile table, with the usual
conventions for inexact calls (``~`` novel allele, ``?`` partial, ``-`` absent).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..records import MlstCall
from ..seqio import read_fasta, write_fasta
from ..utils import LOG, HydraError
from ..engines.blast import blast, check_db_exists, merge_hsps
from ..engines.nucl import QueryBatch

#: An allele must reach these to count as an exact match.
EXACT_IDENTITY = 100.0
EXACT_COVERAGE = 100.0
#: Below this identity a locus is not considered present at all.
MIN_LOCUS_IDENTITY = 90.0
MIN_LOCUS_COVERAGE = 60.0
#: Read depth tapers at the ends of a reference, so a locus may still be called
#: exact with this fraction of its length uncovered.
MAX_UNCALLED_FRACTION = 0.02
#: Intermediate-frequency sites in one locus before it is treated as mixed.
MIXED_SITES_PER_LOCUS = 3


@dataclass
class LocusHit:
    locus: str
    allele: str
    identity: float
    coverage: float
    bitscore: float

    @property
    def exact(self) -> bool:
        return self.identity >= EXACT_IDENTITY and self.coverage >= EXACT_COVERAGE

    @property
    def full_length(self) -> bool:
        return self.coverage >= EXACT_COVERAGE

    @property
    def allele_number(self) -> int:
        """Numeric allele id, or a large value for non-numeric ids."""
        try:
            return int(self.allele)
        except ValueError:
            return 1 << 30


class SchemeProfiles:
    """Lazy access to a scheme's ST profile table."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self._cache: dict[str, tuple[list[str], dict[tuple[str, ...], str], dict[str, str]]] = {}

    def load(self, scheme: str):
        """Return (loci, profile table, per-ST metadata) for one scheme.

        Which header columns are loci is decided by the scheme's own allele
        files, not by guessing which trailing column names look like metadata:
        PubMLST adds columns such as ``MLST_cluster`` between the loci and
        ``clonal_complex``, and a name-based guess shifts every allele by one,
        silently making the scheme unable to assign any ST.
        """
        if scheme in self._cache:
            return self._cache[scheme]
        scheme_dir = self.root / "pubmlst" / scheme
        path = scheme_dir / f"{scheme}.txt"
        known_loci = {p.stem for p in scheme_dir.glob("*.tfa")} if scheme_dir.is_dir() else set()
        loci: list[str] = []
        table: dict[tuple[str, ...], str] = {}
        extra: dict[str, str] = {}
        if path.exists():
            with open(path) as handle:
                header = handle.readline().rstrip("\n").split("\t")
                if known_loci:
                    columns = [(i, name) for i, name in enumerate(header)
                               if i > 0 and name in known_loci]
                else:
                    # No allele files to consult: fall back to dropping the
                    # metadata column names PubMLST is known to use.
                    meta_names = {"clonal_complex", "species", "mlst_clade", "lineage",
                                  "cc", "clade", "mlst_cluster", "subspecies"}
                    columns = [(i, name) for i, name in enumerate(header)
                               if i > 0 and name.lower() not in meta_names]
                loci = [name for _i, name in columns]
                indexes = [i for i, _name in columns]
                meta_index = None
                for i, name in enumerate(header):
                    if name.lower() in ("clonal_complex", "cc"):
                        meta_index = i
                        break
                for line in handle:
                    if not line.strip():
                        continue
                    row = line.rstrip("\n").split("\t")
                    if len(row) <= max(indexes, default=0):
                        continue
                    st = row[0]
                    table[tuple(row[i] for i in indexes)] = st
                    if meta_index is not None and len(row) > meta_index:
                        value = row[meta_index].strip()
                        if value:
                            extra[st] = value
        self._cache[scheme] = (loci, table, extra)
        return self._cache[scheme]

    def loci(self, scheme: str) -> list[str]:
        return self.load(scheme)[0]


class MlstTyper:
    """Runs one BLAST pass over all samples and assigns a scheme + ST to each."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config
        entry = store.require_installed("pubmlst")
        self.root = store.root / entry["path"]
        self.profiles = SchemeProfiles(self.root)
        self.fasta = self.root / "blast" / "mlst.fna"
        self._scheme_genus: dict[str, str] = {}
        self._allele_cache: dict[tuple[str, str], dict[str, str]] = {}

    def type_batch(self, batch: QueryBatch, workdir: Path, threads: int | None = None,
                   force_scheme: str | None = None,
                   exclude_schemes: set[str] | None = None,
                   genus_by_sample: dict[str, str] | None = None,
                   scheme_genus: dict[str, str] | None = None) -> dict[str, MlstCall]:
        check_db_exists(self.fasta, "nucl")
        threads = threads or self.config.threads
        out_tab = workdir / "mlst.blastn.tsv"
        LOG.info("MLST: searching %d contigs against all installed schemes", batch.n_contigs)
        hsps = blast(
            "blastn", batch.path, self.fasta, out_tab,
            threads=threads, evalue=1e-20, task="megablast", word_size=20,
            perc_identity=MIN_LOCUS_IDENTITY - 5.0, max_target_seqs=10000,
        )
        merged = merge_hsps(hsps)

        # sample -> scheme -> locus -> best LocusHit
        found: dict[str, dict[str, dict[str, LocusHit]]] = defaultdict(lambda: defaultdict(dict))
        for hit in merged:
            if hit.identity_pct < MIN_LOCUS_IDENTITY or hit.coverage_pct < MIN_LOCUS_COVERAGE:
                continue
            piece = batch.id_map.get(hit.qseqid)
            if piece is None:
                continue
            parts = hit.sseqid.split(".")
            if len(parts) < 3:
                continue
            scheme, locus, allele = parts[0], parts[1], ".".join(parts[2:])
            if exclude_schemes and scheme in exclude_schemes:
                continue
            if force_scheme and scheme != force_scheme:
                continue
            sample = piece.sample
            candidate = LocusHit(locus=locus, allele=allele, identity=hit.identity_pct,
                                 coverage=hit.coverage_pct, bitscore=hit.bitscore)
            current = found[sample][scheme].get(locus)
            if current is None or _better(candidate, current):
                found[sample][scheme][locus] = candidate

        results: dict[str, MlstCall] = {}
        genus_by_sample = genus_by_sample or {}
        self._scheme_genus = scheme_genus or {}
        for sample in sorted(batch.samples()):
            results[sample] = self._call_sample(found.get(sample, {}), force_scheme,
                                                genus_by_sample.get(sample, ""))
        return results

    def _call_sample(self, by_scheme: dict[str, dict[str, LocusHit]],
                     force_scheme: str | None, genus: str = "") -> MlstCall:
        if not by_scheme:
            return MlstCall(note="no MLST loci detected")
        ranked: list[tuple[tuple, str, dict[str, LocusHit]]] = []
        for scheme, loci_hits in by_scheme.items():
            scheme_loci = self.profiles.loci(scheme)
            total = len(scheme_loci) or len(loci_hits)
            exact = sum(1 for h in loci_hits.values() if h.exact)
            full = sum(1 for h in loci_hits.values() if h.full_length)
            mean_identity = sum(h.identity for h in loci_hits.values()) / len(loci_hits)
            # Rank by exact loci first, then full-length loci, then mean identity.
            # The scheme name breaks any remaining tie so the result does not
            # depend on the order BLAST happened to report its hits in: the
            # salmonella and senterica_achtman_2 schemes, for one, hold
            # byte-identical alleles at every locus and always tie.
            ranked.append(((exact, full, mean_identity, -abs(total - len(loci_hits))),
                           scheme, loci_hits))
        ranked.sort(key=lambda item: (item[0], _reverse_name(item[1])), reverse=True)
        candidates = [
            {"scheme": s, "exact_loci": sc[0], "loci": len(h), "mean_identity": round(sc[2], 2)}
            for sc, s, h in ranked[:5]
        ]
        best_score, best_scheme, best_hits = ranked[0]
        best_loci = self.profiles.loci(best_scheme) or sorted(best_hits)
        # A forced scheme still needs enough loci to mean anything; two exact
        # loci out of seven is not a sequence type.
        if best_score[0] < max(2, len(best_loci) // 2):
            return MlstCall(scheme="-", sequence_type="-", loci_found=best_score[0],
                            loci_total=len(best_loci), candidates=candidates,
                            note=f"no scheme reached half its loci (best: {best_scheme}, "
                                 f"{best_score[0]}/{len(best_loci)} exact)")

        # Several schemes can fit the same genome: neighbouring genera share
        # housekeeping alleles, and some species have both a legacy and a current
        # scheme over the same loci. Evaluate the plausible ones, then prefer the
        # scheme that agrees with the species call and actually resolves an ST.
        shortlist = [item for item in ranked[:6]
                     if item[0][0] >= best_score[0] - 1 or item[1] == force_scheme]
        evaluated = [self._evaluate_scheme(scheme, loci_hits)
                     for _score, scheme, loci_hits in shortlist]

        def preference(call: MlstCall) -> tuple:
            scheme_genus = self._scheme_genus.get(call.scheme, "")
            genus_match = 1 if (genus and scheme_genus and genus == scheme_genus) else 0
            genus_conflict = 1 if (genus and scheme_genus and genus != scheme_genus) else 0
            return (genus_match, -genus_conflict, 1 if call.sequence_type != "-" else 0,
                    call.loci_found, call.similarity)

        chosen = max(evaluated, key=preference)
        chosen.candidates = candidates
        if chosen.scheme != best_scheme:
            reason = ("species" if genus and
                      self._scheme_genus.get(chosen.scheme, "") == genus else "ST resolution")
            LOG.debug("preferred scheme %s over %s on %s", chosen.scheme, best_scheme, reason)
        return chosen

    def _evaluate_scheme(self, scheme: str, loci_hits: dict[str, LocusHit]) -> MlstCall:
        """Build the allele profile and ST for one candidate scheme."""
        scheme_loci = self.profiles.loci(scheme) or sorted(loci_hits)
        alleles: dict[str, str] = {}
        novel: list[str] = []
        for locus in scheme_loci:
            hit = loci_hits.get(locus)
            if hit is None:
                alleles[locus] = "-"
            elif hit.exact:
                alleles[locus] = hit.allele
            elif hit.full_length:
                alleles[locus] = f"~{hit.allele}"
                novel.append(locus)
            else:
                alleles[locus] = f"{hit.allele}?"
                novel.append(locus)

        profile_key = tuple(alleles[locus] for locus in scheme_loci)
        _loci, table, extra = self.profiles.load(scheme)
        st = table.get(profile_key, "-")
        note = ""
        if st == "-":
            st, resolved, why = self._resolve_profile(table, scheme_loci, alleles)
            if st != "-":
                note = why
                alleles = resolved
                profile_key = tuple(alleles[locus] for locus in scheme_loci)
        if st == "-" and not note:
            absent = [locus for locus in scheme_loci if alleles[locus] == "-"]
            inexact = [locus for locus in scheme_loci
                       if alleles[locus] not in ("-",) and not alleles[locus].isdigit()]
            if absent and inexact:
                note = (f"{', '.join(absent)} not found and {', '.join(inexact)} "
                        f"inexact; no ST assigned")
            elif absent:
                note = f"{', '.join(absent)} not found; no ST assigned"
            elif inexact:
                note = f"novel allele(s) at {', '.join(inexact)}; no ST assigned"
            else:
                note = "allele combination not in profile table (novel ST)"
        clonal_complex = extra.get(st, "")
        if clonal_complex:
            note = (note + "; " if note else "") + f"clonal complex {clonal_complex}"
        # Count loci actually recovered from the genome: an allele filled in from
        # the profile table is not evidence, and inflating this would feed both
        # the species confidence and the scheme tie-break.
        exact_loci = sum(1 for locus in scheme_loci
                         if (loci_hits.get(locus) is not None and loci_hits[locus].exact))
        mean_identity = (sum(h.identity for h in loci_hits.values()) / len(loci_hits)
                         if loci_hits else 0.0)
        return MlstCall(scheme=scheme, sequence_type=st, alleles=alleles,
                        similarity=mean_identity, loci_found=exact_loci,
                        loci_total=len(scheme_loci), novel_alleles=novel, note=note)

    @staticmethod
    def _resolve_profile(table: dict[tuple[str, ...], str], scheme_loci: list[str],
                         alleles: dict[str, str]) -> tuple[str, dict[str, str], str]:
        """Second chance at an ST when the direct profile lookup missed.

        PubMLST records a deleted locus as allele ``0`` - *pstS* in the
        *Enterococcus faecium* CC17 lineages is the well known case - so a locus
        Hydra could not find is retried as ``0``.

        A missing locus is treated as a wildcard and the ST is only assigned when
        exactly one profile matches. Anything else would invent a type: in
        *Haemophilus influenzae* seventeen STs share the same six alleles and
        differ only at *fucK*, so a genome whose *fucK* fell on a contig break
        would otherwise be typed as whichever of them happens to carry a zero.
        """
        missing = [locus for locus in scheme_loci if alleles[locus] == "-"]
        if not missing or len(missing) > 2:
            return "-", alleles, ""
        resolved = dict(alleles)
        for locus in missing:
            resolved[locus] = "0"
        st = table.get(tuple(resolved[locus] for locus in scheme_loci))
        if st is None:
            return "-", alleles, ""

        # The ST stands, but say how much of it rests on the absent locus: the
        # same six alleles can fit several profiles that differ only there, and
        # the caller cannot tell a real deletion from a locus lost to a contig
        # break.
        known = [(index, alleles[locus]) for index, locus in enumerate(scheme_loci)
                 if locus not in missing]
        rivals = sorted({other for profile, other in table.items()
                         if other != st and all(profile[i] == v for i, v in known)})
        note = f"{', '.join(missing)} not found, typed as allele 0"
        if rivals:
            shown = ", ".join(rivals[:5]) + ("..." if len(rivals) > 5 else "")
            note += (f"; the alleles found also fit ST {shown}, which differ only at "
                     f"{', '.join(missing)}")
        return st, resolved, note

    # ----------------------------------------------------------------- reads
    def alleles_of(self, scheme: str, locus: str) -> dict[str, str]:
        """Every allele sequence of one locus, keyed by allele id."""
        cached = self._allele_cache.get((scheme, locus))
        if cached is not None:
            return cached
        path = self.root / "pubmlst" / scheme / f"{locus}.tfa"
        alleles: dict[str, str] = {}
        if path.exists():
            for header, sequence in read_fasta(path):
                name = header.split()[0]
                allele = name[len(locus):].lstrip("_-") if name.startswith(locus) else name
                alleles[allele or name] = sequence.upper()
        self._allele_cache[(scheme, locus)] = alleles
        return alleles

    def representative_reference(self, schemes: list[str], out_path: Path,
                                 share_groups: dict[str, str] | None = None) -> dict[str, list]:
        """Write one allele per locus for *schemes*, to map reads against.

        Mapping reads to every allele of a scheme is self-defeating: the alleles
        of a locus differ by a handful of bases, so each read maps equally well
        to hundreds of them and the depth splits. One representative per locus
        gives every read a single home, and the sample's actual allele is then
        read off the consensus.
        """
        records: list[tuple[str, str]] = []
        index: dict[str, list[tuple[str, str]]] = {}
        by_locus: dict[str, str] = {}
        for scheme in schemes:
            scheme_dir = self.root / "pubmlst" / scheme
            if not scheme_dir.is_dir():
                continue
            for locus in sorted(p.stem for p in scheme_dir.glob("*.tfa")):
                alleles = self.alleles_of(scheme, locus)
                if not alleles:
                    continue
                # Schemes for the same species share one representative per
                # locus: two Escherichia schemes both define adk, and giving
                # each its own near-identical reference would split the reads
                # between them so neither locus is fully covered. Schemes for
                # *different* species keep their own, since a sister species'
                # allele can be too divergent for the reads to map to.
                group = (share_groups or {}).get(scheme, scheme)
                key = f"{group}|{locus}"
                name = by_locus.get(key)
                if name is None:
                    # The lowest-numbered allele is the scheme's original reference.
                    chosen = min(alleles,
                                 key=lambda a: (not a.isdigit(), int(a) if a.isdigit() else 0, a))
                    name = f"{scheme}.{locus}"
                    by_locus[key] = name
                    records.append((name, alleles[chosen]))
                index.setdefault(name, []).append((scheme, locus))
        if not records:
            raise HydraError("no MLST alleles available for the requested schemes")
        write_fasta(out_path, records, wrap=0)
        return index

    def type_from_reads(self, consensus_by_reference: dict, index: dict[str, tuple],
                        genus: str = "", force_scheme: str | None = None,
                        min_breadth: float = 90.0) -> MlstCall:
        """Assign an ST from read consensus sequences, one per locus.

        The consensus is matched against every allele of its locus: an exact
        match gives the allele number, otherwise the closest allele is reported
        as novel. The call is marked as read-derived, because it rests on a
        mapped consensus rather than on assembled sequence.
        """
        by_scheme: dict[str, dict[str, LocusHit]] = defaultdict(dict)
        caveats: dict[str, list[str]] = {}
        for reference, consensus in consensus_by_reference.items():
            entries = index.get(reference)
            if not entries or consensus.breadth < min_breadth:
                continue
            called = consensus.sequence
            # Positions the reads did not cover deeply enough carry no evidence
            # either way. Counting them as mismatches would make almost every
            # locus inexact, because depth always tapers at a reference's ends.
            covered = [i for i, base in enumerate(called) if base != "-"]
            if not covered:
                continue
            complete = len(covered) >= (1.0 - MAX_UNCALLED_FRACTION) * len(called)
            # A consensus is always the representative's length, so an allele
            # that differs by an indel cannot be recovered from it and would be
            # silently rounded to the nearest same-length allele - a confident
            # wrong ST. The same applies to a mixed sample, whose consensus is a
            # chimera. In both cases the locus is reported as inexact instead.
            if consensus.has_indel_evidence:
                complete = False
                caveats.setdefault("indel", []).append(reference)
            if len(consensus.mixed_sites) >= MIXED_SITES_PER_LOCUS:
                complete = False
                caveats.setdefault("mixed", []).append(reference)
            for scheme, locus in entries:
                alleles = self.alleles_of(scheme, locus)
                best_allele = ""
                best_matches = -1
                for allele, sequence in alleles.items():
                    if len(sequence) != len(called):
                        continue
                    matches = sum(1 for i in covered if sequence[i] == called[i])
                    if matches > best_matches or (matches == best_matches and
                                                  _allele_rank(allele) < _allele_rank(best_allele)):
                        best_matches, best_allele = matches, allele
                if not best_allele:
                    continue
                identity = 100.0 * best_matches / len(covered)
                by_scheme[scheme][locus] = LocusHit(
                    locus=locus, allele=best_allele, identity=identity,
                    coverage=100.0 if complete else consensus.breadth,
                    bitscore=consensus.depth,
                )
        call = self._call_sample(dict(by_scheme), force_scheme, genus)
        call.source = "reads"
        notes = [call.note, "typed from mapped reads"]
        if caveats.get("indel"):
            loci = ", ".join(sorted({r.split(".")[-1] for r in caveats["indel"]}))
            notes.append(f"{loci} carry insertion/deletion evidence, so their alleles could "
                         f"not be read off a fixed-length consensus and are not exact")
        if caveats.get("mixed"):
            loci = ", ".join(sorted({r.split(".")[-1] for r in caveats["mixed"]}))
            notes.append(f"{loci} show intermediate allele fractions; this looks like a mixed "
                         f"or contaminated sample and its consensus may be a chimera")
        call.note = "; ".join(part for part in notes if part)
        return call

    def available_schemes(self) -> list[str]:
        pubmlst = self.root / "pubmlst"
        if not pubmlst.is_dir():
            return []
        return sorted(p.name for p in pubmlst.iterdir() if p.is_dir())


def _allele_rank(allele: str) -> tuple:
    """Numeric allele ids sort numerically and before any non-numeric id."""
    if not allele:
        return (2, 0, "")
    return (0, int(allele), "") if allele.isdigit() else (1, 0, allele)


def _reverse_name(name: str) -> tuple:
    """Sort key that puts the alphabetically first scheme on top of a descending sort."""
    return tuple(-ord(ch) for ch in name)


def _better(candidate: LocusHit, current: LocusHit) -> bool:
    """Prefer exact matches, then higher coverage and identity.

    Alleles that tie on all of those are separated by allele id, lowest first.
    PubMLST assigns ids in order of discovery, so the lower id is the canonical
    one, and nested alleles of different lengths (common in *ddl* of
    *Enterococcus faecium*) would otherwise be decided by the longer sequence's
    higher bit score and never match a profile.
    """
    return (candidate.exact, candidate.coverage, candidate.identity,
            -candidate.allele_number) > \
           (current.exact, current.coverage, current.identity, -current.allele_number)
