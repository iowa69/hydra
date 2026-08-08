"""Whole-element typing: SCCmec, and the capsule loci that work the same way.

Some typing schemes are not allele profiles over housekeeping loci. They ask which
of a set of reference *elements* a genome carries -- the SCCmec cassette in
staphylococci, the capsule locus in pneumococcus, meningococcus and *Haemophilus*.
The element is tens of kilobases, it is either present or it is not, and the answer
is the reference it covers.

That makes the measurement coverage of the reference, not identity to an allele.
On the validation panel a genome carrying SCCmec II covers the type II reference
100% and a methicillin-susceptible genome covers the best-matching reference 21% --
the shared *orfX* flank that every *S. aureus* carries. The gap is wide, and the
threshold sits in the middle of it rather than at the edge of either side.

Every scheme is scoped to the genera it was defined in. Reporting an SCCmec type
for something that is not a staphylococcus would be reading the flanking region and
calling it a cassette.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..config import Config
from ..db.manager import DatabaseStore
from ..records import SpeciesCall, TypingResult
from ..utils import LOG
from ..engines.blast import blast, check_db_exists
from ..engines.nucl import QueryBatch


@dataclass(frozen=True)
class ElementScheme:
    """One whole-element typing scheme."""

    name: str
    #: Installed database holding the reference elements.
    database: str
    #: Genera the scheme is defined for. Empty means any, which no scheme uses.
    genera: frozenset
    #: Fraction of the reference that must be covered before a type is called.
    min_coverage: float
    #: Identity floor for the alignments that count toward coverage.
    min_identity: float
    #: What a genome that carries no element is reported as.
    absent: str = "-"
    note: str = ""


#: `>SCCmec_type_II(2A)|gb|D86934.2` and
#: `>SCCmec_type_IV(2B)|SCCmec_type_IVa(2B)|gb|AB063172.2` -- the second form names
#: the type and then its subtype, so both are kept and reported together.
SCHEMES: dict[str, ElementScheme] = {
    "SCCmec": ElementScheme(
        name="SCCmec", database="sccmec",
        genera=frozenset({"Staphylococcus"}),
        # Midway between the 97.6% a real cassette covers and the 21% of shared
        # flanking sequence a methicillin-susceptible genome covers.
        min_coverage=60.0, min_identity=90.0,
        absent="none",
        note="mecA-carrying cassette; 'none' means no cassette reached the coverage floor",
    ),
}


def _coverage(intervals: list[tuple[int, int]], length: int) -> float:
    """Fraction of a reference covered by a set of alignments, merged."""
    if not length:
        return 0.0
    covered = end = 0
    for start, stop in sorted(intervals):
        if start > end:
            covered += stop - start + 1
            end = stop
        elif stop > end:
            covered += stop - end
            end = stop
    return 100.0 * covered / length


class CassetteTyper:
    """Types every applicable whole-element scheme in one BLAST pass per database."""

    def __init__(self, store: DatabaseStore, config: Config):
        self.store = store
        self.config = config

    def applicable(self, species: SpeciesCall | None) -> list[ElementScheme]:
        """Schemes worth running for a species call.

        A scheme whose database is not installed is skipped rather than warned
        about on every sample, and a sample with no genus is typed by none of
        them: the element is only interpretable in the organisms it was defined in.
        """
        genus = (species.genus if species else "") or ""
        return [s for s in SCHEMES.values()
                if genus in s.genera and self.store.is_installed(s.database)]

    def type_batch(self, batch: QueryBatch, workdir: Path, threads: int = 1,
                   species_by_sample: dict[str, SpeciesCall] | None = None,
                   ) -> dict[str, list[TypingResult]]:
        species_by_sample = species_by_sample or {}
        wanted: dict[str, set] = {}
        for sample, species in species_by_sample.items():
            for scheme in self.applicable(species):
                wanted.setdefault(scheme.name, set()).add(sample)
        if not wanted:
            return {}

        out: dict[str, list[TypingResult]] = {}
        for scheme_name, samples in wanted.items():
            scheme = SCHEMES[scheme_name]
            handle = self.store.handle(scheme.database)
            check_db_exists(handle.fasta, "nucl")
            LOG.info("%s typing: %d samples against %d reference elements",
                     scheme.name, len(samples), handle.n_sequences)
            hsps = blast(
                "blastn", batch.path, handle.fasta,
                workdir / f"{scheme.database}.blastn.tsv",
                threads=threads, evalue=1e-20, task="megablast",
                perc_identity=scheme.min_identity, max_target_seqs=10000,
            )
            # reference -> intervals, per sample
            per_sample: dict[str, dict[str, list]] = {}
            lengths: dict[str, int] = {}
            meta: dict[str, dict] = {}
            for hsp in hsps:
                piece = batch.id_map.get(hsp.qseqid)
                if piece is None or piece.sample not in samples:
                    continue
                lengths[hsp.sseqid] = hsp.slen
                meta[hsp.sseqid] = handle.meta_for(hsp.sseqid) or {}
                per_sample.setdefault(piece.sample, {}).setdefault(hsp.sseqid, []).append(
                    (min(hsp.sstart, hsp.send), max(hsp.sstart, hsp.send)))

            for sample in samples:
                out.setdefault(sample, []).append(
                    self._call(scheme, per_sample.get(sample, {}), lengths, meta))
        return out

    def _call(self, scheme: ElementScheme, by_reference: dict[str, list],
              lengths: dict[str, int], meta: dict[str, dict]) -> TypingResult:
        ranked = sorted(
            ((_coverage(iv, lengths.get(ref, 0)), ref) for ref, iv in by_reference.items()),
            reverse=True)
        if not ranked or ranked[0][0] < scheme.min_coverage:
            best = f"{ranked[0][0]:.0f}% of the closest reference" if ranked else "no alignment"
            return TypingResult(scheme=scheme.name, call=scheme.absent, lineage="-",
                                note=f"below the {scheme.min_coverage:.0f}% coverage "
                                     f"floor ({best})")
        coverage, reference = ranked[0]
        entry = meta.get(reference) or {}
        # The importer puts the type in `gene` and the subtype, where the reference
        # is subtyped, in `product`.
        call = entry.get("gene") or "-"
        subtype = entry.get("product") or "-"
        accession = entry.get("accession") or reference
        if call == "-":
            return TypingResult(scheme=scheme.name, call="-", lineage="-",
                                note=f"reference {reference} carries no readable type")
        return TypingResult(
            scheme=scheme.name, call=call, lineage=subtype,
            loci_found=1, loci_total=1,
            note=f"{coverage:.0f}% of {accession} covered",
        )
