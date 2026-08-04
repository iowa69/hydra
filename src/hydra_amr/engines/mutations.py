"""Resistance-conferring point-mutation detection.

Two catalogues are supported, both taken from the AMRFinderPlus data directory:

* **protein mutations** (``AMRProt-mutation.tsv``) — checked against translated
  alignments of the assembly to the reference protein;
* **DNA mutations** (``AMR_DNA-<organism>.tsv``) — checked against nucleotide
  alignments to organism-specific reference loci such as 23S rRNA, *gyrA*
  promoters or *pbp4*.

Reference coordinates are 1-based; a negative DNA position counts back from the
end of the reference record, which is how promoter mutations such as
``pbp4_T-266A`` are catalogued.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from ..db.registry import protein_dir
from ..utils import HydraError


#: ``gyrA_S83L`` -> ref ``S``, position ``83``, alt ``L``; also ``pbp4_T-266A``,
#: ``rplD_WR65del`` and ``mgrB_Q30STOP``.
_SYMBOL_RE = re.compile(r"^(?P<ref>[A-Za-z*]*?)(?P<pos>-?\d+)(?P<alt>[A-Za-z*]*)$")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _complement(base: str) -> str:
    return base.translate(_COMPLEMENT)


@dataclass(frozen=True)
class MutationEntry:
    """One catalogued resistance mutation."""

    kind: str            # "prot" | "dna"
    taxgroup: str
    key: str             # protein accession, or DNA reference header
    position: int
    symbol: str
    reported_symbol: str
    ref: str
    alt: str
    drug_class: str
    subclass: str
    name: str            # e.g. Escherichia_quinolone_resistant_GyrA
    gene: str            # e.g. gyrA

    @property
    def is_deletion(self) -> bool:
        return self.alt.lower() in ("del", "deletion")

    @property
    def is_stop(self) -> bool:
        return self.alt == "*"

    @property
    def change_type(self) -> str:
        """How the alternate allele differs from the reference.

        ``ampC_G-15GG`` is an inserted G, not a G-to-G substitution: treating
        every entry as a single-base substitution would make the alternate
        allele identical to the reference and call it fixed in every wild-type
        sample.
        """
        if self.is_deletion:
            return "deletion"
        if self.is_stop:
            return "nonsense"
        if not self.ref or not self.alt:
            return "complex"
        if len(self.ref) == 1 and len(self.alt) == 1:
            return "substitution" if self.ref.upper() != self.alt.upper() else "complex"
        if len(self.alt) > len(self.ref) and self.alt.upper().startswith(self.ref.upper()):
            return "insertion"
        if len(self.alt) < len(self.ref) and self.ref.upper().startswith(self.alt.upper()):
            return "deletion"
        return "complex"

    @property
    def inserted_bases(self) -> str:
        """The bases an insertion adds after the reference position."""
        if self.change_type != "insertion":
            return ""
        return self.alt[len(self.ref):].upper()


def _parse_symbol(symbol: str, position: int) -> tuple[str, str, str]:
    """Split a mutation symbol into (gene, ref allele, alt allele)."""
    if "_" in symbol:
        gene, _, core = symbol.rpartition("_")
    else:
        gene, core = "", symbol
    match = _SYMBOL_RE.match(core)
    if not match:
        return gene or symbol, "", ""
    ref = match.group("ref")
    alt = match.group("alt")
    # The catalogue writes nonsense mutations as 'Ter'; older releases used
    # 'STOP'. Both mean the same stop codon.
    if alt.upper() in ("STOP", "TER", "*"):
        alt = "*"
    return gene or symbol, ref, alt


def _require_columns(path: Path, header: list[str], needed: tuple[str, ...]) -> dict[str, int]:
    """Column indexes, with a clear error when the upstream table changes shape."""
    idx = {name: i for i, name in enumerate(header)}
    missing = [name for name in needed if name not in idx]
    if missing:
        raise HydraError(
            f"{path.name} is missing the column(s) {', '.join(missing)}; the AMRFinderPlus "
            f"data format has changed. Re-run 'hydra db import --force protein'.")
    return idx


def _read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    if not path.exists():
        return [], []
    rows: list[list[str]] = []
    with open(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))
            rows.append(row)
    return header, rows


class MutationCatalog:
    """All mutation knowledge for one organism (or all organisms when unset)."""

    def __init__(self, db_root: Path, organism: str | None = None,
                 all_organisms: bool = False):
        self.db_root = Path(db_root)
        self.organism = organism
        #: Only an explicit request loads every taxgroup at once; see
        #: :meth:`_taxgroup_matches`.
        self.all_organisms = all_organisms
        self.protein: dict[str, list[MutationEntry]] = {}
        self.dna: dict[str, list[MutationEntry]] = {}
        self.suppress: set[str] = set()
        self.susceptible: dict[str, dict] = {}
        self.dna_reference: Path | None = None
        self._load()

    # --------------------------------------------------------------- loading
    def _taxgroup_matches(self, taxgroup: str) -> bool:
        """A catalogue row applies when its taxgroup is the organism or a parent of it.

        ``Escherichia`` covers ``Escherichia_coli``; ``Campylobacter`` covers all
        campylobacters.

        With no organism, nothing matches. Point mutations are only meaningful
        against the species they were catalogued for: applying every taxgroup at
        once to an unidentified genome reports Neisseria penA, Haemophilus ftsI
        and Staphylococcus fusA resistance in the same sample.
        """
        if self.organism is None:
            return self.all_organisms
        if not taxgroup:
            return False
        org = self.organism
        return org == taxgroup or org.startswith(taxgroup + "_") or taxgroup.startswith(org + "_")

    def _load(self) -> None:
        prot_dir = protein_dir(self.db_root)
        header, rows = _read_tsv(prot_dir / "AMRProt-mutation.tsv")
        if header:
            idx = _require_columns(prot_dir / "AMRProt-mutation.tsv", header,
                                   ("taxgroup", "accession_version", "mutation_position",
                                    "standard_mutation_symbol", "class", "subclass"))
            for row in rows:
                taxgroup = row[idx["taxgroup"]]
                if not self._taxgroup_matches(taxgroup):
                    continue
                try:
                    position = int(row[idx["mutation_position"]])
                except ValueError:
                    continue
                symbol = row[idx["standard_mutation_symbol"]]
                gene, ref, alt = _parse_symbol(symbol, position)
                if not alt:
                    continue
                accession = row[idx["accession_version"]]
                self.protein.setdefault(accession, []).append(MutationEntry(
                    kind="prot", taxgroup=taxgroup, key=accession, position=position,
                    symbol=symbol, reported_symbol=(row[idx["reported_mutation_symbol"]]
                                     if "reported_mutation_symbol" in idx else "") or symbol,
                    ref=ref, alt=alt, drug_class=row[idx["class"]], subclass=row[idx["subclass"]],
                    name=row[idx["mutated_protein_name"]] if "mutated_protein_name" in idx else "",
                    gene=gene,
                ))

        header, rows = _read_tsv(prot_dir / "AMRProt-suppress.tsv")
        if header:
            idx = _require_columns(prot_dir / "AMRProt-suppress.tsv", header,
                                   ("taxgroup", "protein_accession"))
            for row in rows:
                if self._taxgroup_matches(row[idx["taxgroup"]]):
                    self.suppress.add(row[idx["protein_accession"]])

        header, rows = _read_tsv(prot_dir / "AMRProt-susceptible.tsv")
        if header:
            idx = _require_columns(prot_dir / "AMRProt-susceptible.tsv", header,
                                   ("taxgroup", "accession_version", "resistance_cutoff",
                                    "gene_symbol", "class", "subclass"))
            for row in rows:
                if not self._taxgroup_matches(row[idx["taxgroup"]]):
                    continue
                try:
                    cutoff = float(row[idx["resistance_cutoff"]].strip())
                except ValueError:
                    continue
                self.susceptible[row[idx["accession_version"]]] = {
                    "gene": row[idx["gene_symbol"]], "cutoff": cutoff,
                    "class": row[idx["class"]], "subclass": row[idx["subclass"]],
                    "name": (row[idx["resistance_protein_name"]]
                             if "resistance_protein_name" in idx else ""),
                }

        if self.organism:
            mut_dir = self.db_root / "mutation" / "dna"
            reference = mut_dir / f"{self.organism}.fna"
            tsv = mut_dir / f"{self.organism}.tsv"
            if reference.exists():
                self.dna_reference = reference
            header, rows = _read_tsv(tsv)
            if header:
                for row in rows:
                    key = row[0]
                    try:
                        position = int(row[1])
                    except ValueError:
                        continue
                    symbol = row[2]
                    gene, ref, alt = _parse_symbol(symbol, position)
                    if not alt:
                        continue
                    self.dna.setdefault(key, []).append(MutationEntry(
                        kind="dna", taxgroup=self.organism, key=key, position=position,
                        symbol=symbol, reported_symbol=row[3] if len(row) > 3 else symbol,
                        ref=ref, alt=alt,
                        drug_class=row[4] if len(row) > 4 else "",
                        subclass=row[5] if len(row) > 5 else "",
                        name=row[6] if len(row) > 6 else "", gene=gene,
                    ))

    @property
    def n_protein(self) -> int:
        return sum(len(v) for v in self.protein.values())

    @property
    def n_dna(self) -> int:
        return sum(len(v) for v in self.dna.values())

    def describe(self) -> str:
        return (f"{self.n_protein} protein and {self.n_dna} DNA mutations"
                f"{f' for {self.organism}' if self.organism else ''}")


def reference_index(position: int, ref_length: int) -> int:
    """Convert a catalogue position to a 0-based index into the reference.

    Positive positions are 1-based from the start; negative positions count back
    from the end (promoter coordinates).
    """
    return position - 1 if position > 0 else ref_length + position


@dataclass
class AlignedObservation:
    """What the query carries at one reference position."""

    observed: str
    reference: str
    aligned: bool
    query_position: int | None = None


def walk_alignment(qseq: str, sseq: str, sstart: int, send: int,
                   targets: dict[int, None], *, nucleotide: bool,
                   qstart: int = 0, qend: int = 0) -> dict[int, AlignedObservation]:
    """Read the query residues aligned to a set of subject positions.

    ``qseq``/``sseq`` are BLAST's gapped alignment strings. For a minus-strand
    nucleotide alignment BLAST prints the subject in descending coordinates and
    complemented, so both strings are complemented back to subject orientation
    before comparison.
    """
    step = 1 if send >= sstart else -1
    minus = step == -1 and nucleotide
    out: dict[int, AlignedObservation] = {}
    subject_pos = sstart
    query_offset = 0
    for q_char, s_char in zip(qseq, sseq):
        if s_char == "-":
            if q_char != "-":
                query_offset += 1
            continue
        if subject_pos in targets:
            observed = q_char
            reference = s_char
            if minus:
                observed = _complement(observed) if observed != "-" else "-"
                reference = _complement(reference)
            out[subject_pos] = AlignedObservation(
                observed=observed.upper(), reference=reference.upper(),
                aligned=True, query_position=query_offset,
            )
        if q_char != "-":
            query_offset += 1
        subject_pos += step
    for position in targets:
        if position not in out:
            out[position] = AlignedObservation(observed="", reference="", aligned=False)
    return out


def multi_residue(qseq: str, sseq: str, sstart: int, send: int, position: int,
                  length: int, *, nucleotide: bool) -> AlignedObservation:
    """Observation spanning *length* reference residues starting at *position*."""
    targets = {position + offset: None for offset in range(length)}
    seen = walk_alignment(qseq, sseq, sstart, send, targets, nucleotide=nucleotide)
    observed = "".join(seen[position + offset].observed for offset in range(length))
    reference = "".join(seen[position + offset].reference for offset in range(length))
    aligned = all(seen[position + offset].aligned for offset in range(length))
    first = seen[position]
    return AlignedObservation(observed=observed, reference=reference, aligned=aligned,
                              query_position=first.query_position)


def evaluate(entry: MutationEntry, obs: AlignedObservation) -> tuple[bool, str]:
    """Decide whether *obs* satisfies *entry*. Returns (is_mutant, note)."""
    if not obs.aligned:
        return False, "position not covered by alignment"
    change = entry.change_type
    if change == "deletion":
        if obs.observed.replace("-", "") == "":
            return True, "deletion"
        return False, ""
    expected_ref = entry.ref.upper()
    if expected_ref and obs.reference and obs.reference != expected_ref:
        # The alignment did not land on the catalogued residue - do not call.
        return False, f"reference mismatch (expected {expected_ref}, aligned {obs.reference})"
    if change == "insertion":
        # An inserted base sits between two reference positions, which a
        # position-indexed alignment walk cannot see; calling it would need the
        # query gap structure. Never guess.
        return False, "insertion not callable from a position-indexed alignment"
    if change == "complex":
        return False, f"unsupported change '{entry.ref}->{entry.alt}'"
    if obs.observed.upper() == entry.alt.upper():
        return True, ""
    return False, ""
