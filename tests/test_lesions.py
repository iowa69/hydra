"""What exactly mutated: hydra_amr.engines.lesions on synthetic genes (no BLAST).

Each test builds a reference protein and a coding sequence for it, edits the
coding sequence the way a real lesion would, embeds it in a contig on either
strand, and writes the BLAST alignment (qseq/sseq, coordinates) that blastx or
blastn would report. The description must then name the change and point at the
right bases of the contig, which the tests check by reading the contig directly.
"""

from __future__ import annotations

import random

from hydra_amr.engines.blast import Hsp
from hydra_amr.engines.lesions import (CODE, describe_disruption, describe_dna,
                                       describe_protein, infer_reference_codon)
from hydra_amr.engines.mutations import MutationEntry
from hydra_amr.engines.nucl import QueryPiece
from hydra_amr.records import Hit
from hydra_amr.utils import revcomp

RNG = random.Random(7)
#: One fixed codon per residue, avoiding homopolymer-prone choices where possible.
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTC", "G": "GGT", "H": "CAC",
         "I": "ATC", "K": "AAG", "L": "CTG", "M": "ATG", "N": "AAC", "P": "CCG", "Q": "CAG",
         "R": "CGT", "S": "TCG", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAC"}
PROTEIN = ("MSKLTEVRAYDHIGNQWFPSKLTEVRAYDHIGNQWFPCSKLTEVRAYDHIGNQWFPEVRAYDHIGNQWFPSKLTE"
           "VRAYDHIGNQWFPSKLTEVRAYDHIGNQWFPGDGDTYAKL")


def cds_for(protein: str) -> str:
    return "".join(CODON[a] for a in protein) + "TAA"


def flank(n: int) -> str:
    # Stop codons in all frames would make translation past a gene end trivial;
    # random sequence is the honest stand-in for intergenic DNA.
    return "".join(RNG.choice("ACGT") for _ in range(n))


def place(gene: str, strand: str, left: int = 300, right: int = 300):
    """Embed *gene* (5'->3') in a contig; return contig and the gene's first-base coordinate."""
    a, b = flank(left), flank(right)
    if strand == "+":
        return a + gene + b, left + 1
    return a + revcomp(gene) + b, left + len(gene)


def protein_hsp(isolate_cds: str, gene_top: int, strand: str, qseq: str, sseq: str,
                sstart: int, first_codon: int = 0, slen: int | None = None) -> Hsp:
    """A blastx HSP whose first aligned codon is codon *first_codon* of the isolate gene."""
    n_query = sum(c != "-" for c in qseq)
    if strand == "+":
        qstart = gene_top + 3 * first_codon
        qend = qstart + 3 * n_query - 1
    else:
        qstart = gene_top - 3 * first_codon
        qend = qstart - 3 * n_query + 1
    send = sstart + sum(c != "-" for c in sseq) - 1
    return Hsp(qseqid="q0_0", sseqid="ref", pident=99.0, length=len(qseq), mismatch=0,
               gapopen=0, qstart=qstart, qend=qend, sstart=sstart, send=send, evalue=0.0,
               bitscore=500.0, slen=slen or len(PROTEIN), qlen=10_000, gaps=0, nident=0,
               strand=strand, qseq=qseq, sseq=sseq)


def read(contig: str, lo: int, hi: int, strand: str) -> str:
    seq = contig[lo - 1:hi]
    return revcomp(seq) if strand == "-" else seq


def mutation_hit(gene: str, contig_name: str = "contig_1", **kw) -> Hit:
    return Hit(sample="s", database="protein", gene=gene, drug_class="QUINOLONE",
               subclass="QUINOLONE", sequence=contig_name, method="POINTX",
               resolution="POINT", **kw)


def entry(symbol: str, position: int, ref: str, alt: str, kind: str = "prot") -> MutationEntry:
    return MutationEntry(kind=kind, taxgroup="X", key="ref", position=position, symbol=symbol,
                         reported_symbol=symbol, ref=ref, alt=alt, drug_class="QUINOLONE",
                         subclass="QUINOLONE", name="", gene=symbol.split("_")[0])


# ------------------------------------------------------------ substitutions
def _substitution(strand: str, upstream_insertion: bool = False):
    position = 20                                    # PROTEIN[19] == "S"
    assert PROTEIN[position - 1] == "S"
    isolate = PROTEIN[:position - 1] + "L" + PROTEIN[position:]
    gene = list(cds_for(isolate))
    # Ser TCG -> Leu TTG: the second base changes.
    gene[3 * (position - 1):3 * position] = list("TTG")
    qseq, sseq = isolate, PROTEIN
    if upstream_insertion:
        # An extra codon (GCT, Ala) inserted after residue 5: the isolate is one
        # codon longer upstream of the mutation, the alignment carries a gap.
        gene[15:15] = list("GCT")
        qseq = isolate[:5] + "A" + isolate[5:]
        sseq = PROTEIN[:5] + "-" + PROTEIN[5:]
    gene = "".join(gene)
    contig, top = place(gene, strand)
    hsp = protein_hsp(gene, top, strand, qseq, sseq, sstart=1)
    hit = mutation_hit("gyrA", start=min(hsp.qstart, hsp.qend), end=max(hsp.qstart, hsp.qend),
                       strand=strand)
    evidence = {"kind": "protein", "hsp": hsp, "offset": 0,
                "entry": entry("gyrA_S20L", position, "S", "L"), "position": position}
    return describe_protein(hit, evidence, contig, None), contig


def _check_substitution(detail, contig):
    assert detail.symbol == "gyrA_S20L"
    assert detail.call_type == "substitution"
    assert detail.protein_change == "S20L" and detail.protein_hgvs == "p.Ser20Leu"
    assert detail.observed == "TTG"
    # The codon really is at the reported place on the contig, on the reported strand.
    assert read(contig, detail.start, detail.end, detail.strand) == "TTG"
    assert CODE["TTG"] == "L"
    assert detail.reference_inferred == "TCG"           # labelled inferred, not observed
    assert detail.reference == ""
    assert detail.nucleotide_change == "c.59C>T (inferred)"
    assert detail.base_change == "C>T at codon position 2"
    assert detail.nucleotide_basis.startswith("inferred")
    assert not detail.unresolved


def test_substitution_plus_strand():
    _check_substitution(*_substitution("+"))


def test_substitution_minus_strand():
    _check_substitution(*_substitution("-"))


def test_indel_upstream_of_the_mutation_is_followed():
    for strand in "+-":
        detail, contig = _substitution(strand, upstream_insertion=True)
        _check_substitution(detail, contig)


def test_position_in_an_alignment_gap_is_unresolved():
    position = 20
    gene = cds_for(PROTEIN[:position - 1] + PROTEIN[position:])        # residue 20 deleted
    contig, top = place(gene, "+")
    qseq = PROTEIN[:position - 1] + "-" + PROTEIN[position:]
    hsp = protein_hsp(gene, top, "+", qseq, PROTEIN, sstart=1)
    hit = mutation_hit("gyrA", start=hsp.qstart, end=hsp.qend, strand="+")
    detail = describe_protein(hit, {"kind": "protein", "hsp": hsp, "offset": 0,
                                    "entry": entry("gyrA_S20L", position, "S", "L"),
                                    "position": position}, contig, None)
    assert "gap" in detail.unresolved
    assert detail.observed == "" and detail.location == ""


def test_codon_past_the_contig_end_is_unresolved():
    gene = cds_for(PROTEIN)
    contig, top = place(gene, "+", right=0)
    contig = contig[:top - 1 + 3 * 19 + 1]               # cut inside codon 20
    hsp = protein_hsp(gene, top, "+", PROTEIN[:20], PROTEIN[:20], sstart=1)
    hit = mutation_hit("gyrA", start=hsp.qstart, end=hsp.qend, strand="+")
    detail = describe_protein(hit, {"kind": "protein", "hsp": hsp, "offset": 0,
                                    "entry": entry("gyrA_S20L", 20, "S", "L"), "position": 20},
                              contig, None)
    assert "end of the contig" in detail.unresolved


def test_parsimony_is_never_presented_as_fact():
    assert infer_reference_codon("TTG", "S") == ["TCG"]
    # Ile ATC from Ser: AGC (one step) beats TCN (two steps).
    assert infer_reference_codon("ATC", "S") == ["AGC"]
    # Two equally close candidates: reported as ambiguous, no base change claimed.
    assert infer_reference_codon("ACT", "S") == ["AGT", "TCT"]


def test_in_frame_insertion_codons_are_read_from_the_assembly():
    # GD duplication after D79 in ...G D G D T Y...: the catalogue names it D81DGD
    # after shifting right; the aligner put the gap on the left.
    base = PROTEIN.index("GDGDTY")
    isolate = PROTEIN[:base + 2] + "GD" + PROTEIN[base + 2:]
    gene = cds_for(isolate)
    contig, top = place(gene, "-")
    qseq = isolate
    sseq = PROTEIN[:base + 2] + "--" + PROTEIN[base + 2:]
    hsp = protein_hsp(gene, top, "-", qseq, sseq, sstart=1)
    anchor = base + 4                                    # D at 1-based base+4, 3'-most
    hit = mutation_hit("ompK36", start=min(hsp.qstart, hsp.qend),
                       end=max(hsp.qstart, hsp.qend), strand="-")
    e = entry(f"ompK36_D{anchor}DGD", anchor, "D", "DGD")
    detail = describe_protein(hit, {"kind": "protein", "hsp": hsp, "offset": 0, "entry": e,
                                    "position": anchor}, contig, "Klebsiella_pneumoniae")
    assert detail.call_type == "insertion"
    assert detail.observed == "GGT GAT"
    assert read(contig, detail.start, detail.end, "-") == "GGTGAT"
    assert detail.nucleotide_change == f"c.{3 * anchor}_{3 * anchor + 1}insGGTGAT"
    assert detail.kleborate.startswith("OmpK36:p.") and detail.kleborate.endswith("insGlyAsp")
    assert not detail.unresolved


# ------------------------------------------------------- nucleotide mutations
def test_nucleotide_catalogue_mutation_both_strands():
    reference = flank(400)
    position = 201
    ref_base = reference[position - 1]
    alt = {"A": "G", "C": "T", "G": "A", "T": "C"}[ref_base]
    isolate = reference[:position - 1] + alt + reference[position:]
    for strand in "+-":
        left = flank(150)
        if strand == "+":
            contig = left + isolate + flank(150)
            qseq, sseq = isolate, reference
            hsp = Hsp("q", "NC_1@23S_ribosomal_RNA@23S:1-400", 99.0, 400, 1, 0, 151, 550, 1, 400,
                      0.0, 700.0, 400, len(contig), 0, 399, "+", qseq, sseq)
            expected = 150 + position
        else:
            contig = left + revcomp(isolate) + flank(150)
            qseq, sseq = revcomp(isolate), revcomp(reference)
            hsp = Hsp("q", "NC_1@23S_ribosomal_RNA@23S:1-400", 99.0, 400, 1, 0, 151, 550, 400, 1,
                      0.0, 700.0, 400, len(contig), 0, 399, "-", qseq, sseq)
            expected = 150 + (400 - position + 1)
        hit = Hit(sample="s", database="protein", gene="23S", drug_class="MACROLIDE",
                  subclass="AZITHROMYCIN", sequence="c1", start=151, end=550, strand=strand,
                  method="POINTN", resolution="POINT")
        e = MutationEntry("dna", "X", "ref", position, f"23S_{ref_base}{position}{alt}", "",
                          ref_base, alt, "MACROLIDE", "AZITHROMYCIN", "", "23S")
        detail = describe_dna(hit, {"kind": "dna", "hsp": hsp, "offset": 0, "entry": e,
                                    "position": position, "reference": reference}, contig)
        assert detail.call_type == "rRNA"
        assert detail.reference == ref_base                    # from the reference sequence
        assert detail.reference_inferred == ""
        assert detail.observed == alt
        assert detail.nucleotide_change == f"n.{position}{ref_base}>{alt}"
        assert (detail.start, detail.strand) == (expected, strand)
        base = contig[expected - 1]
        assert (revcomp(base) if strand == "-" else base) == alt
        assert "reference sequence" in detail.nucleotide_basis
        assert detail.resistance == "azithromycin resistance (macrolide class)"
        assert not detail.unresolved


# ----------------------------------------------------------------- disruptions
LOF = "MKAILVTSGLNQEWRYDPHFAVLKTSGENQIRWYDLHPAVKTLSGENQIWRYDFHPALKTVSGEQNIRWYDLHPA" \
      "KVLTSGNEQIRWYDHLPAKTVLSGENQIRWYDLHPAKTVSLGENQIRWYDHLPAKTVSGLENQIRWYDHLPAKV"


def lof_hit(contig_name: str, lo: int, hi: int, strand: str) -> Hit:
    return Hit(sample="s", database="protein", gene="ompK35", drug_class="CARBAPENEM",
               subclass="CARBAPENEM", sequence=contig_name, start=lo, end=hi, strand=strand,
               element_subtype="DISRUPTION", method="DISRUPTX", resolution="PARTIAL")


def lof_evidence(hsps, best, reference=LOF, contig_name="c1"):
    piece = QueryPiece("s", contig_name, 0)
    return {"kind": "disruption", "hsp": best, "offset": 0, "contig": contig_name,
            "reason": "truncated", "reference": reference,
            "pieces": [(h, piece) for h in hsps]}


def lof_hsp(top, strand, sub_lo, sub_hi, isolate_protein, first_codon, slen=len(LOF)):
    q = isolate_protein[first_codon:first_codon + (sub_hi - sub_lo + 1)]
    return protein_hsp("", top, strand, q, LOF[sub_lo - 1:sub_hi], sstart=sub_lo,
                       first_codon=first_codon, slen=slen)


def test_nonsense_is_named_with_its_codon():
    for strand in "+-":
        gene = list(cds_for(LOF))
        gene[3 * 11:3 * 12] = list("TAG")                   # codon 12 Q -> stop (CAG -> TAG)
        assert LOF[11] == "Q"
        gene = "".join(gene)
        contig, top = place(gene, strand)
        isolate = LOF[:11] + "*" + LOF[12:]
        hsp = lof_hsp(top, strand, 1, len(LOF), isolate, 0)
        lo, hi = sorted((hsp.qstart, hsp.qend))
        detail = describe_disruption(lof_hit("c1", lo, hi, strand),
                                     lof_evidence([hsp], hsp), {"c1": contig}, None)
        assert detail.call_type == "disruption_stop"
        assert detail.symbol == "ompK35_Q12Ter" and detail.protein_hgvs == "p.Gln12Ter"
        assert detail.observed == "TAG"
        assert read(contig, detail.start, detail.end, strand) == "TAG"
        assert detail.reference_inferred == "CAG"
        assert detail.nucleotide_change == "c.34C>T (inferred)"
        assert detail.intact.startswith("residues 1-11 of")


def _frameshift(strand: str, delete_at: int):
    gene = cds_for(LOF)
    isolate_cds = gene[:delete_at] + gene[delete_at + 1:]
    contig, top = place(isolate_cds, strand)
    # Upstream part in the original frame, downstream part one base earlier.
    codon = delete_at // 3                                  # 0-based codon of the deletion
    up = lof_hsp(top, strand, 1, codon, LOF, 0)
    down_first = codon + 1                                  # first intact codon downstream
    down_top = (top + 3 * down_first - 1) if strand == "+" else (top - 3 * down_first + 1)
    down = protein_hsp("", down_top, strand, LOF[down_first:], LOF[down_first:],
                       sstart=down_first + 1, slen=len(LOF))
    lo = min(up.qstart, up.qend, down.qstart, down.qend)
    hi = max(up.qstart, up.qend, down.qstart, down.qend)
    detail = describe_disruption(lof_hit("c1", lo, hi, strand),
                                 lof_evidence([up, down], down), {"c1": contig}, None)
    return detail, contig, isolate_cds


def test_frameshift_is_named_as_amrfinderplus_does():
    for strand in "+-":
        delete_at = 3 * 50 + 1                              # second base of codon 51
        detail, contig, isolate_cds = _frameshift(strand, delete_at)
        assert detail.call_type == "disruption_frameshift"
        assert detail.base_change == "net 1-bp deletion"
        # The protein change is what translating the isolate gene gives.
        from hydra_amr.utils import translate
        product = translate(isolate_cds, False)
        first = next(i for i, (a, b) in enumerate(zip(product, LOF)) if a != b)
        novel = product.index("*", first) - first
        assert detail.protein_change == f"{LOF[first]}{first + 1}{product[first]}fsTer{novel}"
        assert detail.symbol == f"ompK35_{detail.protein_change}"
        assert detail.protein_hgvs.endswith(f"fsTer{novel + 1}")
        assert detail.intact.startswith(f"residues 1-{first} of {len(LOF)}")
        # The deleted base is inferred from the reference protein and says so.
        assert detail.nucleotide_change.endswith("(inferred)") or \
            detail.nucleotide_change.endswith("(inferred, ambiguous)")
        assert "inferred" in detail.nucleotide_basis


def test_interruption_by_inserted_sequence():
    gene = cds_for(LOF)
    cut = 3 * 60
    insert = flank(800) + "TAA"                             # 803 bp, carries a stop
    isolate_cds = gene[:cut] + insert + gene[cut:]
    contig, top = place(isolate_cds, "+")
    up = lof_hsp(top, "+", 1, 60, LOF, 0)
    down_top = top + cut + len(insert)
    down = protein_hsp("", down_top, "+", LOF[60:], LOF[60:], sstart=61, slen=len(LOF))
    detail = describe_disruption(lof_hit("c1", up.qstart, down.qend, "+"),
                                 lof_evidence([up, down], up), {"c1": contig}, None)
    assert detail.call_type == "disruption_interrupted"
    assert detail.base_change == f"{len(insert)} bp inserted"
    assert "not identified" in detail.nucleotide_basis
    assert detail.start <= top + cut <= detail.end + 1


def test_gene_cut_by_a_contig_end_is_not_called_a_disruption():
    gene = cds_for(LOF)
    for strand in "+-":
        # Only the first 70 codons are on the contig; it ends there.
        if strand == "+":
            contig, top = place(gene[:210], "+", right=0)
        else:
            contig, top = place(gene[:210], "-", left=0)
        hsp = lof_hsp(top, strand, 1, 70, LOF, 0)
        lo, hi = sorted((hsp.qstart, hsp.qend))
        detail = describe_disruption(lof_hit("c1", lo, hi, strand),
                                     lof_evidence([hsp], hsp), {"c1": contig}, None)
        assert detail.call_type == "contig_end"
        assert "NOT evidence" in detail.type_label
        assert "not evidence" in detail.description
        assert "PARTIAL_CONTIG_END" in detail.description


def test_missing_start_codon_region_is_located_not_guessed():
    gene = cds_for(LOF)
    contig, top = place(flank(400) + gene[3 * 40:], "+", left=0)   # codons 1-40 replaced
    first_top = 400 + 1
    hsp = protein_hsp("", first_top, "+", LOF[40:], LOF[40:], sstart=41, slen=len(LOF))
    detail = describe_disruption(lof_hit("c1", hsp.qstart, hsp.qend, "+"),
                                 lof_evidence([hsp], hsp), {"c1": contig}, None)
    assert detail.call_type == "disruption_5prime"
    assert detail.symbol == f"ompK35_M1_{LOF[39]}40del"
    assert "not identified" in detail.unresolved


# ------------------------------------------------------------- reads and reports
def test_read_site_is_described_from_the_reference_and_the_counts():
    from hydra_amr.engines.lesions import describe_read_site
    from hydra_amr.engines.reads import SiteCall

    call = SiteCall(reference="NC_004431.1@23S_ribosomal_RNA@23S:237160-240064", position=2059,
                    symbol="23S_A2059G", gene="23S", ref_base="A", alt_base="G", depth=200,
                    ref_count=150, alt_count=50, allele_fraction=0.25, p_value=1e-30,
                    drug_class="MACROLIDE", subclass="AZITHROMYCIN", status="heteroresistant")
    detail = describe_read_site(call)
    assert detail.call_type == "rRNA" and detail.nucleotide_change == "n.2059A>G"
    assert detail.reference == "A" and detail.reference_inferred == ""
    assert detail.observed == "G in 50/200 reads (allele fraction 0.250)"
    assert "no contig" in detail.location
    assert "heteroresistant" in detail.description


def test_read_variant_base_change_comes_from_the_reference_gene():
    from hydra_amr.engines.lesions import describe_read_variant
    from hydra_amr.engines.reads import Variant

    variant = Variant(reference="gyrA_ref", position=248, ref_base="C", alt_base="T", depth=80,
                      alt_count=80, allele_fraction=1.0, codon_position=83, ref_aa="S",
                      alt_aa="L")
    entry = MutationEntry("prot", "Escherichia", "WP_1", 83, "gyrA_S83L", "gyrA_S83L", "S",
                          "L", "QUINOLONE", "QUINOLONE", "", "gyrA")
    detail = describe_read_variant(variant, entry, "gyrA", "gyrA_ref", "QUINOLONE", "QUINOLONE")
    assert detail.symbol == "gyrA_S83L" and detail.protein_hgvs == "p.Ser83Leu"
    assert detail.nucleotide_change == "c.248C>T" and detail.reference == "C"
    assert detail.base_change == "C>T at codon position 2"
    assert detail.resistance == "quinolone resistance"


def test_mutation_table_and_long_table_carry_the_detail():
    from hydra_amr.records import HIT_COLUMNS, MutationDetail, SampleResult
    from hydra_amr.report.tables import long_table, mutation_table

    hit = Hit(sample="s1", database="protein", gene="gyrA", drug_class="QUINOLONE",
              sequence="c1", start=10, end=2600, method="POINTX", resolution="POINT",
              mutation=MutationDetail(symbol="gyrA_S83I", type_label="Amino-acid substitution",
                                      protein_change="S83I", protein_hgvs="p.Ser83Ile",
                                      observed="ATC", reference_inferred="AGC",
                                      location="c1:2350-2352 (-)",
                                      resistance="quinolone resistance"))
    other = Hit(sample="s1", database="ncbi", gene="blaKPC-2", sequence="c2", start=1, end=882)
    result = SampleResult(sample="s1", hits=[hit, other])
    table = mutation_table([result])
    assert list(table["mutation"]) == ["gyrA_S83I"]
    assert table.loc[0, "observed_codon"] == "ATC"
    assert table.loc[0, "reference_codon_inferred"] == "AGC"
    assert table.loc[0, "mutation_location"] == "c1:2350-2352 (-)"
    long = long_table([result])
    assert list(long.columns) == list(HIT_COLUMNS)
    assert set(long.loc[long["gene"] == "blaKPC-2", "mutation"]) == {""}
    assert hit.as_row()["protein_change_hgvs"] == "p.Ser83Ile"
