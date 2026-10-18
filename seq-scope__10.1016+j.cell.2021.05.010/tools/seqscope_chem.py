"""Oligos printed in the defining Seq-Scope paper (Cho et al., Cell 2021, PMC8238917).

The paper's seed backbone labels the 5' segment "P5 sequence" and the 3' segment
"P7 sequence" (JATS Methods line 217).  The printed bases are the Illumina P7 adapter
and the reverse complement of P5; the constants below are named for their bases and
the construction asserts that identity against ``illumina``.
"""
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, oligo, revcomp
from protocol_model import role


SEED_5P = "CAAGCAGAAGACGGCATACGAGAT"      # printed; paper label "P5 sequence"
TR1 = "TCTTTCCCTACACGACGCTCTTCCGATCT"
HDMI20 = "NNVNNVNNVNNVNNVNNNNN"
HDMI32 = "NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN"
HR1 = "TCTTGTGACTACAGCACCCTCGACTCTCGC"
POLYT = "T" * 27
DRAI = "TTTAAA"
DRAI_ADAPTER = "GACTTTCACCAGTCCATGAT"
SEED_3P = "GTGTAGATCTCGGTGGTCGCCGTATCATT"   # printed; paper label "P7 sequence"
if SEED_5P != il.P7 or SEED_3P != il.P5_RC:
    raise AssertionError("Seq-Scope seed end segments no longer equal P7 / revcomp(P5)")
READ1_DRAI = "ATCATGGACTGGTGAAAGTCTTTAAA" + "A" * 27 + "GCGAGAGTCGAGGGTGCTGTAGTCACAAGA"
TR2_RANDOM = "TCAGACGTGTGCTCTTCCGATCT" + "N" * 9

HDMI_ECORI = "HNNBNBNBNBNBNBNBNNNN"
HR1B = "CCCGTTCGCAACATGTCTGGCGTCATA"
ECORI = "GAATTC"
ECORI_ADAPTER = "CGCAGTCCAG"
READ1_ECORI = "CTGGACTGCGGAATTCTATGACGCCAGACATGTTGCGAACGGG"
UMI_OLIGO = "A" * 30 + "N" * 8 + "C" + "TATGACGCCAGACATGTTGCGAACGGG"


def _hdmi(identifier, bases):
    return Segment("HDMI spatial barcode", bases, "cbc", placeholder=True,
                   length_bp=len(bases),
                   feature=feature(identifier, "spatial_barcode", "random"))


def _seed(name, identifier, bases, hr):
    return oligo(name, [
        Segment("P7 adapter (paper label: P5)", SEED_5P), Segment("TruSeq Read 1", TR1),
        _hdmi(identifier, bases), Segment("HDMI read arm", hr),
        Segment("oligo-dT capture", POLYT), Segment("DraI site", DRAI),
        Segment("DraI adapter", DRAI_ADAPTER),
        Segment("P5 reverse complement (paper label: P7)", SEED_3P),
    ])


def oligos():
    sample_index = feature("seqscope_i7", "sample_index", "unknown")
    return (
        _seed("HDMI-DraI", "seqscope_hdmi20", HDMI20, HR1),
        _seed("HDMI32-DraI", "seqscope_hdmi32", HDMI32, HR1),
        oligo("Read1-DraI custom sequencing primer",
              [Segment("Read1-DraI", READ1_DRAI)]),
        oligo("TruSeq Read2-conjugated random primer", [
            Segment("TruSeq Read 2", TR2_RANDOM[:-9]),
            Segment("9-nt random UMI", TR2_RANDOM[-9:], "umi",
                    feature=feature("seqscope_random_primer_umi", "umi", "random")),
        ]),
        oligo("first-round PCR forward", [
            Segment("forward-primer body", "TCTTTCCCTACACGACGC"),
            Segment("phosphorothioate-protected 3-prime bases", "*T*C",
                    placeholder=True, length_bp=2),
        ]),
        oligo("first-round PCR reverse", [
            Segment("reverse-primer body", "TCAGACGTGTGCTCTTCC"),
            Segment("phosphorothioate-protected 3-prime bases", "*G*A",
                    placeholder=True, length_bp=2),
        ]),
        oligo("second-round PCR forward", [
            Segment("full forward-primer body", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCT"),
            Segment("phosphorothioate-protected 3-prime bases", "*T*C",
                    placeholder=True, length_bp=2),
        ]),
        oligo("second-round PCR reverse", [
            Segment("P7-side primer body", "CAAGCAGAAGACGGCATACGAGAT"),
            Segment("8-mer sample index", "N" * 8, "cbc", placeholder=True,
                    feature=sample_index),
            Segment("Read-2-side primer body", "GTGACTGGAGTTCAGACGTGTGCTCTTCC"),
            Segment("phosphorothioate-protected 3-prime bases", "*G*A",
                    placeholder=True, length_bp=2),
        ]),
        oligo("HDMI-TruEcoRI", [
            Segment("P7 adapter (paper label: P5)", SEED_5P), Segment("TruSeq Read 1", TR1),
            _hdmi("seqscope_hdmi_ecori", HDMI_ECORI),
            Segment("HDMI read arm B", HR1B), Segment("EcoRI site", ECORI),
            Segment("EcoRI adapter", ECORI_ADAPTER),
            Segment("P5 reverse complement (paper label: P7)", SEED_3P),
        ]),
        oligo("Read1-EcoRI custom sequencing primer",
              [Segment("Read1-EcoRI", READ1_ECORI)]),
        oligo("UMI-oligo for the EcoRI array", [
            Segment("oligo-dA", "A" * 30),
            Segment("8-nt UMI", "N" * 8, "umi",
                    feature=feature("seqscope_ecori_umi", "umi", "random")),
            Segment("junction C", "C"), Segment("HR1B-binding arm", UMI_OLIGO[-28:]),
        ]),
    )


# PCR primers as printed (Methods, Library Construction, JATS lines 279-281); the
# phosphorothioate marks are dropped only for site location.
PCR1_FORWARD = "TCTTTCCCTACACGACGCTC"
PCR1_REVERSE = "TCAGACGTGTGCTCTTCCGA"
PCR2_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTC"
PCR2_REVERSE_5P = "CAAGCAGAAGACGGCATACGAGAT"
PCR2_REVERSE_3P = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGA"


def _merge(left, right):
    """Join two printed sequences across their longest exact overlap."""
    for n in range(min(len(left), len(right)), 7, -1):
        if left.endswith(right[:n]):
            return left + right[n:]
    raise AssertionError("printed oligos do not overlap")


def _core(*, p7_graft=False):
    """Top-strand core after DraI-exposed oligo-dT capture and random-primed copy.

    The 3' T run is the printed 27 T plus the TTT half of the blunt DraI TTT|AAA cut;
    the random primer's N9 is copied as the bottom-strand 5' tract (read 2 cycles 1-9).
    """
    tr2_adapter = TR2_RANDOM[:-9]
    left = [Segment("P7 adapter carried from the seed 5' end", SEED_5P)] if p7_graft else []
    return left + [
        Segment("TruSeq Read 1", TR1),
        _hdmi("seqscope_hdmi20", HDMI20), Segment("HR1", HR1),
        Segment("oligo-dT (27 T + DraI TTT)", POLYT + DRAI[:3]),
        role("copied transcript cDNA upstream of N9 priming footprint; variable length"),
        Segment("random-primer N9 UMI", "N" * 9, "umi", placeholder=True,
                feature=feature("seqscope_random_primer_umi", "umi", "random")),
        Segment("TruSeq Read 2 adapter (random-primer tail)", revcomp(tr2_adapter)),
    ]


def graph_oligo(name, *parts):
    """Whole-strand graph state for a printed oligo (phosphorothioates named, not drawn)."""
    return Construct([Segment(name, "".join(parts))], name)


def graph_states():
    """Exact whole-strand states for the graph's oligo inputs (HDMI windows excluded)."""
    return {
        "read1_drai": graph_oligo("Read1-DraI", READ1_DRAI),
        "pcr1_forward": graph_oligo("first-round forward (PS before last two bases)", PCR1_FORWARD),
        "pcr1_reverse": graph_oligo("first-round reverse (PS before last two bases)", PCR1_REVERSE),
        "pcr2_forward": graph_oligo("second-round forward (PS before last two bases)", PCR2_FORWARD),
    }


def second_strand_duplex():
    """Surface cDNA strand paired with the eluted second strand (PCR1 template)."""
    return Construct(_core(p7_graft=True), "Seq-Scope second-strand duplex")


def pcr1_product():
    return Construct(_core(), "Seq-Scope first-round PCR product")


def final_library():
    """Second-round PCR product: P5/TR1 forward and P7/i7/TR2 reverse (printed)."""
    if not PCR2_FORWARD.startswith(il.P5) or PCR2_REVERSE_5P != il.P7:
        raise AssertionError("second-round PCR primer arms no longer equal P5/P7")
    tr2_arm = revcomp(_merge(PCR2_REVERSE_3P, TR2_RANDOM[:-9]))
    if tr2_arm != revcomp(il.TRUSEQ_READ2):
        raise AssertionError("printed reverse primers no longer rebuild TruSeq Read 2")
    core = _core()[:-1]
    lib = Construct([
        Segment("P5", il.P5, "p5"), *core,
        Segment("TruSeq Read 2 / Index 1 arm", tr2_arm, "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("seqscope_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], "Seq-Scope final library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers
