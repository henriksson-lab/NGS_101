"""Source-exact stLFR oligos from Cheng et al. Protocol Exchange supplement.

The 2019 defining paper points to doi:10.1038/protex.2018.116 for the bench
protocol.  Its public PDF links an XLSX workbook containing these common oligos
and all 1,536 members of each of the three barcode sets.  The catalogue models
the barcode architecture here and treats that workbook, rather than a copied
third-party whitelist, as the authoritative member list.
"""
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Segment, feature, oligo

BEAD_COMMON_T = "AAAAAAAAAATGTGAGCCAAGGAGTTG"
BEAD_COMMON_B = "CCAGAGCAACTCCTTGGCTCACA"
BRIDGE_RNA = "GCACUGACGACAUGAUCACCAAGGAUCGCCAUAGUCCAUGCUA"

TRANSPOSON_1_T = "CGATCCTTGGTGATCATGTCGTCAGTGCTTGTCTTCCTAAGATGTGTATAAGAGACAG"
TRANSPOSON_2_T = "GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG"
TRANSPOSON_B = "CTGUCTCUTATACACAUCT"
PCR1 = "TGTGAGCCAAGGAGTTG"
PCR2 = "GCCTCCCTCGCGCCATCAG"

BGI_R1 = "GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG"
BGI_BARCODE = "CGAGAACGTCTTGTGAGCCAAGGAGTTGCTCTGG"
BGI_R2 = "CGTCAGTGCTTGTCTTCCTAAGATGTGTATAAGAGACAG"
BGI_MDA1 = "TGATCACCAAGGATCGCCATAGTCCATGCTA"
BGI_MDA2 = "CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGA"

BRANCH_T = "CTGATGGCGCGAGGGAGGC"
BRANCH_B = "TCGCGCCATCA"
R1_GAP = "CAACTCCTTGGCTCACACGGAGGGAGCGCGGTAGTC"

ILLUMINA_TN1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ILLUMINA_TN2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
ILLUMINA_BRIDGE = "CCGAGCCCACGAGACCCAUAGUCCAUGCUA"
ILLUMINA_P5_GAP_F = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
ILLUMINA_P5_GAP_R = "AGATGTGTATAAGAGACA"
ILLUMINA_PCR1 = "CAAGCAGAAGACGGCATACGAGAT"
ILLUMINA_PCR2 = "AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC"
ILLUMINA_R1 = ILLUMINA_TN1
ILLUMINA_I7_STLFR = "CCGAGCCCACGAGACCCATAGTCCATGCTA"
ILLUMINA_R2 = ILLUMINA_TN2


def barcode_transfer_geometry():
    """Source-known 3′ end of the assembled bead oligo at the splinted nick."""
    return reported_adapter_end(
        "assembled bead barcode adapter", kind="splinted", polarity="3-prime",
        terminal_sequence="G", body_sequence="TAGCATGGACTATG",
        statement=("The workbook prints the bead-adapter 3′ flank and RNA bridge "
                   "spanning it and Transposon1T. This declares the reported adapter "
                   "end, not proof of the complete ligation. The figure's different "
                   "common-primer region is disclosed separately."))


def branch_adapter_geometry():
    return reported_adapter_end(
        "BranchT 5-prime end", kind="blunt", polarity="5-prime",
        terminal_sequence=BRANCH_T[0], body_sequence=BRANCH_T[1:],
        statement="The paper describes joining the blunt adapter 5′ phosphate to a recessed genomic 3′ hydroxyl. BranchT is phosphorylated and BranchB terminates in ddG in the source workbook. No complete recessed-substrate junction is constructed here.")


def barcode_read():
    """Source 42-position barcode string, without an assigned bead identity."""
    return Construct([_barcode("barcode 1", "N" * 10), Segment("BC1-BC2 linker", "TCTGCG"),
                      _barcode("barcode 2", "N" * 10), Segment("BC2-BC3 linker", "CCTTCC"),
                      _barcode("barcode 3", "N" * 10)], name="Three-set stLFR barcode string")


def assembled_bead_oligo():
    return Construct([Segment("bead common; workbook variant", BEAD_COMMON_T),
                      Segment("BC1 left flank", "CTCTGG"), *barcode_read().segments,
                      Segment("capture flank", "TAGCATGGACTATGG")],
                     name="Assembled bead oligo; ordered-workbook common sequence")


def _barcode(name: str, bases: str) -> Segment:
    return Segment(name, bases, "cbc", placeholder=True, length_bp=10,
                   feature=feature(name.lower().replace(" ", "_"),
                                   "linked_read_barcode", "whitelist",
                                   whitelist="stLFR Protocol Exchange XLSX set"))


def oligos():
    return (
        oligo("BeadCommonT", [Segment("poly(A) linker", "A" * 10),
                              Segment("bead common", BEAD_COMMON_T[10:])],
              five="5′-/52-Bio/-", three="-3′"),
        oligo("BeadCommonB", [Segment("bead common complement", BEAD_COMMON_B)]),
        oligo("Bridge (RNA)", [Segment("capture bridge", BRIDGE_RNA)]),
        oligo("Barcode1T family", [Segment("BC1 left flank", "CTCTGG"),
                                    _barcode("barcode 1", "N" * 10),
                                    Segment("BC1 right flank", "TCTGCG")],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Barcode2T family", [_barcode("barcode 2", "N" * 10)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Barcode3T family", [Segment("BC3 left flank", "CCTTCC"),
                                    _barcode("barcode 3", "N" * 10),
                                    Segment("BC3 right flank", "TAGCATGGACTATGG")],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Transposon1T", [Segment("capture transposon top", TRANSPOSON_1_T)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Transposon2T", [Segment("second transposon top", TRANSPOSON_2_T)]),
        oligo("TransposonB", [Segment("uracil transposon bottom", TRANSPOSON_B)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("PCR1", [Segment("bead-side PCR primer", PCR1)]),
        oligo("PCR2", [Segment("branch/transposon-side PCR primer", PCR2)]),
        oligo("BGI R1 primer (two-transposon table)", [Segment("BGI Read 1 primer", BGI_R1)]),
        oligo("BGI stLFR barcode sequencing primer",
              [Segment("BGI barcode primer", BGI_BARCODE)]),
        oligo("BGI R2 sequencing primer", [Segment("BGI Read 2 primer", BGI_R2)]),
        oligo("BGI MDA primer 1", [Segment("BGI circle-amplification primer 1", BGI_MDA1)]),
        oligo("BGI MDA primer 2", [Segment("BGI circle-amplification primer 2", BGI_MDA2)]),
        oligo("BranchT", [Segment("3-prime branch adapter top", BRANCH_T)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("BranchB", [Segment("3-prime branch adapter bottom", BRANCH_B)],
              five="5′-", three="-/3′dd/G-3′"),
        oligo("R1 Seq Primer Gap (3-prime branch table)", [Segment("gap-ligation Read 1 primer", R1_GAP)]),
        oligo("Illumina Nextera Tn1", [Segment("Illumina Tn1", ILLUMINA_TN1)]),
        oligo("Illumina Nextera Tn2", [Segment("Illumina Tn2", ILLUMINA_TN2)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Illumina bridge", [Segment("Illumina capture bridge", ILLUMINA_BRIDGE)]),
        oligo("Illumina P5 gap adapter F",
              [Segment("P5 gap adapter forward", ILLUMINA_P5_GAP_F)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Illumina P5 gap adapter R",
              [Segment("P5 gap adapter reverse", ILLUMINA_P5_GAP_R)],
              five="5′-", three="-/3′dd/G-3′"),
        oligo("Illumina PCR1", [Segment("Illumina P7 PCR primer", ILLUMINA_PCR1)]),
        oligo("Illumina PCR2 family", [Segment("source P5/i5 family", ILLUMINA_PCR2, placeholder=True)]),
        oligo("Illumina R1 sequencing primer", [Segment("Illumina Read 1 primer", ILLUMINA_R1)]),
        oligo("Illumina stLFR i7 sequencing primer",
              [Segment("Illumina stLFR index primer", ILLUMINA_I7_STLFR)]),
        oligo("Illumina R2 sequencing primer", [Segment("Illumina Read 2 primer", ILLUMINA_R2)]),
    )
