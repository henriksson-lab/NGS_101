"""Published ChIL-seq oligos from Handa et al. Supplementary Table 4."""
from chemdraw import Segment, oligo

CHIL_V1_FW = ("CCCATGCATTAATACGACTCACTATAGGGTAGATCGCTCG"
              "TCGGCAGCGTCAGATGTGTATAAGAGACAG")
CHIL_V1_RV = ("CTGTCTCTTATACACATCTGACGCTGCCGACGAGCGATCT"
              "ACCCTATAGTGAGTCGTATTAATGCATGGG")
MEDSB_FW = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
MEDSB_RV = "CTGTCTCTTATACACATCT"
READ2_N11 = MEDSB_FW + "N" * 11
AD1_NOMX = "AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG"
AD2_INDEX = "CAAGCAGAAGACGGCATACGAGAT" + "N" * 8 + "GTCTCGTGGGCTCGGAGATGT"


def oligos():
    return (
        oligo("ChILv1-Fw", [Segment("ChIL probe forward", CHIL_V1_FW)],
              five="5′-/5Azido/-", three="-3′"),
        oligo("ChILv1-Rv", [Segment("ChIL probe reverse", CHIL_V1_RV)],
              five="5′-/5Phos/-", three="-/3TMR/-3′"),
        oligo("MEDSB-Fw", [Segment("MEDS-B forward", MEDSB_FW)]),
        oligo("MEDSB-Rv", [Segment("MEDS-B reverse", MEDSB_RV)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("Read2+N11 first-strand synthesis primer",
              [Segment("Read2 plus N11", READ2_N11)]),
        oligo("Ad1_noMX PCR primer", [Segment("Ad1 no index", AD1_NOMX)]),
        oligo("Ad2_index PCR primer", [Segment("Ad2 index", AD2_INDEX)]),
    )
