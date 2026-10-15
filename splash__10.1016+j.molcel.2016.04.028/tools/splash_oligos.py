"""Published SPLASH adapter and RT-primer oligos."""
from chemdraw import Segment, oligo
import illumina as il


RNA_LINKER = "CTGTAGGCACCATCAAT"
RT_LEFT = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC"
RT_MIDDLE = "CACTCA"
RT_RIGHT = "TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG"
UNIVERSAL_PCR_PRIMER = il.NEBNEXT_UNIVERSAL_PRIMER
INDEX1_PCR_PRIMER = il.P7 + "CGTGAT" + il.TRUSEQ_READ2


def oligos():
    return (
        oligo("3-prime RNA adapter", [Segment("3-prime adapter", RNA_LINKER)],
              five="5′-/5rApp/-", three="-/3ddC/-3′"),
        oligo("96-nt RT primer", [
            Segment("RT primer 5-prime arm", RT_LEFT),
            Segment("internal spacer 1", "[iSp18]", placeholder=True),
            Segment("central tract", RT_MIDDLE),
            Segment("internal spacer 2", "[iSp18]", placeholder=True),
            Segment("adapter-binding 3-prime arm", RT_RIGHT),
        ]),
        oligo("NEB E7335L Universal PCR Primer", [
            Segment("universal PCR primer", UNIVERSAL_PCR_PRIMER[:-1]),
        ], three="-s-T-3′"),
        oligo("NEB E7335L Index 1 PCR Primer", [
            Segment("index 1 PCR primer", INDEX1_PCR_PRIMER[:-1]),
        ], three="-s-T-3′"),
    )
