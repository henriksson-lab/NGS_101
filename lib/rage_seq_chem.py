"""Published RAGE-Seq cDNA amplification primers."""
from chemdraw import Segment, oligo


TSO_PRIMER = "AAGCAGTGGTATCAACGCAGAGT"
R1_PRIMER = "CTACACGACGCTCTTCCGATCT"


def oligos():
    return (
        oligo("TSO primer", [Segment("10x template-switch PCR primer", TSO_PRIMER)]),
        oligo("R1 primer", [Segment("10x Read-1-side PCR primer", R1_PRIMER)]),
    )
