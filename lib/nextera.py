"""
Nextera / Tn5 tagmentation -- the workhorse adapter chemistry.

Tn5 appears in a large fraction of NGS protocols (SMART-seq onwards, ATAC-seq, scATAC,
Nextera DNA/XT, and as the back end of many single-cell methods), always with the same
19-bp mosaic end and the same two primer entry points. Defined once here.

Sources
-------
Illumina  "Illumina Adapter Sequences", Document # 1000000002694 v22, Sept 2025
          (Nextera / Nextera XT sections)
Reznikoff WS, Ann Rev Genet 2008;42:269-286 -- Tn5 transposition mechanism
Picelli S et al., Genome Res 2014;24(12):2033-40 -- Tn5 production and tagmentation
"""

from __future__ import annotations

from chemdraw import revcomp
from illumina import P5, P7

# --------------------------------------------------------------- core sequences
ME = "AGATGTGTATAAGAGACAG"            # 19-bp mosaic end: the Tn5 binding site
ME_RC = revcomp(ME)                   # CTGTCTCTTATACACATCT -- what you read through into
S5 = "TCGTCGGCAGCGTC"                 # N/S5xx primer entry point (14 nt)
S7 = "GTCTCGTGGGCTCGG"                # N7xx primer entry point (15 nt)
S5_RC = revcomp(S5)                   # GACGCTGCCGACGA
S7_RC = revcomp(S7)                   # CCGAGCCCACGAGAC

# A loaded Tn5 dimer carries two adaptors, each = entry point + ME.
ADAPTOR_S5 = S5 + ME
ADAPTOR_S7 = S7 + ME

# Tagmentation leaves a 9-bp gap on each strand, because the two Tn5 monomers of the
# dimer cut the two strands 9 bp apart. It is filled at 72 C -- the first step of the
# Nextera PCR programme, before any denaturation.
TAGMENTATION_GAP = 9
GAP_FILL_TEMP_C = 72

# ------------------------------------------------------------------- primers
READ1_PRIMER = S5 + ME                        # 33 nt
READ2_PRIMER = S7 + ME                        # 34 nt
INDEX1_PRIMER = ME_RC + S7_RC                # i7 read: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
INDEX2_PRIMER = ME_RC + S5_RC                # i5 read: CTGTCTCTTATACACATCTGACGCTGCCGACGA


def n5xx_primer(i5_index: str) -> str:
    """Nextera (XT) N/S5xx index primer, 5'->3'."""
    return P5 + i5_index + S5


def n7xx_primer(i7_index: str) -> str:
    """Nextera (XT) N7xx index primer, 5'->3'."""
    return P7 + i7_index + S7


# ------------------------------------------------------- what tagmentation produces
# A Tn5 reaction loaded with both adaptors inserts them at random, so a fragment ends up
# with one of three end combinations. Only the heteroduplex amplifies: the others carry
# the same flowcell adapter at both ends and so cannot bridge P5 to P7. This is why the
# scg_lib_structs pages draw every product and label which one survives -- the suppression
# is the point, not an afterthought.
TAGMENTATION_OUTCOMES = (
    ("s5", "s5", False, "P5 at both ends; no P7 site, so it cannot amplify"),
    ("s7", "s7", False, "P7 at both ends; no P5 site, so it cannot amplify"),
    ("s5", "s7", True, "the only amplifiable product"),
)


def amplifiable(end_a: str, end_b: str) -> bool:
    """Does a fragment with these two ends amplify under N5xx + N7xx primers?"""
    return {end_a, end_b} == {"s5", "s7"}
