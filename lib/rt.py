"""
Reverse transcription and template switching -- the shared RNA-side building blocks.

Two ideas recur across nearly every mRNA protocol:

1. **Priming.** An oligo-dT primer anneals to the poly(A) tail and MMLV-family reverse
   transcriptase copies the mRNA into cDNA. The 3'-most two bases are usually anchored
   (`VN`: V = A/C/G, N = any) so the primer sits at the poly(A) junction rather than
   sliding along the tail. Random priming (`NNNNNN`) is used instead when the target has
   no poly(A) -- bacterial mRNA, or to capture internal sequence.

2. **Template switching.** On reaching the 5' end of the mRNA, MMLV's terminal
   transferase activity adds a few untemplated nucleotides -- predominantly **C**, usually
   three. A template-switching oligo (TSO) ending in **rGrGrG** base-pairs with that CCC
   overhang, and the polymerase switches template and copies the TSO. The result is a
   known handle at the 5' end of the cDNA, so the full transcript is flanked by defined
   sequence without any ligation step.

   The 3'-terminal G's are ribonucleotides (`rG`), or in SMART-seq2 an LNA (`+G`), both of
   which raise the duplex stability of the short G:C contact enough to make the switch
   efficient.

Sources
-------
Zhu YY et al., Biotechniques 2001;30(4):892-7 -- template switching / SMART
Picelli S et al., Nat Methods 2013;10(11):1096-8 -- SMART-seq2, LNA TSO
Hagemann-Jensen M et al., Nat Biotechnol 2020;38(6):708-14 -- SMART-seq3, UMI in the TSO
"""

from __future__ import annotations

# The SMART / ISPCR handle. Shared by SMART-seq, SMART-seq2, SPLiT-seq, Drop-seq and the
# Takara SMARTer kits -- if a protocol says "ISPCR primer" or "SMART handle", it is this.
SMART_HANDLE = "AAGCAGTGGTATCAACGCAGAGT"      # 23 nt

# How MMLV ends a first strand: untemplated nucleotides, predominantly C, usually 3.
UNTEMPLATED_TAIL = "CCC"
TSO_G_TAIL = "rGrGrG"                          # pairs with the CCC overhang
TSO_G_TAIL_LNA = "rGrG+G"                      # SMART-seq2: last G is a locked nucleic acid


def oligo_dt(length: int = 30, anchor: str = "VN", handle: str = "") -> str:
    """An anchored oligo-dT primer, 5'->3'. `anchor=""` gives an unanchored tail."""
    return f"{handle}{'T' * length}{anchor}"


def tso(handle: str, tail: str = TSO_G_TAIL, umi: str = "") -> str:
    """A template-switching oligo: handle, optional UMI, then the 3' G tail."""
    return f"{handle}{umi}{tail}"
