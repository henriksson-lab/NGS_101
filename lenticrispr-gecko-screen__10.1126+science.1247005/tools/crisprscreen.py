"""
Pooled CRISPR screening readout primers, as data.

The single most important fact about this field: **different protocols use different
primers, and they are not interchangeable.** Different labs put the index on different
ends, anneal to different parts of the vector, and disagree about whether a given pair
even produces a usable product on a given backbone.

Every claim here is checkable against the real maps in ../ref/plasmids/ -- see selftest.py,
which computes each amplicon rather than quoting it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import revcomp
import crispr
from crispr import SCAFFOLD_V1, U6_3PRIME, U6_PLUS1

# ------------------------------------------------------------- Illumina blocks
# No sequence of Illumina's is retyped here: the blocks come from lib/illumina.py.
# The full-length P5-side adapter is P5 + TruSeq read 1, overlapping by the 4 nt
# (ACAC) that both constants carry -- 58 nt, as every published readout primer prints it.
_P5_R1_OVERLAP = 4
P5_READ1 = il.P5 + il.TRUSEQ_READ1[_P5_R1_OVERLAP:]                         # 58 nt
# == il.TRUSEQ_P5_FULL; selftest.py pins the composition against it.
P7 = il.P7                                                                  # 24 nt
READ2_IDX = il.TRUSEQ_READ2                                                 # 34 nt

# ------------------------------------------------------- vector annealing sites
# Which part of the vector a primer grabs decides which backbones it works on.
# The U6 sites are all nested 3' ends of one stretch of promoter -- so it is written
# once (the part 5' of lib/crispr.py's U6_3PRIME) and each published primer's site is
# cut from it by length. selftest.py checks the stretch is really in all three maps.
U6_UPSTREAM_OF_3PRIME = crispr.U6_UPSTREAM[1:]           # the 22 nt the GPP sites need
U6_REGION = U6_UPSTREAM_OF_3PRIME + U6_3PRIME          # 33 nt, ends ...CACC (no +1 G)
_U6_WITH_G = U6_REGION + U6_PLUS1
SITE_U6_SHORT = _U6_WITH_G[-22:]                     # 22 nt, GPP P5_ARGON
SITE_U6_LONG = U6_REGION                             # 33 nt, Zhang/Joung, no +1 G
SITE_U6_24 = _U6_WITH_G[-24:]                        # 24 nt, Shalem / Zhang F01-F12
SITE_CPPT = "TCTACTATTCTTTCCCCTGCACTGT"              # HIV-1 cPPT/CTS  (GPP "KERMIT")
SITE_EFS = "CCAATTCCCACTCCTTTCAAGACCT"               # EFS/EF-1a leader (GPP "BEAKER")
SITE_SCAFFOLD = revcomp(SCAFFOLD_V1[50:73])          # sgRNA scaffold  (Joung KO-Rev)

# ---------------------------------------------------------------- demo spacers
# One pair of example spacers for the whole page, so a figure and a size can never be
# drawn from different ones. The leading letter is the point: it decides the +1 G.
DEMO_SPACER_G = "GTCGCTGAGTACTTCGAAAT"      # starts with G: no +1 appended
DEMO_SPACER = "ATCGATCGATCGATCGATCG"        # does not: +1 G appended (Brunello-style)

# ------------------------------------------------------------------ staggers
# Every lab solves the same problem -- the bases between the sequencing primer and the
# variable spacer are identical in every cluster -- and they solve it differently.
GPP_STAGGERS = ("", "C", "GC", "AGC", "CAAC", "TGCACC", "ACGCAAC", "GAAGACCC")  # no 5
# NB: SITE_U6_LONG starts with G, so the stagger stops one base short of where the
# published primer strings appear to split. Getting this wrong makes every stagger 1 nt
# too long -- the selftest catches it.
JOUNG_STAGGERS = ("TAAGTAGAG", "ATCATGCTTA", "GATGCACATCT", "CGATTGCTCGAC",
                  "TCGATAGCAATTC", "ATCGATAGTTGCTT", "GATCGATCCAGTTAG",
                  "CGATCGATTTGAGCCT", "ACGATCGATACACGATC", "TACGATCGATGGTCCAGA")


def gpp_p5(stagger: str = "") -> str:
    """Broad GPP P5_ARGON. Unindexed; the index lives on the P7 side."""
    return P5_READ1 + stagger + SITE_U6_SHORT


def gpp_p7(index: str, site: str = SITE_CPPT) -> str:
    """Broad GPP P7. `site` picks the vector family: KERMIT (cPPT) or BEAKER (EFS)."""
    return P7 + index + READ2_IDX + site


def joung_fwd(stagger: str) -> str:
    """Zhang/Joung one-step forward. Unindexed; staggers are 9-18 nt."""
    return P5_READ1 + stagger + SITE_U6_LONG


def joung_rev(index: str) -> str:
    """Zhang/Joung one-step reverse. Anneals in the SCAFFOLD, so it is vector-independent."""
    return P7 + index + READ2_IDX + SITE_SCAFFOLD


# Shalem 2014's two-step readout. Step 1 carries no Illumina sequence at all.
SHALEM_F1 = "AATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG"
# Its 3' 15 nt are the start of the same cPPT site KERMIT and the v2 adaptor use, so they
# are taken from SITE_CPPT rather than retyped -- which is also why Shalem's pair fails on
# lentiCRISPRv2 for exactly the reason KERMIT does.
SHALEM_R1 = "CTTTAGTTTGTATGTCTGTTGCTATTATG" + SITE_CPPT[:15]

# The Zhang "v2 adaptor": lentiCRISPRv2 has no cPPT downstream of the guide, so this
# grafts the missing priming site on before the standard F01-F12 primers can be used.
V2ADAPTOR_F = SHALEM_F1
V2ADAPTOR_GRAFT = SITE_CPPT                 # non-templated 5' tail: the site it installs
V2ADAPTOR_ANNEAL = "TGTGGGCGATGTGCGCTCTG"   # 3' half, the only part present in the vector
V2ADAPTOR_R = V2ADAPTOR_GRAFT + V2ADAPTOR_ANNEAL

GPP_P7_INDEX_A01 = "CGGTTCAA"       # carried as-is in the oligo; the i7 read reports its
                                    # reverse complement -- computed, not quoted (the
                                    # sequencing-primer table on the page shows it)
JOUNG_REV_INDEX_1 = "TCGCCTTG"

# "Works" means a short, amplifiable product. A pair whose sites sit in the wrong order
# on the circle still finds both primers, but the product runs the long way round -- which
# is why selftest.py checks LENGTH, not merely that a product exists.
VECTOR_READOUT = {
    "lentiCRISPR v1":  {"gpp": "KERMIT", "shalem": True,  "joung": True},
    "lentiGuide-Puro": {"gpp": "KERMIT", "shalem": True,  "joung": True},
    "lentiCRISPRv2":   {"gpp": "BEAKER", "shalem": False, "joung": True},
}

MAX_USABLE_AMPLICON = 1000   # longer than this is a wrap-around artefact, not a product

# ------------------------------------------------------- sequencing primers
# Nothing custom is loaded for a screen: every readout primer above ends in a stock
# TruSeq handle, so all four reads are primed by the stock Illumina primers. The sites
# are not asserted here -- lib/seqprimers.py locates each one on the finished amplicon
# and raises if it is not there, so the page cannot claim a read it cannot prime.
_GPP = "Broad GPP sequencing protocol for pooled screens"
SEQ_PRIMERS = [
    sp.TRUSEQ["R1"],
    sp.TRUSEQ["I1"],
    sp.custom("Index 2 (i5)", "TruSeq Index 2, reverse-complement workflow",
              il.INDEX2_PRIMER_RC, _GPP,
              "These libraries are SINGLE-indexed -- the index sits on the P7 side only. "
              "The i5 primer still has its site, but there is no i5 index to report: an "
              "i5 cycle would read into P5, so index 2 is not run."),
    sp.TRUSEQ["R2"],
]
