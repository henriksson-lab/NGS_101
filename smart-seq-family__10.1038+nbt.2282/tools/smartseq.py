"""
Segment definitions for the SMART-seq family.

Two distinct chemistries share this file because they share a back end (Nextera
tagmentation) but differ completely at the front:

  SMART-seq / SMART-seq2   one handle at both ends of the cDNA (the SMART/ISPCR handle),
                           so amplification is single-primer and semi-suppressive. No UMI,
                           no 5' marker -- every fragment is anonymous within a transcript.

  SMART-seq3 / 3xpress     asymmetric handles: the TSO carries the Nextera s5+ME sequence
  / FLASH-seq              plus an 11-bp tag and an 8-bp UMI, while the oligo-dT carries a
                           different handle. A read that contains the tag is known to come
                           from the transcript's 5' end, which is what lets SMART-seq3
                           combine full-length coverage with UMI counting.

Source: https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html
        Ramskold 2012 (Nat Biotechnol 30:777), Picelli 2013 (Nat Methods 10:1096),
        Hagemann-Jensen 2020 (Nat Biotechnol 38:708), Hahaut 2022 (Nat Biotechnol 40:1447).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Segment
from illumina import P5, P7_RC

# =============================================================== SMART-seq / SMART-seq2
# One handle, used for priming, for template switching and for amplification.
DT_LEN, DT_ANCHOR = 30, "VN"
SS2_DT_LINKER = "AC"                                        # between handle and (T)30
SS2_TSO_LINKER = "ACAT"                                     # between handle and G tail
SS2_OLIGO_DT = rt.oligo_dt(DT_LEN, DT_ANCHOR, rt.SMART_HANDLE + SS2_DT_LINKER)
SS2_TSO = rt.SMART_HANDLE + SS2_TSO_LINKER + rt.TSO_G_TAIL_LNA      # ...ACATrGrG+G
SS2_ISPCR = rt.SMART_HANDLE

# =========================================== SMART-seq3 / SMART-seq3xpress / FLASH-seq
SS3_OLIGO_DT_HANDLE = "ACGAGCATCAGCAGCATACGA"               # 5'-biotinylated
SS3_OLIGO_DT = rt.oligo_dt(DT_LEN, DT_ANCHOR, SS3_OLIGO_DT_HANDLE)
SS3_TSO_TAG = "ATTGCGCAATG"                                 # the 11-bp 5' marker
SS3_TSO_ME3 = nx.ME[-8:]                                    # AGAGACAG: the 3' end of ME
SS3_TSO_HANDLE = SS3_TSO_ME3 + SS3_TSO_TAG                  # 3' end of ME, then the tag
SS3_UMI_LEN = 8

# The three family members differ only in what sits between the UMI and the G tail.
SS3_TSO_SPACERS = {
    "SMART-seq3": "",
    "SMART-seq3xpress": "WW",
    "FLASH-seq": "CTAAC",
}
SS3_FWD_PCR = nx.S5 + nx.ME + SS3_TSO_TAG                   # ...AGATGTGTATAAGAGACAGATTGCGCAATG
SS3_REV_PCR = SS3_OLIGO_DT_HANDLE


def ss3_tso(variant: str = "SMART-seq3") -> str:
    """The TSO for one family member, with the UMI written as N's."""
    return rt.tso(SS3_TSO_HANDLE, umi="N" * SS3_UMI_LEN + SS3_TSO_SPACERS[variant])


def _seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def _nextera_i5():
    return [_seg("Illumina P5", P5, "p5"),
            _seg("i5", "N" * 8, None, placeholder=True),
            _seg("s5", nx.S5, "s5"),
            _seg("ME", nx.ME, "me")]


def _nextera_i7():
    return [_seg("ME ", nx.ME_RC, "me"),
            _seg("s7", nx.S7_RC, "s7"),
            _seg("i7", "N" * 8, None, placeholder=True),
            _seg("Illumina P7", P7_RC, "p7")]


def smartseq2_library() -> Construct:
    """SMART-seq / SMART-seq2 final library. Nextera at both ends; no UMI, no 5' marker."""
    return Construct([*_nextera_i5(),
                      _seg("cDNA", "XXXXXXXX...XXXXXXXX", None, placeholder=True),
                      *_nextera_i7()],
                     name="SMART-seq2 library")


def smartseq3_library(variant: str = "SMART-seq3", five_prime: bool = True) -> Construct:
    """SMART-seq3 final library.

    `five_prime=True` gives a fragment that retained the TSO -- it carries the 11-bp tag
    and the UMI, and is identifiable as coming from the transcript's 5' end.
    `five_prime=False` gives an internal fragment, which is indistinguishable from a
    SMART-seq2 fragment and is used only for coverage.
    """
    middle = []
    if five_prime:
        middle = [_seg("5' tag", SS3_TSO_TAG, "r1"),
                  _seg("UMI", "N" * SS3_UMI_LEN, "umi", placeholder=True)]
        spacer = SS3_TSO_SPACERS[variant]
        if spacer:
            middle.append(_seg("", spacer, "r1", placeholder=True))
        middle.append(_seg("", "GGG", "tso"))
    return Construct([*_nextera_i5(), *middle,
                      _seg("cDNA", "XXXXXXXX...XXXXXXXX", None, placeholder=True),
                      *_nextera_i7()],
                     name=f"{variant} library ({'5-prime' if five_prime else 'internal'})")


# ------------------------------------------------------------------ sequencing primers
# By reference: every sequence lives in lib/ (nextera.py), and lib/seqprimers.py computes
# where each lands on the final library. Both chemistries read with stock Nextera primers.
# The SMART-seq3 point is that the TSO supplied SS3_TSO_ME3, the last 8 nt of ME, so the
# stock Read 1 primer (s5 + ME) ends exactly where the 11-bp tag begins.
SS2_SEQ_PRIMERS = [sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2")]
SS3_SEQ_PRIMERS = SS2_SEQ_PRIMERS          # same primers; what Read 1 reports differs


def seq_primer_sets():
    """(final library, declared primers) for every final library this page draws."""
    out = [(smartseq2_library(), SS2_SEQ_PRIMERS)]
    for v in SS3_TSO_SPACERS:
        out.append((smartseq3_library(v), SS3_SEQ_PRIMERS))
        out.append((smartseq3_library(v, five_prime=False), SS3_SEQ_PRIMERS))
    return out
