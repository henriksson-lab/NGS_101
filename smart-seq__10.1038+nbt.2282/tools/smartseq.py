"""
Segment definitions for the SMART-seq family.

Two distinct chemistries share this file because they share a back end (Nextera
tagmentation) but differ completely at the front:

  SMART-seq / SMART-seq2   one handle at both ends of the cDNA (the SMART/ISPCR handle),
                           so amplification is single-primer and semi-suppressive. No UMI,
                           no 5' marker -- every fragment is anonymous within a transcript.

  SMART-seq3 / 3xpress     asymmetric handles: the TSO carries the Nextera s5+ME sequence
                           plus an 11-bp tag and an 8-bp UMI, while the oligo-dT carries a
                           different handle. A read that contains the tag is known to come
                           from the transcript's 5' end, which is what lets SMART-seq3
                           combine full-length coverage with UMI counting.

Sources
-------
Every base below is either a canonical constant from lib/ or a transcription of one of
these primary documents; no sequence here comes from a secondary redrawing of the family.

Ramskold 2012 (Nat Biotechnol 30:777), PMC3467340 Methods
    the CDS oligo-dT primer, 5'-AAGCAGTGGTATCAACGCAGAGTACT(30)VN-3'. The SMARTer II A
    oligo and the kit IS PCR primer are named but never printed.
Picelli 2013 (Nat Methods 10:1096), Supplementary Tables 2 and 4
    the thirteen template-switch oligos tested, and which of them Smart-seq2 uses.
Hagemann-Jensen 2020 (Nat Biotechnol 38:708), Supplementary Table 1 + protocols.io
    smart-seq3-protocol-36wgq5rjxgk5
    the complete Smart-seq3 oligo set, and (condition 1) the Nature Protocols Smart-seq2
    oligo-dT, TSO and forward PCR primer.
Hagemann-Jensen 2022 (Nat Biotechnol 40:1452), PMC9546772 Methods + protocols.io
    smart-seq3xpress-bwh4pb8w
    the Smart-seq3xpress WW-spacer TSO and the ordered oligo table.
Hahaut 2022 (Nat Biotechnol 40:1447), Supplementary Note 2
    FS dT30VN, FS TSO, ISPCR and the FS-UMI TSO series.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "smart-seq2__10.1038+nmeth.2639" / "tools"))

import nextera as nx
import rt
import smartseq2 as s2
import seqprimers as sp
from chemdraw import Construct, Segment, feature
from illumina import P5, P7_RC

UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UMI": UMI_FEATURE, "i5": I5_FEATURE, "i7": I7_FEATURE}

# =============================================================== SMART-seq / SMART-seq2
# One handle, used for priming, for template switching and for amplification.
# =============================================================== FLASH-seq
# Hahaut et al. Supplementary Note 2: the main FLASH-seq branch retains the common
# SMART handle and has no UMI. The ordered oligos are 5'-biotinylated. Every string
# below is a verbatim transcription of that table.
FS_OLIGO_DT = rt.oligo_dt(s2.DT_LEN, s2.DT_ANCHOR, rt.SMART_HANDLE + s2.DT_LINKER)
FS_TSO = rt.SMART_HANDLE + s2.DT_LINKER + rt.TSO_G_TAIL
FS_ISPCR = rt.SMART_HANDLE
# The optional FS-UMI branch, not the main chemistry this page draws: Supplementary
# Note 2 lists nine spacer variants and the paper selects CTAAC.
FS_UMI_TSO_CTAAC = rt.SMART_HANDLE + "N" * 8 + "CTAAC" + rt.TSO_G_TAIL

# =========================================== SMART-seq3 / SMART-seq3xpress
# Hagemann-Jensen 2020 Supplementary Table 1 and the author protocol both print
# /5Biosg/ACGAGCATCAGCAGCATACGA(T)30VN and /5Biosg/AGAGACAGATTGCGCAATG(N)8rGrGrG.
SS3_OLIGO_DT_HANDLE = "ACGAGCATCAGCAGCATACGA"               # 5'-biotinylated
SS3_OLIGO_DT = rt.oligo_dt(s2.DT_LEN, s2.DT_ANCHOR, SS3_OLIGO_DT_HANDLE)
SS3_TSO_TAG = "ATTGCGCAATG"                                 # the 11-bp 5' marker
SS3_TSO_ME3 = nx.ME[-8:]                                    # AGAGACAG: the 3' end of ME
SS3_TSO_HANDLE = SS3_TSO_ME3 + SS3_TSO_TAG                  # 3' end of ME, then the tag
SS3_UMI_LEN = 8

# The SMART-seq3 derivatives differ in what sits between the UMI and the G tail.
# Hagemann-Jensen 2022 Methods (PMC9546772) contrasts the "Original Smartseq3 TSO"
# (...NNNNNNNNrGrGrG) with the xpress "Improved TSO" (...NNNNNNNNWWrGrGrG).
SS3_TSO_SPACERS = {
    "SMART-seq3": "",
    "SMART-seq3xpress": "WW",
}
SS3_FWD_PCR = nx.S5 + nx.ME + SS3_TSO_TAG                   # ...AGATGTGTATAAGAGACAGATTGCGCAATG
SS3_REV_PCR = SS3_OLIGO_DT_HANDLE


def ss3_tso(variant: str = "SMART-seq3") -> str:
    """The TSO for one family member, with the UMI written as N's."""
    return rt.tso(SS3_TSO_HANDLE, umi="N" * SS3_UMI_LEN + SS3_TSO_SPACERS[variant])


def _seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
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
        middle.append(_seg("template-switch GGG", "GGG", "tso"))
    return Construct([*_nextera_i5(), *middle,
                      _seg("cDNA", "XXXXXXXX...XXXXXXXX", None, placeholder=True),
                      *_nextera_i7()],
                     name=f"{variant} library ({'5-prime' if five_prime else 'internal'})")


# ------------------------------------------------------------------ sequencing primers
# By reference: every sequence lives in lib/ (nextera.py), and lib/seqprimers.py computes
# where each lands on the final library. Both chemistries read with stock Nextera primers.
# The SMART-seq3 point is that the TSO supplied SS3_TSO_ME3, the last 8 nt of ME, so the
# stock Read 1 primer (s5 + ME) ends exactly where the 11-bp tag begins.
SS3_SEQ_PRIMERS = s2.SEQ_PRIMERS          # same primers; what Read 1 reports differs


def seq_primer_sets():
    """(final library, declared primers) for every final library this page draws."""
    out = [(s2.final_library(), s2.SEQ_PRIMERS)]
    for v in SS3_TSO_SPACERS:
        out.append((smartseq3_library(v), SS3_SEQ_PRIMERS))
        out.append((smartseq3_library(v, five_prime=False), SS3_SEQ_PRIMERS))
    return out
