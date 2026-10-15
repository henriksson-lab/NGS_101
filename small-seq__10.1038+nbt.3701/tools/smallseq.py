"""
Segment definitions and constants for Small-seq -- single-cell small-RNA sequencing.

Source
------
Hagemann-Jensen M, Abdullayev I, Sandberg R, Faridani OR.
"Small-seq for single-cell small-RNA sequencing." Nat Protoc 2018;13(10):2407-2424.
doi:10.1038/s41596-018-0049-y          (extracted text: pdf/hagemann.txt, gitignored)
Original method: Faridani OR et al., Nat Biotechnol 2016;34(12):1264-1266.

Why this chemistry is in the repo
--------------------------------
It is the only protocol here that is *not* built on TruSeq DNA adapters and not built on
either template switching or tagmentation. Four things are new:

  1. **TruSeq Small RNA adapters** (RA5 / RA3) rather than the TruSeq DNA adapter, so the
     read-1 landing site is a completely different sequence (see 02_sequencing_primer.md).
  2. **Sequential ligation** -- 3' adapter first, then a digestion step that destroys the
     unligated 3' adapter, then the 5' adapter. No polymerase is involved in adapter
     attachment at all.
  3. **A UMI in the 5' adapter**, not in a TSO and not on a probe. It is the first thing
     sequenced, before any of the insert.
  4. **A masking oligonucleotide** that occludes the 3' end of 5.8S rRNA so that the 3'
     adapter cannot be ligated to it. Depletion by hybridisation and *inaction*, rather
     than by pulldown.

Evidence marking (repo convention)
----------------------------------
  [green]  verbatim from the paper
  [yellow] derived or inferred here
  [red]    not published -- not in the paper's main text at all
  [check]  computed and asserted in tools/selftest.py

Every real sequence below is [green] from Reagent setup (p. 2414), which states
"All oligonucleotides are listed in the 5' to 3' direction". Line-wrap spaces in the
extracted text have been closed up; nothing else was changed.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, revcomp
from illumina import P5, P7, P7_RC, TRUSEQ_SMALL_RNA_RA5, TRUSEQ_SMALL_RNA_RA3
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Small-seq — single-cell small-RNA sequencing"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nbt.3701">Faridani et al. (2016)</a>; detailed protocol: <a href="https://doi.org/10.1038/s41596-018-0049-y">Hagemann-Jensen et al. (2018)</a>.'

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "whitelist",
                     whitelist="small-seq-srx-192-index-set")
FEATURES = {"UMI": UMI_FEATURE, "i7": I7_FEATURE}

# ======================================================================= oligos
# The stock RA5 body and RA3 are imported from Illumina's adapter-sequences document;
# Small-seq's RTP differs from the stock kit primer and remains local.
# All [green] verbatim, Reagent setup. Modifications are carried separately as _MODS
# strings so that the base sequence stays computable (tm(), revcomp(), len()).

# --- 5' adapter: RNA, 5'-amino-blocked, free 3'-OH, UMI + CA at the 3' end -----------
# RA5 (NH2-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrCrHrHrHrHrHrHrHrHrCrA)
RA5_HANDLE = TRUSEQ_SMALL_RNA_RA5                  # 26 nt, written as DNA for arithmetic
RA5_UMI_LEN = 8                                 # [green] "a stretch of eight 'H'"
RA5_UMI_BASE = "H"                              # [green] rH = rA, rU or rC
RA5_UMI_DIVERSITY = 3 ** RA5_UMI_LEN            # [green] "6,561 different UMIs"
RA5_LINKER = "CA"                               # [green] rCrA, separates UMI from insert
RA5_MODS = "NH2 (aminolink C6) at the 5' end; 3'-OH; every base a ribonucleotide"

# --- 3' adapter: DNA, 5'-adenylated, 3'-blocked ---------------------------------------
# RA3 (rAppTGGAATTCTCGGGTGCCAAGG-ddC)
RA3 = TRUSEQ_SMALL_RNA_RA3                         # 21 nt
RA3_MODS = "5'-rApp (pre-adenylated); 3'-ddC"

# --- RT primer, which is also the reverse primer of PCR 1 -----------------------------
# RTP (biotin-CCTTGGCACCCGAGAATTCCrA)
RTP = "CCTTGGCACCCGAGAATTCCA"                   # 21 nt
# Kept verbatim although it equals revcomp(RA3): the paper prints the two independently, and
# that equality (asserted in selftest) is what the digestion step depends on.
RTP_MODS = "5'-biotin (blocks lambda exonuclease); 3'-terminal rA"

# --- forward primer of both PCRs ------------------------------------------------------
# RP1 (AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA)   <- as printed, 50 nt
# Built from its parts so that P5 (lib/illumina.py) and RA5 are not typed twice; the
# selftest asserts the result equals the printed string in the comment above.
RP1_ANNEAL_NT = 21          # [check] its 3' 21 nt == RA5_HANDLE[:21]; the 5' 29 are P5
RP1 = P5 + RA5_HANDLE[:RP1_ANNEAL_NT]           # 50 nt

# --- 5.8S rRNA masking oligonucleotide ------------------------------------------------
# [green] ATCGGCAAGCGACGCTCAGA CAGGCGTAGCCCCGGG AGGAACCCGGGGC CGCAAGTGCGTTCGAAGTGTCGATGAT-biotin
MASK_58S = ("ATCGGCAAGCGACGCTCAGA"
            "CAGGCGTAGCCCCGGG"
            "AGGAACCCGGGGC"
            "CGCAAGTGCGTTCGAAGTGTCGATGAT")      # 76 nt
MASK_58S_MODS = "3'-biotin, to block extension"
# [yellow] what it must be annealing to. The paper says only that it "targets the 3' end
# of the rRNA molecule" and "binds to the 3' end of the highly abundant 5.8S rRNA"; it
# never prints the rRNA. Checking this base-by-base needs a reference sequence, which is
# data and therefore does not live in this repo (README, "Data policy"). Recorded so the
# next person can BLAST it rather than retype it.
MASK_58S_TARGET = revcomp(MASK_58S)
RRNA_58S_NT = 156           # [yellow] textbook length of mature human 5.8S; NOT in the paper

# --- index primers: detailed-protocol Supplementary Table 1 --------------------------
# The supplement prints all 192 full primers. SR001 is the representative sequence used
# in the molecular drawing; the other family members change only the eight-base window.
SRX_INDEX_NT = 8                                # [green] "barcodes of 8 bp"
SRX_N_INDICES = 192                             # [green] "for 192 samples" (two 96 plates)
SRX_CONC_UM = 50                                # [green] dissolved to 50 uM, 1 ul per well
SRX_001_PRIMER_INDEX = "TTCGTGAT"               # sequence inside the ordered primer
SRX_001_INDEX_READ = "ATCACGAA"                 # table's reported barcode, revcomp above
SRX_COMMON_LINK = "GTGACTGGAGTT"                # Supplementary Table 1, between index and RTP foot
SRX_SEQUENCE = P7 + SRX_001_PRIMER_INDEX + SRX_COMMON_LINK + RTP
SRX_P7_NT = len(P7)                             # 24 [check] = len(illumina.P7)
SRX_LINKER_NT = len(SRX_COMMON_LINK)             # 12, source-written SR001 primer
SRX_RIGHT_ARM_NT = SRX_P7_NT + SRX_INDEX_NT + SRX_LINKER_NT      # 44

# ================================================================= protocol constants
LYSIS_TEMP_C, LYSIS_MIN = 72, 20                # [green] Step 9
LIG3_TEMP_C, LIG3_HOURS = 30, 6                 # [green] Step 12, then 4 C for 10 h
LIG5_TEMP_C, LIG5_MIN = 37, 60                  # [green] Step 18
RT_TEMP_C, RT_MIN = 42, 60                      # [green] Step 21, SuperScript II
PCR1_ANNEAL_C, PCR1_CYCLES = 60, 13             # [green] Step 24, cycles 2-14
PCR2_ANNEAL_C, PCR2_CYCLES = 67, 13             # [green] Step 31, cycles 2-14
PCR1_PRIMERS = ("RP1", "RTP (leftover from the digestion step)")   # [green]
PCR2_PRIMERS = ("RP1", "SRX index primer")                         # [green]

READ_LEN = 51                                   # [green] Step 38, "a single read of 51 bp"
TRIM_5P = RA5_UMI_LEN + len(RA5_LINKER)         # 10 nt removed before mapping [green]
SMALLRNA_MAX_NT = 40                            # [green] maxRlen, "40 bp or shorter"
PRECURSOR_MIN_NT = 41                           # [green] minRlen, reads of maximum length
SMALLRNA_MIN_NT = 18                            # [green] cutadapt --minimum-length 18

PIPPIN_WINDOW_BP = (130, 160)                   # [green] Step 35A(i) start/end
GEL_CUT_BP = (120, 200)                         # [green] Step 35B(iv)
SIZE_SELECT_TARGET_BP = (145, 155)              # [green] "aiming to cut targets in"
PIPPIN_RESULT_BP = (140, 180)                   # [green] Anticipated results, Fig. 3d
RRNA_PEAK_BP = (270, 290)                       # [green] Fig. 3b, mask omitted
LIBRARY_PROFILE_BP = (100, 300)                 # [green] Anticipated results, Fig. 3a

TYPICAL_MIRNA_NT = 22                           # [yellow] not stated; used for sizing only

# Enzymes, in order. [green] Reagent tables of Steps 10, 13, 16, 19, 22.
ENZYMES = [
    ("3' ligation", "T4 RNA ligase 2, truncated KQ",
     "no ATP, so only the pre-adenylated adapter can be ligated -- small RNAs cannot be "
     "ligated to each other"),
    ("digestion", "5' deadenylase",
     "strips the 5'-rApp from unligated RA3, leaving a 5'-phosphate"),
    ("digestion", "lambda exonuclease",
     "digests the 5'-phosphorylated strand of a duplex: RTP anneals to free RA3 and the "
     "exonuclease eats the RA3, while RTP's 5'-biotin protects it"),
    ("5' ligation", "T4 RNA ligase 1 + ATP",
     "needs a 5'-phosphate on the acceptor, which capped mRNA does not have"),
    ("RT", "SuperScript II, in Taq buffer",
     "Taq buffer rather than RT buffer because excess MgCl2 is carried over"),
    ("PCR", "Phusion Hot Start II", "high fidelity, both rounds"),
]


# ======================================================================= constructs
def _seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def ra5_oligo() -> str:
    """RA5 as ordered, written as DNA (H = rA/rU/rC placeholder)."""
    return RA5_HANDLE + RA5_UMI_BASE * RA5_UMI_LEN + RA5_LINKER


def insert_placeholder(n: int) -> str:
    """A stand-in insert of n characters -- never complemented as if it were bases."""
    if n <= 8:
        return "X" * n
    return "X" * (n - 6) + "...XXX"


def ligated_rna(insert_nt: int = TYPICAL_MIRNA_NT) -> Construct:
    """Small RNA after both ligations: RA5 + UMI + CA + insert + RA3. Single-stranded."""
    return Construct([
        _seg("RA5", RA5_HANDLE, "r1"),
        _seg("UMI", RA5_UMI_BASE * RA5_UMI_LEN, "umi", placeholder=True),
        _seg("CA", RA5_LINKER, "tso"),
        _seg("small RNA", insert_placeholder(insert_nt), None, placeholder=True),
        _seg("RA3", RA3, "r2"),
    ], name="ligated small RNA")


def free_ra3_rtp_duplex() -> Scene:
    """The source-written RTP binds the complete free RA3 after rApp removal.

    Hagemann-Jensen et al. 2018, p. 2410 and Reagent setup p. 2414: deadenylase
    leaves a 5′ phosphate on RA3; the 5′-biotinylated RTP is its reverse
    complement and protects itself while lambda exonuclease removes RA3.
    """
    sc = Scene()
    sc.strand("free RA3", [Segment("RA3", RA3, "r2")], mod5="p", mod3="ddC")
    sc.anneal("RTP", [Segment("RTP", RTP, "r2")], to="free RA3",
              pair=("RTP", "RA3"), mod5="biotin")
    return sc


def library(insert_nt: int = TYPICAL_MIRNA_NT) -> Construct:
    """The final, indexed library.

    SR001 is representative of the 192-member indexed-primer family. The top-strand
    sequence contains the reverse complements of the SRX tail and primer index because
    SRX extends leftwards from its RTP foot on RA3.
    """
    return Construct([
        _seg("P5", P5, "p5"),
        _seg("RA5", RA5_HANDLE, "r1"),
        _seg("UMI", RA5_UMI_BASE * RA5_UMI_LEN, "umi", placeholder=True),
        _seg("CA", RA5_LINKER, "tso"),
        _seg("small RNA", insert_placeholder(insert_nt), None, placeholder=True),
        _seg("RA3", RA3, "r2"),
        _seg("SRX common link", revcomp(SRX_COMMON_LINK), "r3"),
        _seg("i7", SRX_001_INDEX_READ, "cbc",
             note="SR001 representative; 192 full primers are published"),
        _seg("P7", P7_RC, "p7"),
    ], name="Small-seq library")


def pcr1_product(insert_nt: int = TYPICAL_MIRNA_NT) -> Construct:
    """Exact first-round amplicon before the indexed SRX tail is installed."""
    return Construct([
        _seg("P5", P5, "p5"),
        _seg("RA5", RA5_HANDLE, "r1"),
        _seg("UMI", RA5_UMI_BASE * RA5_UMI_LEN, "umi", placeholder=True),
        _seg("CA", RA5_LINKER, "tso"),
        _seg("small RNA", insert_placeholder(insert_nt), None, placeholder=True),
        _seg("RA3", RA3, "r2"),
    ], name="Small-seq PCR-1 amplicon")


def library_bp(insert_nt: int) -> int:
    """Library length for a given insert, from the segment arithmetic alone."""
    return (len(P5) + len(RA5_HANDLE) + RA5_UMI_LEN + len(RA5_LINKER)
            + insert_nt + len(RA3) + SRX_RIGHT_ARM_NT)


LEFT_ARM_NT = len(P5) + len(RA5_HANDLE) + RA5_UMI_LEN + len(RA5_LINKER) + len(RA3)   # 86
ADAPTER_DIMER_BP = library_bp(0)                 # [yellow] 130 bp -- see 01_small-seq.md


# ============================================================== sequencing primers
# By reference: sequences live in lib/ (TruSeq DNA) or above (the Small RNA oligos). The
# section on the page computes where each lands on library() and what it reads first.
#
# TODO(lib): move to illumina.py -- the TruSeq Small RNA block (RA5, RA3, RTP, RP1, RPI
# and the Small RNA Sequencing / Index primers), verbatim from the Illumina adapter
# document. Until then the read-1 entry below is the site the CONSTRUCT requires (RA5, 26
# nt, ending where the UMI begins), not a vendor primer sequence.
READ1_SITE = RA5_HANDLE

SEQ_PRIMERS = [
    sp.custom("Read 1", "TruSeq Small RNA read-1 site (RA5) -- required by the construct",
              READ1_SITE, "derived: [green] \"Sequencing starts from the UMI\"",
              "The paper names no sequencing primer. Any usable one is a suffix of P5 + RA5 "
              "ending on RA5's last base; the vendor Small RNA read-1 primer is not yet in "
              "lib/ to check against."),
    sp.custom("Read 1", "RP1 re-used as read-1 primer (counter-example)", RP1,
              "Reagent setup, p. 2414",
              "Matches, but stops 5 nt short of RA5's 3' end, so read 1 would start CGATC "
              "instead of on the UMI -- contradicting the paper. Not what was used."),
    sp.mismatching(sp.TRUSEQ["R1"], "no TruSeq DNA adapter in this library; a run loaded "
                                    "only with it reads nothing"),
    sp.custom("Index 1 (i7)", "TruSeq Index 1 (13-nt 3′ annealing site)",
              sp.TRUSEQ["I1"].seq, sp.TRUSEQ["I1"].source,
              "Small-seq's 12-nt SRX common tail and the next RA3 base pair with "
              "the primer's 13-nt 3′ end; the index read is printed in the supplement.",
              min_anneal=13),
    sp.mismatching(sp.TRUSEQ["I2"], "single-indexed: RP1 puts P5 directly on RA5, so there "
                                    "is no i5 index and no TruSeq i5 primer site"),
    sp.mismatching(sp.TRUSEQ["R2"], "single-read protocol (Step 38, \"a single read of 51 "
                                    "bp\"); no TruSeq DNA read-2 site"),
]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("single_cell_rna","total RNA from one lysed cell",
                    MaterialKind.BIOLOGICAL,primary=True)
    mask=graph.input("rrna_mask","3′-biotinylated 5.8S rRNA masking oligo",
                     MaterialKind.OLIGO)
    masked=graph.product("masked_rna","RNA with the 5.8S 3′ end occluded",
                         MaterialKind.MOLECULAR)
    graph.react("mask_rrna","Hybridize the masking oligo to 5.8S rRNA",
                StepKind.HYBRIDIZATION,
                ((rna,InputRole.TEMPLATE),(mask,InputRole.STRUCTURAL)),(masked,))
    ra3=graph.input("ra3","pre-adenylated, 3′-blocked RA3 adapter",MaterialKind.ADAPTER)
    rnl2=graph.input("rnl2","T4 RNA ligase 2 truncated KQ",MaterialKind.CATALYST)
    lig3=graph.product("ligated_3p","RA3-ligated small RNA",MaterialKind.MOLECULAR)
    free_ra3=graph.product("free_ra3","unligated RA3 adapter",MaterialKind.MOLECULAR)
    graph.react("ligate_3p","Ligate RA3 to accessible RNA 3′ ends",StepKind.LIGATION,
                ((masked,InputRole.TEMPLATE),(ra3,InputRole.STRUCTURAL),
                 (rnl2,InputRole.CATALYST)),(lig3,free_ra3),
                geometry=AdapterLigationGeometry.unavailable(
                    "This is a pre-adenylated single-strand RNA ligation; no duplex overhang exists."))
    rtp=graph.input("rtp","biotin-blocked RTP oligo",MaterialKind.OLIGO)
    deadenylase=graph.input("deadenylase","5′ deadenylase",MaterialKind.CATALYST)
    deadenylated=graph.product("deadenylated_ra3","deadenylated free RA3–RTP duplex",
                               MaterialKind.MOLECULAR,state=free_ra3_rtp_duplex())
    graph.react("deadenylate","Remove rApp from unligated RA3",StepKind.DEPHOSPHORYLATION,
                ((free_ra3,InputRole.TEMPLATE),(rtp,InputRole.STRUCTURAL),
                 (deadenylase,InputRole.CATALYST)),(deadenylated,))
    lambda_exo=graph.input("lambda_exonuclease","lambda exonuclease",MaterialKind.CATALYST)
    adapter_waste=graph.product("adapter_waste","digested unligated RA3",MaterialKind.WASTE)
    graph.react("digest_free_adapter","Digest free RA3 annealed to protected RTP",
                StepKind.DIGESTION,
                ((deadenylated,InputRole.TEMPLATE),(lambda_exo,InputRole.CATALYST)),
                (adapter_waste,))
    graph.discard(adapter_waste,"unligated 3′ adapter is enzymatically destroyed")
    ra5=graph.input("ra5","UMI-bearing RA5 adapter",MaterialKind.ADAPTER)
    rnl1=graph.input("rnl1","T4 RNA ligase 1",MaterialKind.CATALYST)
    doubly=graph.product("doubly_ligated","RA5/RA3-ligated small RNA",MaterialKind.MOLECULAR,
                         state=ligated_rna())
    graph.react("ligate_5p","Ligate the UMI-bearing RA5 adapter",StepKind.LIGATION,
                ((lig3,InputRole.TEMPLATE),(ra5,InputRole.STRUCTURAL),
                 (rnl1,InputRole.CATALYST)),(doubly,),
                geometry=AdapterLigationGeometry.unavailable(
                    "This is a single-strand RNA ligation; no duplex overhang exists."))
    rt=graph.input("reverse_transcriptase","SuperScript II reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","small-RNA first-strand cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Reverse-transcribe from RTP",StepKind.REVERSE_TRANSCRIPTION,
                ((doubly,InputRole.TEMPLATE),(rtp,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    rp1=graph.input("rp1","RP1 forward PCR primer",MaterialKind.OLIGO)
    pol=graph.input("polymerase","Phusion Hot Start II polymerase",MaterialKind.CATALYST)
    pcr1_state=pcr1_product()
    pcr1=graph.product("pcr1_product","first-round small-RNA amplicon",MaterialKind.MOLECULAR,
                       state=pcr1_state)
    graph.react("pcr1","Amplify with RP1 and leftover RTP",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (rtp,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(pcr1,),
                geometry=PCRGeometry.exact(
                    pcr1_state, PCRPrimer("RP1", RP1, anneal_nt=RP1_ANNEAL_NT),
                    PCRPrimer("RTP", RTP),
                    statement="The detailed protocol prints RP1 and RTP; their inward-facing sites are computed on the first-round amplicon."))
    srx=graph.input("srx","SRX indexed reverse primer",MaterialKind.OLIGO,
                    state=SRX_SEQUENCE,
                    source=SOURCE+" (Supplementary Table 1, SR001 representative)")
    unselected=graph.product("unselected_library","indexed small-RNA library before size selection",
                             MaterialKind.LIBRARY,state=library(),source=SOURCE)
    graph.react("pcr2","Install the indexed P7 arm with RP1 and SRX",StepKind.PCR,
                ((pcr1,InputRole.TEMPLATE),(rp1,InputRole.STRUCTURAL),
                 (srx,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(unselected,),
                geometry=PCRGeometry.exact(
                    pcr1_state, PCRPrimer("RP1", RP1, anneal_nt=RP1_ANNEAL_NT),
                    PCRPrimer("SRX (SR001)", SRX_SEQUENCE, anneal_nt=len(RTP)),
                    statement="Detailed-protocol Supplementary Table 1 prints every full SRX primer; SR001 is the representative geometry."))
    size=graph.input("size_selection","130–160-bp Pippin or 120–200-bp gel selection",
                     MaterialKind.SELECTION)
    final=library()
    selected=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                           source=SOURCE)
    offsize=graph.product("offsize","off-size library and some adapter dimers",MaterialKind.WASTE)
    graph.react("size_select","Select the small-RNA library size range",StepKind.SELECTION,
                ((unselected,InputRole.TEMPLATE),(size,InputRole.SELECTION)),(selected,offsize))
    graph.discard(offsize,"off-size molecules are excluded, although the overlapping adapter dimer is not fully removed")
    graph.sequence(selected,SequencingPlan(
        library=final,primers=(SEQ_PRIMERS[0], SEQ_PRIMERS[3]),
        acquired_roles=("Read 1", "Index 1 (i7)"),
        read_lengths={"Read 1":READ_LEN, "Index 1 (i7)":SRX_INDEX_NT},
        limitation="The detailed protocol establishes the SRX index read; it does not print the distinct TruSeq Small RNA Read 1 sequencing-primer oligo."))
    graph.declare_reviewed_scope(("single_cell_rna","rrna_mask","mask_rrna","ra3","rnl2",
      "ligate_3p","rtp","deadenylase","deadenylate","lambda_exonuclease",
      "digest_free_adapter","ra5","rnl1","ligate_5p","reverse_transcriptase",
      "reverse_transcribe","rp1","polymerase","pcr1","srx","pcr2","size_selection",
      "size_select","library"),SOURCE)
    return graph.require_publishable()
