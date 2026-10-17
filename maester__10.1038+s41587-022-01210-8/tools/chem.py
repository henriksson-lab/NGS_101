"""MAESTER 10x Genomics 3' v3 branch: source-owned declarations and reaction graph.

Source passages (published article, PMC9288977 author manuscript JATS text twin
``s41587-022-01210-8_PMC9288977.xml.txt``): Methods "10x Genomics" library
line 133, "MAESTER for 10x Genomics 3' v3" lines 157-171, MAESTER data
processing line 193.  Supplementary Table 3 (MOESM5 workbook, sheet
"MAESTER 10x", text lines 4-8, 11-77 primer mixes, 78-125 sample-index primers
and footnotes).  Supplementary Figure 3 (MOESM1 PDF p. 3, rendered) draws the
final 10x product.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "MAESTER — mitochondrial variant enrichment from scRNA-seq"
NOTES = "01_maester.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-022-01210-8">'
          'Miller et al., <i>Nature Biotechnology</i> (2022)</a>.')
SUMMARY = ("Part of the 10x Genomics 3′ v3 full-length cDNA is split into 12 multiplex PCRs: tiled "
           "mitochondrial primers carrying a common handle pair with an indexed P5 primer on the "
           "10x Read 1 end, keeping each cell barcode and UMI. A second PCR adds the indexed P7 end. "
           "The schematic shows the 10x branch; the Seq-Well/Drop-seq branch (biotinylated primers, "
           "streptavidin capture, Nextera handle) uses different oligos.")
CAVEAT = ("The mitochondrial primer panel (65 primers in 12 mixes, Supplementary Table 3) is drawn as one "
          "named target-site region, not a single amplicon. The Read 2 sequencing primer is not printed; "
          "the Methods state only that no custom primer is needed on the NovaSeq forward-strand workflow.")

# Supplementary Table 3 (MOESM5, sheet "MAESTER 10x"), verbatim 5'->3'.
PARTIAL_READ1 = "CTACACGACGCTCTTCCGATCT"      # line 4, optional cDNA amplification primer
P5_GENERIC = "AATGATACGGCGACCACCGAGATCTACAC"  # line 6
CI7P = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"    # line 7, Index 1 spike-in (reverse-complement workflow)
CI5P = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA"        # line 8, Index 2 spike-in (reverse-complement workflow)
HANDLE = "CACCCGAGAATTCCA"                   # common 5' part of every mitochondrial primer
ND1_4009 = HANDLE + "AACACCCTCACCACTACAATCT"  # line 11, PvG1218_MT-ND1_4009 (representative)
GOT_P5_TAIL = "TCTTTCCCTACACGACGCTC"          # GoT-P5-i5-BCXX after P5 + 8-nt i5 (lines 101-124)
GOT_P7_ARM = "GTGACTGGAGTTCCTTGG" + HANDLE    # GoT-P7-i7-BCXX after P7 + 8-nt i7 (lines 78-100)
GOT_P5 = P5_GENERIC + "N" * 8 + GOT_P5_TAIL
GOT_P7 = il.P7 + "N" * 8 + GOT_P7_ARM
if not GOT_P5.startswith(il.P5) or not GOT_P7_ARM.endswith(HANDLE):
    raise AssertionError("Supplementary Table 3 primer transcription changed")

# GoT-P5-i5 anneals to the 10x cDNA Read 1 end through its 3' 13 nt.
_OVERLAP = 13
if GOT_P5_TAIL[-_OVERLAP:] != PARTIAL_READ1[:_OVERLAP]:
    raise AssertionError("GoT-P5-i5 tail does not overlap PartialRead1")
READ1_ARM = GOT_P5_TAIL + PARTIAL_READ1[_OVERLAP:]

CELL = feature("maester_10x_cell_barcode", "cell_barcode", "whitelist",
               whitelist="10x Genomics 3′ v3 gel-bead barcode")
UMI = feature("maester_10x_umi", "umi", "random")
I5 = feature("maester_i5", "sample_index", "fixed",
             note="GoT-P5-i5-BC01–BC24, Supplementary Table 3; added in PCR1")
I7 = feature("maester_i7", "sample_index", "fixed",
             note="GoT-P7-i7-BC01–BC24, Supplementary Table 3; added in PCR2")


def _cell_core():
    # Supplementary Figure 3 draws a 16-position CB and 12-position UMI before poly(dT).
    return [Segment("cell barcode", "J" * 16, "cbc", placeholder=True, feature=CELL),
            Segment("UMI", "N" * 12, "umi", placeholder=True, feature=UMI),
            Segment("poly(dT)", "T" * 10, placeholder=True)]


MT_SITES = "mitochondrial target-primer sites"

def cdna_template():
    return Construct([Segment("10x partial Read 1", PARTIAL_READ1, "r1"), *_cell_core(),
                      Segment("transcript 3′ portion (antisense)", "X" * 20, placeholder=True),
                      Segment(MT_SITES, "X" * 22, placeholder=True),
                      Segment("upstream transcript and TSO end", "X" * 12, placeholder=True)],
                     name="10x 3′ v3 full-length cDNA")


def pcr1_product():
    return Construct([Segment("P5", P5_GENERIC, "p5"),
                      Segment("i5", "N" * 8, "cbc", placeholder=True, feature=I5),
                      Segment("Read 1 arm", READ1_ARM, "r1"), *_cell_core(),
                      Segment("transcript 3′ portion (antisense)", "X" * 20, placeholder=True),
                      Segment("target-primer site (reverse complement)", "X" * 22, placeholder=True),
                      Segment("handle reverse complement", revcomp(HANDLE), "r2")],
                     name="PCR1 mitochondrial amplicon")


FINAL_LIBRARY = Construct([
    Segment("P5", P5_GENERIC, "p5"),
    Segment("i5", "N" * 8, "cbc", placeholder=True, feature=I5),
    Segment("Read 1 arm", READ1_ARM, "r1"), *_cell_core(),
    Segment("mitochondrial transcript (antisense)", "X" * 20, placeholder=True),
    Segment("target-primer site (reverse complement)", "X" * 22, placeholder=True),
    Segment("GoT-P7 arm reverse complement", revcomp(GOT_P7_ARM), "r2"),
    Segment("i7 reverse complement", "N" * 8, "cbc", placeholder=True, feature=I7),
    Segment("P7 reverse complement", il.P7_RC, "p7")], name="MAESTER 10x library")

_T3 = "Miller et al. 2022, Supplementary Table 3"
SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "10x-Ci7P", CI7P, _T3,
              "Printed for the NovaSeq reverse-complement workflow; forward-strand runs use standard primers."),
    sp.custom("Index 2 (i5)", "10x-Ci5P", CI5P, _T3,
              "Printed for the NovaSeq reverse-complement workflow; forward-strand runs use standard primers."),
)
READ_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 8, "Index 2 (i5)": 8}
SEQUENCING_INTRO = ("NovaSeq SP 300-cycle kit: 28 cycles Read 1, 256 cycles Read 2 and 2 × 8 index cycles "
                    "(Methods). Read 2 starts at the mitochondrial primer's target sequence.")
FINAL_CAPTION = ("Supplementary Figure 3 final product, drawn here from the P5 end. The P7 arm is the "
                 "GoT-P7-i7 primer's small-RNA-type arm and handle, not the TruSeq Read 2 arm of the parent 10x library.")
SHOW_READ_CONTENT = False
SEQUENCING_UNAVAILABLE = ("The Read 2 primer that anneals to the GoT-P7 arm is not printed; the Methods say only that no "
                          "custom sequencing primers are required on the NovaSeq forward-strand workflow.")
READ_SPECS = (
    sp.read_spec("Read 1", 28, "16-nt 10x cell barcode, then 12-nt UMI (Methods line 193: CB + UMI)",
                 "extract"),
    sp.read_spec("Index 1 (i7)", 8, "GoT-P7-i7 sample index (demultiplex with its reverse complement)",
                 "demultiplex"),
    sp.read_spec("Index 2 (i5)", 8, "GoT-P5-i5 sample index (orientation depends on forward-strand vs reverse-complement workflow)",
                 "demultiplex"),
    sp.read_spec("Read 2", 256, "Mitochondrial primer target sequence, then mitochondrial transcript toward the poly(A); "
                 "the first 24 bp are trimmed before alignment (Methods line 193)", "include"),
)


def oligos():
    rows = (("PartialRead1 (optional cDNA amplification)", PARTIAL_READ1),
            ("PvG1218_MT-ND1_4009 (one of 65 mitochondrial primers)", ND1_4009),
            ("GoT-P5-i5-BCXX (PCR1; N = 8-nt i5)", GOT_P5),
            ("P5-generic (PCR2)", P5_GENERIC),
            ("GoT-P7-i7-BCXX (PCR2; N = 8-nt i7)", GOT_P7),
            ("10x-Ci7P (Index 1, reverse-complement workflow)", CI7P),
            ("10x-Ci5P (Index 2, reverse-complement workflow)", CI5P))
    return [f"<b>{name}</b> (Supplementary Table 3, 5′→3′): <code>{seq}</code>" for name, seq in rows]


def sections():
    return []


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    template = cdna_template()
    cdna = graph.input("cdna", "10x Genomics 3′ v3 cDNA, 20 ng per reaction", MaterialKind.BIOLOGICAL,
                       state=template, primary=True)
    mt = graph.input("mt_primers", "12 mitochondrial primer mixes (2–11 handle-tailed primers each, 65 total)",
                     MaterialKind.OLIGO, source=SOURCE + " (Supplementary Table 3)")
    p5 = graph.input("got_p5", "GoT-P5-i5-BCXX indexed P5 primer", MaterialKind.OLIGO,
                     state=Construct([Segment("GoT-P5-i5", GOT_P5)]), source=SOURCE + " (Supplementary Table 3)")
    kapa1 = graph.input("pcr1_polymerase", "KAPA HiFi HotStart ReadyMix; 6 cycles 98/65/72 °C",
                        MaterialKind.CATALYST)
    pcr1 = graph.product("pcr1", "12 PCR1 mitochondrial amplicon reactions", MaterialKind.MOLECULAR,
                         state=pcr1_product())
    graph.react("pcr1_multiplex", "PCR1: amplify mitochondrial transcripts in 12 primer-mix reactions", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL), (mt, InputRole.STRUCTURAL),
                 (kapa1, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.regional(
                    template, PCRPrimer("GoT-P5-i5-BCXX", GOT_P5, anneal_nt=_OVERLAP),
                    reverse_pool="12 mitochondrial primer mixes (65 primers)", reverse_region=MT_SITES,
                    statement="GoT-P5-i5 overlaps the 10x Read 1 end by 13 nt; each mitochondrial primer anneals at its own transcript site."))
    ampure1 = graph.input("ampure_1x", "empirical pooling ratios, then 1× AMPure XP; elute 20 µl H2O",
                          MaterialKind.SELECTION)
    pooled = graph.product("pooled_pcr1", "pooled, purified PCR1 amplicons", MaterialKind.MOLECULAR,
                           state=pcr1.state)
    waste1 = graph.product("pcr1_primers", "free primers and short products", MaterialKind.WASTE)
    graph.react("pool_cleanup", "Pool PCR1 reactions and remove primers", StepKind.CLEANUP,
                ((pcr1, InputRole.TEMPLATE), (ampure1, InputRole.SELECTION)), (pooled, waste1))
    graph.discard(waste1, "is removed by 1× AMPure XP")
    p5g = graph.input("p5_generic", "P5-generic primer", MaterialKind.OLIGO,
                      state=Construct([Segment("P5-generic", P5_GENERIC, "p5")]))
    p7 = graph.input("got_p7", "GoT-P7-i7-BCXX indexed P7 primer", MaterialKind.OLIGO,
                     state=Construct([Segment("GoT-P7-i7", GOT_P7)]))
    kapa2 = graph.input("pcr2_polymerase", "KAPA HiFi HotStart ReadyMix; 6 cycles 98/60/72 °C",
                        MaterialKind.CATALYST)
    library = graph.product("pcr2", "indexed MAESTER library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("pcr2_index", "PCR2: add the indexed P7 end", StepKind.PCR,
                ((pooled, InputRole.TEMPLATE), (p5g, InputRole.STRUCTURAL), (p7, InputRole.STRUCTURAL),
                 (kapa2, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    pooled.state, PCRPrimer("P5-generic", P5_GENERIC),
                    PCRPrimer("GoT-P7-i7-BCXX", GOT_P7, anneal_nt=len(HANDLE)),
                    statement="Supplementary Table 3 and Methods: P5-generic plus GoT-P7-i7 (5 µM mix), KAPA, 6 cycles."))
    ampure2 = graph.input("ampure_08x", "0.8× AMPure XP; elute 20 µl TE", MaterialKind.SELECTION)
    clean = graph.product("library", "purified MAESTER library (300–1,500 bp)", MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    waste2 = graph.product("pcr2_waste", "primers and short products", MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify the library", StepKind.CLEANUP,
                ((library, InputRole.TEMPLATE), (ampure2, InputRole.SELECTION)), (clean, waste2))
    graph.discard(waste2, "is removed by 0.8× AMPure XP")
    graph.sequence(clean, SequencingPlan(library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
                                         read_lengths=READ_LENGTHS))
    graph.annotate(SourceDiscrepancy(
        "library", "Published primer names differ",
        (SourceClaim("Methods, 10x PCR1/PCR2 (PMC9288977 lines 161, 169)",
                     "PCR1 uses GoT-P5-i5-BCXX; PCR2 uses P5-generic and XV-P7-i7-BCXX"),
         SourceClaim("Supplementary Figure 3 (MOESM1 p. 3)",
                     "PCR1 drawn with XV-P5-i5-BCXX; PCR2 with XV-P7-i7-BCXX"),
         SourceClaim("Supplementary Table 3 (MOESM5 lines 78–125)",
                     "only GoT-P7-i7-BC01–24 and GoT-P5-i5-BC01–24 are printed; PCR2 footnote uses GoT-P7-i7-BCXX")),
        "The model uses the Supplementary Table 3 names and sequences.",
        "The final product drawn in Supplementary Figure 3 matches the Table 3 GoT-P5-i5 and GoT-P7-i7 sequences.",
        "No separate XV-series sequence is printed; the XV names are treated as aliases, not distinct oligos."))
    graph.declare_reviewed_scope((
        "cdna", "mt_primers", "got_p5", "pcr1_polymerase", "pcr1_multiplex", "ampure_1x", "pool_cleanup",
        "p5_generic", "got_p7", "pcr2_polymerase", "pcr2_index", "ampure_08x", "final_cleanup", "library"),
        SOURCE)
    return graph.require_publishable()


SPEC = sys.modules[__name__]
