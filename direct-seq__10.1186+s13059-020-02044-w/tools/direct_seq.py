"""Direct-seq programmed-guide capture and nested-PCR model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import (MolecularState, Scene, Segment, Workflow, duplex_rows,
                      feature, revcomp)
from tenx_rna import SEQ_PRIMERS, captured_rna_library, transcript_library


TITLE = "Direct-seq — poly(dT) capture of programmed guide RNA"
NOTES = "01_direct-seq.html"
SOURCE = '<a href="https://doi.org/10.1186/s13059-020-02044-w">Song et al., <i>Genome Biology</i> (2020)</a>.'
SUMMARY = ("A 30-nt A/G capture tract is encoded in the sgRNA scaffold. A normal "
           "10x 3′ v3 bead primer therefore copies guide RNA with the same cell barcode "
           "and UMI as mRNA; exact tRNA-directed nested PCR then makes the guide library.")
CAVEAT = ("The programmed transcript uses the published Tail-8A8G scaffold and human "
          "Gln tRNA. G residues intentionally interrupt the capture tract, so the poly(dT) "
          "drawing marks that duplex as mismatch-tolerant rather than pretending it is a "
          "perfect 30-bp poly(A):poly(dT) helix.")

CAPTURE = "AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA"
HUMAN_TRNA_GLN = "GGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCT"
TRNA_BIND = "CTCTGGACTCTGAATCCAGCGAT"
if HUMAN_TRNA_GLN.find(TRNA_BIND) != 24:
    raise ValueError("Direct-seq tRNA primer site moved from its source-defined position")
TRNA_DOWNSTREAM = HUMAN_TRNA_GLN[24 + len(TRNA_BIND):]

# Table S2/S3 Tail-8A8G programmed scaffold, normalized to DNA capitals.
SCAFFOLD_BEFORE_CAPTURE = ("GTTTGAGAGCTATGCTGGAAACAGCATAGCAAGTTCAAATAAGGCTAGTCCGTTATCA"
                           "ACTTGAAAAAGTGGCACCGAGTCGGTGCTCGG")
POL_III_TERMINATOR = "T" * 7
GUIDE = feature("direct_guide", "guide_barcode", "whitelist",
                whitelist="Direct-seq screen guide set")

# Supplementary Table S3d.
P5_READ1 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
TRNA_READ2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" + TRNA_BIND
P7_READ2_INDEX1 = "CAAGCAGAAGACGGCATACGAGATACTCGCTA" + "GTGACTGGAGTTCAGACGTGT"
P7_READ2_INDEX2 = "CAAGCAGAAGACGGCATACGAGATGGAGCTAC" + "GTGACTGGAGTTCAGACGTGT"


def _programmed_transcript_scene() -> Scene:
    sc = Scene()
    sc.strand("guide RNA", [
        Segment("human Gln tRNA", HUMAN_TRNA_GLN),
        Segment("guide spacer", "G" * 20, "cbc", placeholder=True, feature=GUIDE),
        Segment("programmed sgRNA scaffold", SCAFFOLD_BEFORE_CAPTURE),
        Segment("8A8G capture tract", CAPTURE),
        Segment("Pol III terminator", POL_III_TERMINATOR),
    ], label="programmed Pol III guide transcript")
    sc.labels("guide RNA")
    return sc


def _capture_scene() -> Scene:
    transcript = [Segment("guide body", "X" * 28, placeholder=True),
                  Segment("8A8G capture tract", CAPTURE)]
    primer = [Segment("10x bead handle", "H" * 12, placeholder=True),
              Segment("cell barcode", "B" * 16, "cbc", placeholder=True),
              Segment("UMI", "U" * 12, "umi", placeholder=True),
              Segment("poly(dT)30", "T" * 30)]
    sc = Scene(); sc.strand("guide RNA", transcript, label="programmed guide RNA")
    sc.anneal("bead primer", primer, to="guide RNA",
              pair=("poly(dT)30", "8A8G capture tract"),
              unpaired=("10x bead handle", "cell barcode", "UMI", "poly(dT)30",
                        "8A8G capture tract"),
              label="10x 3′ v3 gel-bead primer")
    sc.arrow("bead primer", "reverse transcription across programmed guide")
    sc.note("guide RNA", "three intentional G interruptions make this a mismatch-tolerant capture duplex")
    sc.labels("bead primer")
    return sc


def _barcoded_guide_cdna_scene() -> Scene:
    sc = Scene(); sc.strand("guide cDNA", [
        Segment("10x cDNA handle", "H" * 12, placeholder=True),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True),
        Segment("UMI", "U" * 12, "umi", placeholder=True),
        Segment("poly(dT)-derived tract", "T" * 20),
        Segment("capture-tract cDNA", revcomp(CAPTURE)),
        Segment("scaffold-derived cDNA", revcomp(SCAFFOLD_BEFORE_CAPTURE)),
        Segment("guide spacer cDNA", "F" * 20, "cbc", placeholder=True, feature=GUIDE),
        Segment("tRNA-derived cDNA", revcomp(HUMAN_TRNA_GLN)),
    ], label="cell-barcoded guide cDNA")
    sc.labels("guide cDNA")
    return sc


def _nested_pcr_scene() -> Scene:
    target = [Segment("guide cDNA upstream", "X" * 30, placeholder=True),
              Segment("tRNA primer-site complement", revcomp(TRNA_BIND))]
    primer = [Segment("Read-2 adapter tail", TRNA_READ2[:-len(TRNA_BIND)], "r2"),
              Segment("tRNA-specific 3′ end", TRNA_BIND)]
    sc = Scene(); sc.strand("guide cDNA", target, label="pre-amplified 10x cDNA")
    sc.anneal("tRNA_Read2", primer, to="guide cDNA",
              pair=("tRNA-specific 3′ end", "tRNA primer-site complement"),
              unpaired=("Read-2 adapter tail",), label="first nested-PCR reverse primer")
    sc.arrow("tRNA_Read2", "guide enrichment PCR")
    sc.labels("tRNA_Read2")
    return sc


RNA_LIB = transcript_library("Direct-seq gene-expression library")
GUIDE_LIB = captured_rna_library([
    Segment("poly(dT)-derived tract", "T" * 20),
    Segment("capture-tract cDNA", revcomp(CAPTURE)),
    Segment("scaffold-derived cDNA", revcomp(SCAFFOLD_BEFORE_CAPTURE)),
    Segment("guide spacer cDNA", "F" * 20, "cbc", placeholder=True, feature=GUIDE),
    Segment("post-primer tRNA cDNA", revcomp(TRNA_DOWNSTREAM)),
    Segment("tRNA primer-site cDNA", revcomp(TRNA_BIND)),
], "Direct-seq enriched guide library")

FINAL_LIBRARIES = (
    ("Gene-expression library", RNA_LIB, SEQ_PRIMERS,
     "Standard 10x 3′ v3 transcript library.",
     "Read 1 records cell barcode and UMI; Read 2 enters endogenous cDNA."),
    ("Enriched guide library", GUIDE_LIB, SEQ_PRIMERS,
     "The exact programmed-guide cDNA is retained between the common 10x identifier end and the tRNA-directed Read-2 end.",
     "Read 1 records the same cell barcode and UMI. Read 2 traverses the tRNA primer site and remaining tRNA sequence before identifying the guide spacer."),
)
READ_LENGTHS = {"Gene-expression library": {"Read 1": 28},
                "Enriched guide library": {"Read 1": 28, "Read 2": 90}}
READ_USES = {"Gene-expression library": {"Read 1": "extract", "Read 2": "include"},
             "Enriched guide library": {"Read 1": "extract", "Read 2": "include"}}


def oligos():
    return (
        f"<b>8A8G capture tract</b> — 5′-{CAPTURE}-3′",
        f"<b>Human Gln tRNA</b> — 5′-{HUMAN_TRNA_GLN}-3′",
        f"<b>P5_Read1</b> — 5′-{P5_READ1}-3′",
        f"<b>tRNA_Read2</b> — 5′-{TRNA_READ2}-3′",
        f"<b>P7_Read2_Index1</b> — 5′-{P7_READ2_INDEX1}-3′",
        f"<b>P7_Read2_Index2</b> — 5′-{P7_READ2_INDEX2}-3′",
    )


def workflow() -> Workflow:
    wf = Workflow(MolecularState("Programmed guide transcript",
                                 tuple(_programmed_transcript_scene().rows())))
    wf.react("Capture guide RNA and mRNA with the same 10x poly(dT) bead primer",
             _capture_scene().rows(), name="Guide RNA annealed to bead primer",
             note="The guide's encoded 8A8G tract replaces a conventional poly(A) tail; no guide-specific RT primer is added.")
    wf.react("Reverse-transcribe and pre-amplify cell-barcoded guide cDNA",
             _barcoded_guide_cdna_scene().rows(), name="Shared pre-amplified cDNA",
             note="The captured guide receives the ordinary 16-nt cell barcode and 12-nt UMI in the GEM.")
    wf.react("Split the cDNA and enrich guide molecules with P5_Read1 and tRNA_Read2",
             _nested_pcr_scene().rows(), name="First nested guide-PCR product",
             note="The 3′ 23 nt of tRNA_Read2 match bases 25–47 of the published human Gln tRNA; its 5′ tail installs the Read-2 site.")
    wf.react("Repeat PCR with P5_Read1 and indexed P7_Read2",
             duplex_rows(GUIDE_LIB, label="guide-index library"),
             name="Final guide-index library",
             note="The third PCR adds P7 and one of the exact eight-base sample indexes. The mRNA half follows the normal 10x 3′ v3 library workflow.")
    return wf
