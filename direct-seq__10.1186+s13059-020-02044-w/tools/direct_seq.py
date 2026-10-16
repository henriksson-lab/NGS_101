"""Direct-seq programmed-guide capture and nested-PCR model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow, duplex_rows,
                      feature, revcomp)
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import illumina as il
import seqprimers as sp


TITLE = "Direct-seq — poly(dT) capture of programmed guide RNA"
NOTES = "01_direct-seq.html"
SOURCE = '<a href="https://doi.org/10.1186/s13059-020-02044-w">Song et al., <i>Genome Biology</i> (2020)</a>.'
SUMMARY = ("A 30-nt A/G capture tract is encoded in the sgRNA scaffold. A normal "
           "10x 3′ v3 bead primer therefore copies guide RNA with the same cell barcode "
           "and UMI as mRNA; exact tRNA-directed nested PCR then makes the guide library.")
CAVEAT = ("The programmed transcript uses the published Tail-8A8G scaffold and human "
          "Gln tRNA. G residues intentionally interrupt the capture tract, so the poly(dT) "
          "drawing marks that duplex as mismatch-tolerant rather than pretending it is a "
          "perfect 30-bp poly(A):poly(dT) helix. The 10x bead handle follows the "
          "paper-cited CG000184 Rev A oligonucleotide appendix; cell-barcode and UMI "
          "members are variable.")

CAPTURE = "AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA"
BEAD_HANDLE = "CTACACGACGCTCTTCCGATCT"  # 10x CG000184 Rev A, oligo appendix.
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


# PCR templates are represented independently of the final sequencing construct.  The
# first template carries the 10x pre-amplification handle on which P5_Read1 lands.  Its
# product carries the complete tRNA_Read2 sequence, which is the landing site for the
# indexed reverse primer in the following PCR.  Keeping these reaction states distinct
# avoids silently treating an intermediate 10x handle as part of the biological insert.
PREAMPLIFIED_GUIDE = Construct([
    Segment("P5_Read1 annealing region", P5_READ1[-20:]),
    Segment("cell barcode", "B" * 16, "cbc", placeholder=True),
    Segment("UMI", "U" * 12, "umi", placeholder=True),
    Segment("captured guide cDNA", "X" * 80, placeholder=True),
    Segment("tRNA primer-site cDNA", revcomp(TRNA_BIND)),
], "pre-amplified Direct-seq guide cDNA")

FIRST_GUIDE_AMPLICON = Construct([
    Segment("P5_Read1 product end", P5_READ1),
    Segment("cell barcode", "B" * 16, "cbc", placeholder=True),
    Segment("UMI", "U" * 12, "umi", placeholder=True),
    Segment("captured guide cDNA", "X" * 80, placeholder=True),
    Segment("tRNA_Read2 product end", revcomp(TRNA_READ2)),
], "first Direct-seq guide-PCR amplicon")

FIRST_GUIDE_PCR = PCRGeometry.exact(
    PREAMPLIFIED_GUIDE, PCRPrimer("P5_Read1", P5_READ1),
    PCRPrimer("tRNA_Read2", TRNA_READ2),
    statement="Supplementary Table S3d gives both first guide-enrichment primers.")
INDEX_GUIDE_PCR = PCRGeometry.exact(
    FIRST_GUIDE_AMPLICON, PCRPrimer("P5_Read1", P5_READ1),
    PCRPrimer("P7_Read2_Index1", P7_READ2_INDEX1),
    statement="Supplementary Table S3d gives the indexed terminal-PCR primers.")


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
    primer = [Segment("10x bead handle", BEAD_HANDLE),
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
        Segment("10x cDNA handle", BEAD_HANDLE),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True),
        Segment("UMI", "U" * 12, "umi", placeholder=True),
        Segment("poly(dT)-derived tract", "T" * 30),
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


# Source-owned 10x v3 / Table S3d architecture. The GEX and guide branches
# share the TruSeq barcode end; tRNA_Read2 supplies the guide branch's other end.
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x Chromium 3-prime v3")
UMI_FEATURE = feature("umi", "umi", "random")
INDEX_FEATURE = feature("sample_index_i7", "sample_index", "unknown")


def require_guide_primer_sites(library):
    """Reject an endpoint incompatible with the independently printed PCR pair."""
    return PCRGeometry.exact(library, PCRPrimer("Table S3d P5_Read1", P5_READ1, anneal_nt=len(P5_READ1)),
                             PCRPrimer("Table S3d P7_Read2_Index1", P7_READ2_INDEX1),
                             statement="Direct-seq Table S3d source primers must land on the declared guide endpoint.")


def _indexed_library(payload, name):
    library = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
        Segment("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
        *payload,
        Segment("Read 2 primer-site complement", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True, feature=INDEX_FEATURE),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name)
    problems = sp.verify(library, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS))
    if problems:
        raise ValueError("Direct-seq endpoint: " + "; ".join(problems))
    require_guide_primer_sites(library)
    return library


RNA_LIB = _indexed_library([
    Segment("poly(dT)-derived tract", "T" * 30),
    Segment("endogenous cDNA", "X" * 40, placeholder=True),
], "Direct-seq gene-expression library")
GUIDE_LIB = _indexed_library([
    Segment("poly(dT)-derived tract", "T" * 30),
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
        f"<b>10x CG000184 Rev A gel-bead primer</b> — 5′-{BEAD_HANDLE}-N16-N12-{'T' * 30}-3′",
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


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input(
        "cells", "cells expressing programmed guide RNA and endogenous mRNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("Direct-seq cells", (role("programmed guide RNA"),
                                             role("poly(A) mRNA"))),
        primary=True)
    gems = g.input("gems", "10x 3-prime v3 GEM partitioning condition",
                   MaterialKind.CONDITION)
    occupied = g.product("occupied", "occupied barcoded GEM reactions",
                         MaterialKind.MOLECULAR, state=cells.state)
    failed = g.product("failed", "empty or unresolved GEM products",
                       MaterialKind.WASTE)
    g.react("partition", "Partition cells with barcoded gel beads",
            StepKind.PARTITION,
            ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
            (occupied, failed), changes_molecule=False)
    g.discard(failed, "does not retain a resolvable cell barcode")

    bead = g.input("bead_primer", "10x cell-barcode/UMI oligo(dT) bead primer",
                   MaterialKind.OLIGO,
                   state=Construct([Segment("10x Read 1 handle", BEAD_HANDLE),
                                    Segment("cell barcode", "B" * 16, placeholder=True),
                                    Segment("UMI", "U" * 12, placeholder=True),
                                    Segment("poly(dT)30", "T" * 30)],
                                   "10x 3′ v3 gel-bead primer"),
                   source="10x CG000184 Rev A, oligonucleotide appendix")
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    shared = g.product("shared_cdna", "cell-barcoded mRNA and guide cDNA",
                       MaterialKind.MOLECULAR,
                       state=molecule("shared Direct-seq cDNA",
                                      (role("cell barcode"), role("UMI"),
                                       role("mRNA cDNA"), role("guide cDNA"))))
    g.react("reverse_transcribe",
            "Capture mRNA and the programmed 8A8G guide with the same bead primer",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (shared,))

    rna_branch = g.product("rna_branch", "gene-expression cDNA branch",
                           MaterialKind.MOLECULAR, state=shared.state)
    guide_branch = g.product("guide_branch", "guide-enrichment cDNA branch",
                             MaterialKind.MOLECULAR, state=PREAMPLIFIED_GUIDE)
    g.react("split", "Split pre-amplified cDNA for gene and guide libraries",
            StepKind.SPLIT, ((shared, InputRole.TEMPLATE),),
            (rna_branch, guide_branch), changes_molecule=False)

    reviewed = ["cells", "gems", "partition", "bead_primer", "rt",
                "reverse_transcribe", "split"]
    finish_pcr_branch(
        g, "rna", rna_branch, RNA_LIB, SEQ_PRIMERS, reviewed,
        label="Direct-seq gene expression",
        read_lengths=READ_LENGTHS["Gene-expression library"],
        read_uses=READ_USES["Gene-expression library"], inferred_primers=True)

    p5 = g.input("p5_read1", "P5_Read1 forward primer", MaterialKind.OLIGO,
                 state=P5_READ1)
    trna = g.input("trna_read2", "tRNA_Read2 reverse primer", MaterialKind.OLIGO,
                   state=TRNA_READ2)
    polymerase_1 = g.input("guide_polymerase_1", "guide-enrichment PCR polymerase",
                           MaterialKind.CATALYST)
    first = g.product("first_guide_amplicon", "first nested guide-PCR amplicon",
                      MaterialKind.MOLECULAR, state=FIRST_GUIDE_AMPLICON)
    g.react("guide_pcr_1", "Enrich guide cDNA with P5_Read1 and tRNA_Read2",
            StepKind.PCR,
            ((guide_branch, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (trna, InputRole.STRUCTURAL), (polymerase_1, InputRole.CATALYST)),
            (first,), geometry=FIRST_GUIDE_PCR)

    indexed = g.input("p7_indexed", "P7_Read2 indexed reverse primer",
                      MaterialKind.OLIGO, state=P7_READ2_INDEX1)
    polymerase_2 = g.input("guide_polymerase_2", "terminal guide-PCR polymerase",
                           MaterialKind.CATALYST)
    guide_library = g.product("guide_library", "indexed guide sequencing library",
                              MaterialKind.LIBRARY, state=GUIDE_LIB)
    g.react("guide_pcr_2", "Add the indexed P7 end to the enriched guide library",
            StepKind.PCR,
            ((first, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (indexed, InputRole.STRUCTURAL), (polymerase_2, InputRole.CATALYST)),
            (guide_library,), geometry=INDEX_GUIDE_PCR)
    g.sequence(guide_library, SequencingPlan(
        library=GUIDE_LIB, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths=READ_LENGTHS["Enriched guide library"],
        read_uses=READ_USES["Enriched guide library"]))
    reviewed.extend(("p5_read1", "trna_read2", "guide_polymerase_1",
                     "guide_pcr_1", "p7_indexed", "guide_polymerase_2",
                     "guide_pcr_2", "guide_library"))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
