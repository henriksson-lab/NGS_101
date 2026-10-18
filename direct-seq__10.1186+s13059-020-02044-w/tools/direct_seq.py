"""Direct-seq programmed-guide capture and nested-PCR model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow, duplex_rows,
                      feature, revcomp)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import molecule, role

SOURCE_OLIGO_EXCEPTIONS={"bead_primer":"Inherited exact scaffold from paper-cited CG000184 Rev A Appendix (source_reviews prior independent reads): sources/10x-chromium-single-cell-3-v3-crispr-screening/scg_CG000184_ChromiumSingleCellSingleCell3v3_FeatureBarcodingtechnology_CRISPR_RevA.pdf.txt lines6215–6235; original is cached in that vendor protocol bundle, not this article bundle."}
import illumina as il
import seqprimers as sp


TITLE = "Direct-seq — poly(dT) capture of programmed guide RNA"
NOTES = "01_direct-seq.html"
SOURCE = '<a href="https://doi.org/10.1186/s13059-020-02044-w">Song et al., <i>Genome Biology</i> (2020)</a>.'
METHODS = ("Song et al. 2020 Methods, 'scRNA-seq on the 10 × 3′ RNA-seq platform' "
           "(PMC7278172 JATS lines 132–144)")
SUMMARY = ("A 30-nt A/G capture tract is encoded in the sgRNA scaffold. A normal "
           "10x 3′ v3 bead primer therefore copies guide RNA with the same cell barcode "
           "and UMI as mRNA; exact tRNA-directed nested PCR then makes the guide library.")
CAVEAT = ("The programmed transcript uses the published Tail-8A8G scaffold and human "
          "Gln tRNA. G residues intentionally interrupt the capture tract, so the poly(dT) "
          "drawing marks that duplex as mismatch-tolerant rather than pretending it is a "
          "perfect 30-bp poly(A):poly(dT) helix. The 10x bead handle follows the "
          "paper-cited CG000184 Rev A oligonucleotide appendix; cell-barcode and UMI "
          "members are variable. INFERRED — outer preamplification ends remain unassigned; the known bead handle is retained without extending it to an assumed PCR-primer landing site, and actual mixed incoming PCR states remain distinct from candidate final libraries.")

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


CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x Chromium 3-prime v3")
UMI_FEATURE = feature("umi", "umi", "random")
INDEX_FEATURE = feature("sample_index_i7", "sample_index", "unknown")


# PCR templates are represented independently of the final sequencing construct.  The
# first template retains the bead-derived handle but does not invent its unknown
# outer pre-amplification end. The first enrichment PCR candidate
# product carries the complete tRNA_Read2 sequence, which is the landing site for the
# indexed reverse primer in the following PCR.  Keeping these reaction states distinct
# avoids silently treating an intermediate 10x handle as part of the biological insert.
PREAMPLIFIED_GUIDE = Construct([
    role("native preamplification 5′ extent unassigned",inferred=True),
    Segment("bead-derived Read1 handle",BEAD_HANDLE),
    Segment("cell barcode", "B" * 16, "cbc", placeholder=True, length_bp=16, feature=CELL_BARCODE),
    Segment("UMI", "U" * 12, "umi", placeholder=True, length_bp=12, feature=UMI_FEATURE),
    Segment("poly(dT)-derived tract", "T"*30),
    Segment("scaffold-derived cDNA",revcomp(SCAFFOLD_BEFORE_CAPTURE)),
    Segment("guide spacer cDNA","F"*20,placeholder=True,length_bp=20,feature=GUIDE),
    Segment("post-primer tRNA cDNA",revcomp(TRNA_DOWNSTREAM)),
    Segment("tRNA primer-site cDNA", revcomp(TRNA_BIND)),
    Segment("remaining tRNA cDNA",revcomp(HUMAN_TRNA_GLN[:24])),
    role("opposite preamplification end/extent unassigned",inferred=True),
], "pre-amplified Direct-seq guide cDNA")

def retained_bead_handle_suffix():
    """Local sequence relation only, not a complete native PCR landing claim.

    Table S3d P5_Read1 ends within the independently reviewed CG000184
    bead handle. The remaining handle must be copied before the cell barcode.
    """
    overlaps=[n for n in range(1,min(len(P5_READ1),len(BEAD_HANDLE))+1)
              if P5_READ1.endswith(BEAD_HANDLE[:n])]
    if not overlaps:
        raise ValueError("source P5_Read1 has no terminal overlap with the reviewed bead handle")
    return BEAD_HANDLE[max(overlaps):]


def require_retained_bead_handle(library):
    prefix=[]
    for segment in library:
        if segment.feature == CELL_BARCODE:
            break
        prefix.append(segment.top)
    else:
        raise ValueError("guide product lacks its cell-barcode boundary")
    if "".join(prefix) != P5_READ1+retained_bead_handle_suffix():
        raise ValueError("guide product must retain the copied bead-handle suffix before its cell barcode")
    return library


FIRST_GUIDE_AMPLICON = Construct([
    Segment("P5_Read1 product end", P5_READ1),
    Segment("retained bead Read1 handle", retained_bead_handle_suffix()),
    Segment("cell barcode", "B" * 16, "cbc", placeholder=True, length_bp=16, feature=CELL_BARCODE),
    Segment("UMI", "U" * 12, "umi", placeholder=True, length_bp=12, feature=UMI_FEATURE),
    Segment("poly(dT)-derived tract", "T"*30),
    Segment("scaffold-derived cDNA",revcomp(SCAFFOLD_BEFORE_CAPTURE)),
    Segment("guide spacer cDNA","F"*20,placeholder=True,length_bp=20,feature=GUIDE),
    Segment("post-primer tRNA cDNA",revcomp(TRNA_DOWNSTREAM)),
    Segment("tRNA_Read2 product end", revcomp(TRNA_READ2)),
], "first Direct-seq guide-PCR amplicon")
require_retained_bead_handle(FIRST_GUIDE_AMPLICON)

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
              Segment("cell barcode", "B" * 16, "cbc", placeholder=True, length_bp=16, feature=CELL_BARCODE),
              Segment("UMI", "U" * 12, "umi", placeholder=True, length_bp=12, feature=UMI_FEATURE),
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
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True, length_bp=16, feature=CELL_BARCODE),
        Segment("UMI", "U" * 12, "umi", placeholder=True, length_bp=12, feature=UMI_FEATURE),
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
    require_retained_bead_handle(library)
    return library


RNA_LIB = Construct([role("10x v3 transcript library; native preparation and adapter/index geometry unassigned", inferred=True)],
                    name="Direct-seq manufacturer-prepared expression library (bounded)")

GUIDE_LIB = _indexed_library([
    Segment("poly(dT)-derived tract", "T" * 30),
    Segment("scaffold-derived cDNA", revcomp(SCAFFOLD_BEFORE_CAPTURE)),
    Segment("guide spacer cDNA", "F" * 20, "cbc", placeholder=True, feature=GUIDE),
    Segment("post-primer tRNA cDNA", revcomp(TRNA_DOWNSTREAM)),
    Segment("tRNA primer-site cDNA", revcomp(TRNA_BIND)),
], "Direct-seq enriched guide library")

FINAL_LIBRARIES = (
    ("Gene-expression library", RNA_LIB, (),
     "Standard 10x 3′ v3 transcript library.",
     "Read 1 records cell barcode and UMI; Read 2 enters endogenous cDNA."),
    ("Enriched guide library", GUIDE_LIB, SEQ_PRIMERS,
     "The exact programmed-guide cDNA is retained between the common 10x identifier end and the tRNA-directed Read-2 end.",
     "Read 1 records the same cell barcode and UMI. Read 2 traverses the tRNA primer site and remaining tRNA sequence before identifying the guide spacer."),
)
# The inspected article, supplements and Methods print no sequencing run lengths;
# acquired cycles are therefore not declared (identifier spans come from features).
READ_LENGTHS = {}
READ_USES = {"Gene-expression library": {"Read 1": "extract", "Read 2": "include"},
             "Enriched guide library": {"Read 1": "extract", "Read 2": "include"}}

READ_SPECS={"Gene-expression library":(
    sp.read_spec("Read 1",None,"native primer and acquired cycles unassigned; 10x v3 identifier readout not reconstructed","undetermined"),sp.read_spec("Read 2",None,"transcript sequence in the represented v3 route; native primer unassigned","include"),)}
SEQUENCING_UNAVAILABLE={"Gene-expression library":"CG000184 manufacturer library construction is not reconstructed here; native primer/index geometry and acquired cycles remain unassigned."}


def independent_populations(label, *groups):
    sc=Scene()
    for index,parts in enumerate(groups):
        sc.strand(f"population {index+1}",parts,label=f"{label}: population {index+1}")
    return sc


def preamplified_mixture():
    sc=Scene.duplex(tuple(PREAMPLIFIED_GUIDE),unpaired=tuple(s.name+"'" for s in PREAMPLIFIED_GUIDE if s.is_role_token()))
    sc.strand("other preamplified cDNA population",[role("preamplified endogenous cDNA; native ends unassigned")],
              label="same tube: endogenous transcript cDNA population")
    return sc


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
             note="The 3rd PCR adds P7 and one of the exact eight-base sample indexes. The mRNA half follows the normal 10x 3′ v3 library workflow.")
    return wf


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cells = g.input(
        "cells", "cells expressing programmed guide RNA and endogenous mRNA",
        MaterialKind.BIOLOGICAL,
        state=independent_populations("Direct-seq cells",(role("programmed guide RNA"),),(role("poly(A) mRNA"),)),
        primary=True)
    gems = g.input("gems", "10x 3-prime v3 GEM partitioning condition",
                   MaterialKind.CONDITION)
    occupied = g.product("occupied", "occupied barcoded GEM reactions",
                         MaterialKind.MOLECULAR, state=cells.state)
    failed = g.product("failed", "empty or unresolved GEM products",
                       MaterialKind.MOLECULAR)
    g.react("partition", "Partition cells with barcoded gel beads",
            StepKind.PARTITION,
            ((cells, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
            (occupied, failed), changes_molecule=False)
    g.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,"Empty/ambiguous GEM outcomes are computational assignments, not a demonstrated physical discard."))

    bead = g.input("bead_primer", "10x cell-barcode/UMI oligo(dT) bead primer",
                   MaterialKind.OLIGO,
                   state=Construct([Segment("10x Read 1 handle", BEAD_HANDLE),
                                    Segment("cell barcode", "B" * 16, placeholder=True, length_bp=16, feature=CELL_BARCODE),
                                    Segment("UMI", "U" * 12, placeholder=True, length_bp=12, feature=UMI_FEATURE),
                                    Segment("poly(dT)30", "T" * 30)],
                                   "10x 3′ v3 gel-bead primer"),
                   source="10x CG000184 Rev A, oligonucleotide appendix")
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    shared = g.product("shared_cdna", "cell-barcoded mRNA and guide cDNA",
                       MaterialKind.MOLECULAR,
                       state=independent_populations("first-strand cDNA mixture",
                         tuple(_barcoded_guide_cdna_scene().strands["guide cDNA"].segs),
                         (role("cell barcode",molecular_feature=CELL_BARCODE),role("UMI",molecular_feature=UMI_FEATURE),role("endogenous first-strand cDNA"))))
    g.react("reverse_transcribe",
            "Capture mRNA and the programmed 8A8G guide with the same bead primer",
            StepKind.REVERSE_TRANSCRIPTION,
            ((occupied, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (shared,))

    tenx_src = "10x CG000184 (cited by Methods; cDNA primer species and bases not declared in this model)"
    tenx_primers = g.input("tenx_cdna_primers", "10x cDNA Primers, Read 1-handle species",
                           MaterialKind.OLIGO, inferred=True, source=tenx_src)
    tenx_primers_2 = g.input("tenx_cdna_primers_2", "10x cDNA Primers, opposite-end species",
                             MaterialKind.OLIGO, inferred=True, source=tenx_src)
    preamp_pol = g.input("preamp_polymerase", "10x cDNA amplification master mix",
                         MaterialKind.CATALYST)
    preamp = g.product("preamp_cdna", "pre-amplified cDNA (1st PCR)",
                       MaterialKind.MOLECULAR, state=preamplified_mixture())
    g.react("preamplify", "INFERRED outer end boundaries — pre-amplify cDNA following CG000184 (the 1st PCR)", StepKind.PCR,
            ((shared, InputRole.TEMPLATE), (tenx_primers, InputRole.STRUCTURAL),
             (tenx_primers_2, InputRole.STRUCTURAL), (preamp_pol, InputRole.CATALYST)), (preamp,),
            source=METHODS,
            geometry=PCRGeometry.unavailable(
                shared.state,
                reason="Methods follow CG000184 up to cDNA pre-amplification; its primers are not declared here."))
    aliquot_a = g.product("aliquot_a", "40 µl cDNA aliquot for size selection",
                          MaterialKind.MOLECULAR, state=preamp.state)
    aliquot_b = g.product("aliquot_b", "40 µl cDNA aliquot for gRNA library II",
                          MaterialKind.MOLECULAR, state=preamp.state)
    g.react("split", "Divide pre-amplified cDNA into two 40 µl aliquots",
            StepKind.SPLIT, ((preamp, InputRole.TEMPLATE),),
            (aliquot_a, aliquot_b), changes_molecule=False, source=METHODS)

    spri_a = g.input("spri_0_6_1_2", "0.6× left-sided then 0.6–1.2× double-sided bead size selection",
                     MaterialKind.SELECTION, source=METHODS)
    rna_branch = g.product("rna_branch", "0.6× left-sided fraction for the mRNA library",
                           MaterialKind.MOLECULAR, state=preamp.state)
    guide_branch = g.product("guide_branch", "0.6–1.2× fraction for index gRNA library I",
                             MaterialKind.MOLECULAR, state=preamp.state)
    small_waste = g.product("below_1_2x", "material not bound at 1.2×",
                            MaterialKind.WASTE)
    g.react("size_select", "Size-select aliquot A into mRNA and gRNA-library-I fractions",
            StepKind.SELECTION,
            ((aliquot_a, InputRole.TEMPLATE), (spri_a, InputRole.SELECTION)),
            (rna_branch, guide_branch, small_waste), changes_molecule=False, source=METHODS)
    g.discard(small_waste, "left in the 1.2× supernatant")
    spri_b = g.input("ampure_1_2", "1.2× AMPure bead purification", MaterialKind.SELECTION,
                     source=METHODS)
    guide_branch_2 = g.product("guide_branch_2", "purified cDNA for index gRNA library II",
                               MaterialKind.MOLECULAR, state=preamp.state)
    b_waste = g.product("aliquot_b_supernatant", "1.2× AMPure supernatant", MaterialKind.WASTE)
    g.react("purify_b", "Purify aliquot B with 1.2× AMPure beads", StepKind.CLEANUP,
            ((aliquot_b, InputRole.TEMPLATE), (spri_b, InputRole.SELECTION)),
            (guide_branch_2, b_waste), changes_molecule=False, source=METHODS)
    g.discard(b_waste, "left in the bead supernatant")
    g.mark_terminal(Terminal("guide_branch_2", Disposition.UNRESOLVED,
                             "gRNA library II receives the same 2nd and 3rd PCRs as library I; that copy is not drawn separately."))

    reviewed = ["cells", "gems", "partition", "bead_primer", "rt",
                "reverse_transcribe", "tenx_cdna_primers", "tenx_cdna_primers_2", "preamp_polymerase",
                "preamplify", "split", "spri_0_6_1_2", "size_select", "ampure_1_2", "purify_b"]
    manufacturer=g.input("rna_manufacturer_route","Chromium CG000184 mRNA-library preparation",MaterialKind.CONDITION,source=METHODS)
    rna_library=g.product("rna_library","manufacturer-prepared mRNA library",MaterialKind.LIBRARY,state=RNA_LIB)
    g.react("rna_library_prep","Prepare mRNA library following the cited Chromium guide",StepKind.MAGIC,
        ((rna_branch,InputRole.TEMPLATE),(manufacturer,InputRole.CONDITION)),(rna_library,),
        reason="Published Methods delegate the mRNA library to CG000184. Its actual fragmentation/end preparation/ligation/index PCR route is not reconstructed by a generic terminal PCR.")
    g.sequence(rna_library,SequencingPlan(library=RNA_LIB,limitation=SEQUENCING_UNAVAILABLE["Gene-expression library"]))
    reviewed.extend(("rna_manufacturer_route","rna_library_prep","rna_library"))

    p5 = g.input("p5_read1", "P5_Read1 forward primer", MaterialKind.OLIGO,
                 state=P5_READ1)
    trna = g.input("trna_read2", "tRNA_Read2 reverse primer", MaterialKind.OLIGO,
                   state=TRNA_READ2)
    polymerase_1 = g.input("guide_polymerase_1", "NEBNext Ultra II Q5 Master Mix; 8 × 50 µl reactions, 14 cycles (60 °C anneal)",
                           MaterialKind.CATALYST, source=METHODS)
    first = g.product("first_guide_amplicon", "first nested guide-PCR amplicon",
                      MaterialKind.MOLECULAR, state=FIRST_GUIDE_AMPLICON)
    first_raw = g.product("first_guide_pcr", "pooled 2nd-PCR reactions",
                          MaterialKind.MOLECULAR, state=FIRST_GUIDE_AMPLICON)
    g.react("guide_pcr_1", "Enrich guide cDNA with P5_Read1 and tRNA_Read2 (the 2nd PCR)",
            StepKind.PCR,
            ((guide_branch, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (trna, InputRole.STRUCTURAL), (polymerase_1, InputRole.CATALYST)),
            (first_raw,), geometry=PCRGeometry.unavailable(guide_branch.state,reason="This is the actual size-selected preamplified mixture. The Table S3d primers are explicit; the native first-PCR handle and its complete pairing on the incoming guide species are not reconstructed. The separately drawn tRNA binding view checks only that local site, not complete first-cycle geometry."))
    spri_c = g.input("spri_0_7_1_0", "0.7–1.0× double-sided size selection, eluted in 80 µl",
                     MaterialKind.SELECTION, source=METHODS)
    pcr1_waste = g.product("pcr2_excluded", "fragments outside the 0.7–1.0× window and primers",
                           MaterialKind.WASTE)
    g.react("guide_pcr_1_cleanup", "Size-select pooled 2nd-PCR products",
            StepKind.SELECTION, ((first_raw, InputRole.TEMPLATE), (spri_c, InputRole.SELECTION)),
            (first, pcr1_waste), changes_molecule=False, source=METHODS)
    g.discard(pcr1_waste, "outside the double-sided bead window")

    indexed = g.input("p7_indexed", "P7_Read2 indexed reverse primer",
                      MaterialKind.OLIGO, state=P7_READ2_INDEX1)
    polymerase_2 = g.input("guide_polymerase_2", "NEBNext Ultra II Q5 Master Mix; 5 × 50 µl reactions, 5 cycles (54 °C anneal, 65 °C extension)",
                           MaterialKind.CATALYST, source=METHODS)
    guide_raw = g.product("guide_pcr_products", "indexed 3rd-PCR products",
                          MaterialKind.LIBRARY, state=GUIDE_LIB)
    g.react("guide_pcr_2", "Add the indexed P7 end to the enriched guide library (the 3rd PCR)",
            StepKind.PCR,
            ((first, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (indexed, InputRole.STRUCTURAL), (polymerase_2, InputRole.CATALYST)),
            (guide_raw,), geometry=INDEX_GUIDE_PCR)
    spri_d = g.input("ampure_0_7_1_0", "0.7–1.0× AMPure XP double-sided size selection",
                     MaterialKind.SELECTION, source=METHODS)
    guide_library = g.product("guide_library", "indexed guide sequencing library",
                              MaterialKind.LIBRARY, state=GUIDE_LIB)
    pcr2_waste = g.product("pcr3_excluded", "fragments outside the 0.7–1.0× window and primers",
                           MaterialKind.WASTE)
    g.react("guide_library_cleanup", "Size-select the index gRNA library",
            StepKind.SELECTION, ((guide_raw, InputRole.TEMPLATE), (spri_d, InputRole.SELECTION)),
            (guide_library, pcr2_waste), changes_molecule=False, source=METHODS)
    g.discard(pcr2_waste, "outside the double-sided bead window")
    g.sequence(guide_library, SequencingPlan(
        library=GUIDE_LIB, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_lengths=READ_LENGTHS.get("Enriched guide library", {}),
        read_uses=READ_USES["Enriched guide library"]))
    reviewed.extend(("p5_read1", "trna_read2", "guide_polymerase_1",
                     "guide_pcr_1", "spri_0_7_1_0", "guide_pcr_1_cleanup", "p7_indexed",
                     "guide_polymerase_2", "guide_pcr_2", "ampure_0_7_1_0",
                     "guide_library_cleanup", "guide_library"))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
