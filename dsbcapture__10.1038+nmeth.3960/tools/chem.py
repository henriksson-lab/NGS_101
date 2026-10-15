"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, da_tailed, molecule, role
from dsb_capture import BreakCaptureDesign, _illumina_break_library, _finish_pcr

DSB_P5 = il.TRUSEQ_P5_FULL                 # modified P5 top strand as drawn, 58 nt

DSB_P7 = il.INDEX1_PRIMER                  # modified P7 top strand as drawn, 33 nt

DSB_P5_BOTTOM = revcomp(DSB_P5)[1:]        # the 3′ T of the P5 top strand stays unpaired

DSB_P7_BOTTOM = revcomp(DSB_P7) + "T"      # the P7 duplex carries its T on the bottom strand

DSB_PCR_F = DSB_P5[:44]                    # drawn as the first 44 bases of the P5 adapter

DSB_AD006 = "GCCAAT"                       # TruSeq sample index AD006, read as the i7 index

DSB_PCR_R_AD06 = il.P7 + revcomp(DSB_AD006) + revcomp(DSB_P7)[:21]

_DSB_FIGURE_ONLY = ("Lensing et al. draw this oligo in Fig. 1a / Supplementary Fig. 1c as"
                    " raster artwork; no base of the DSBCapture oligo set occurs as text in"
                    " the author manuscript, Supplementary Figures or Supplementary Tables,"
                    " so there is no extracted source text to compare against.")

SOURCE_OLIGO_EXCEPTIONS = {
    "modified_p5": _DSB_FIGURE_ONLY,
    "modified_p7": _DSB_FIGURE_ONLY,
    "library_forward": _DSB_FIGURE_ONLY,
    "library_reverse": _DSB_FIGURE_ONLY + " AD006 is one of three drawn example primers.",
}


def workflow():
    return chain(molecule("fixed nuclei", (role("genomic DNA with a native DSB"),)), (
        ("Blunt native break ends and add one 3′ dA",
         molecule("A-tailed native break", (role("repaired DSB"), role("3-prime dA"))),
         "End repair normalizes the break; Klenow exo− then creates a cohesive ligation end."),
        ("Ligate the biotinylated, T-overhang modified P5 adapter in situ",
         molecule("P5-labelled break", (role("biotinylated modified P5 adapter"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated modified P5 adapter", "break-adjacent genomic DNA", "T/A ligation"),)),
         "Because the Read 1 site is already at the native break, R1 starts with the first genomic base."),
        ("Digest free adapter with lambda exonuclease, extract and sonicate",
         molecule("labelled sonication fragment", (role("modified P5 at native break"), role("break-to-random genomic fragment"), role("sonication-created free end"))),
         "Covalently break-labelled fragments retain biotin; free adapters are removed before extraction."),
        ("Streptavidin-capture, blunt and add one 3-prime dA to the distal end",
         da_tailed("distal-end dA-tailed DSBCapture fragment", "modified-P5-labelled break-to-random genomic DNA"),
         "The native-break P5 end is already labelled; only the sonication-created distal end is prepared."),
        ("Ligate modified P7 to the dA-tailed distal end",
         molecule("dual-adapter DSBCapture fragment", (role("modified P5"), role("break-to-random genomic DNA"), role("modified P7")),
                  junctions=(("modified P5", "break-to-random genomic DNA", "native-break ligation"),
                             ("break-to-random genomic DNA", "modified P7", "distal T/A ligation"))),
         "The second adapter is added only to the random sonication end."),
        ("INFERRED — PCR-add P5, i7 and P7; size-select the library",
         molecule("DSBCapture sequencing library", (role("P5 and Read 1 site"), role("DNA beginning at the DSB"), role("Read 2 and Index 1 site"), role("six-base sample index", inferred=True), role("P7"))),
         "The drawn PCR primers complete a single-index Illumina library; the index window's six-base width is measured from Supplementary Fig. 1c artwork and never stated in the text."),
    ))


TITLE = 'DSBCapture — direct in-situ sequencing of DNA breaks'
NOTES = '01_dsbcapture.html'
SOURCE = 'Lensing et al., Nature Methods (2016); <a href="https://doi.org/10.1038/nmeth.3960">doi:10.1038/nmeth.3960</a>. The oligo set below is published only as figure artwork (Fig. 1a, Supplementary Fig. 1c).'
SUMMARY = 'A biotinylated modified P5 adapter is ligated to repaired breaks inside fixed nuclei; after sonication and capture, modified P7 is installed only at the random distal end, so Read 1 begins at the DSB.'
CAVEAT = "<b>Oligo provenance.</b> Lensing et al. print the modified P5, modified P7 and control modified P5 duplexes and the PCR primers only as raster artwork in Fig. 1a and Supplementary Fig. 1c: no base of them occurs as text in the author manuscript, Supplementary Figures or Supplementary Tables, so the bases shown here cannot be re-checked against extracted source text. Every drawn base is a canonical Illumina TruSeq part, so this page derives the set from the shared <code>illumina</code> constants and from each duplex's own top strand; a figure transcription and a reconstruction out of those constants are therefore indistinguishable. The sample-index window is measured from the same artwork — the paper never states an index length, an index read or demultiplexing — and is shaded as inferred."
FINAL_CAPTION = 'Published modified-P5/modified-P7 library; the proximal orientation is shown. Read 1 starts at the repaired native break. The shaded sample-index window is the figure-measured inference described above.'
READ_LENGTHS = {'Read 1': 75, 'Index 1 (i7)': 6, 'Read 2': 75}
READ_USES = {'Index 1 (i7)': 'undetermined'}
SEQUENCING_INTRO = "The paper used 2×75-bp NextSeq runs (and 2×50-bp MiSeq for the EcoRV pilot); downstream break calling retains the proximal Read 1. The indexed PCR R implies an index read, but the paper never describes one, so the Index 1 length is taken from the figure's index window and its handling is left undetermined."
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _illumina_break_library("DSBCapture library", inferred=False, tagged=False, index_inferred=True)
def oligos():
    return (f'modified P5 top — 5′-/5Biosg/{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′', f'modified P5 bottom — 5′-/5Phos/{DSB_P5_BOTTOM}-3′', f'modified P7 top — 5′-/5Phos/{DSB_P7}-3′', f'modified P7 bottom — 5′-{DSB_P7_BOTTOM[:-1]}*{DSB_P7_BOTTOM[-1]}-3′', f'PCR F — 5′-{DSB_PCR_F}-3′', f'PCR R (AD006 example, one of three drawn) — 5′-{DSB_PCR_R_AD06}-3′', 'control modified P5 — Supplementary Fig. 1c also publishes a non-biotinylated control modified P5 adapter, the negative control of Supplementary Fig. 1e (lane C). Its bases are not transcribed here: the artwork is the only source, and whether it differs from the modified P5 above by nothing but the omitted 5′ biotin is not stated anywhere in the text.')

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    nuclei = graph.input("nuclei", "fixed nuclei containing native DNA breaks",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("fixed nuclei", (role("genomic DNA with native DSBs"),)),
                         primary=True)
    repair = graph.input("repair_enzyme", "native-break blunting enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt repaired native breaks",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("repair", "Blunt native break ends", StepKind.END_REPAIR,
                ((nuclei, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed repaired native breaks",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("tail", "Add one 3-prime dA to native break ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    p5 = graph.input("modified_p5", "biotinylated T-overhang modified P5 adapter",
                     MaterialKind.ADAPTER, state=(DSB_P5, DSB_P5_BOTTOM))
    p5_ligase = graph.input("p5_ligase", "modified-P5 adapter ligase", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "modified-P5-labelled native breaks",
                             MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_p5", "Ligate modified P5 at the native break",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p5_ligase, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "DSBCapture modified P5", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The modified P5 adapter carries a single 3′ T complementary to the dA-tailed native break."))
    exonuclease = graph.input("lambda_exonuclease", "lambda exonuclease",
                              MaterialKind.CATALYST)
    cleaned = graph.product("cleaned", "covalently P5-labelled break molecules",
                            MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("remove_free_adapter", "Digest free modified-P5 adapter",
                StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (cleaned,))
    shear = graph.input("sonication_condition", "DNA sonication condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "P5-labelled break-to-sonication fragments",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    graph.react("sonicate", "Sonicate extracted break-labelled DNA", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "captured modified-P5 break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled sonication fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated break fragments",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the biotinylated native-break end")
    distal_repair = graph.input("distal_repair", "distal-end repair enzyme mix",
                                MaterialKind.CATALYST)
    distal_blunt = graph.product("distal_blunt", "blunt sonication-created distal ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the sonication-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_blunt,))
    distal_tailing = graph.input("distal_tailing", "distal single-dA tailing enzyme",
                                 MaterialKind.CATALYST)
    distal_tailed = graph.product("distal_tailed", "dA-tailed distal break fragments",
                                  MaterialKind.MOLECULAR, state=distal_blunt.state)
    graph.react("tail_distal", "Add one 3-prime dA to the distal end", StepKind.TAILING,
                ((distal_blunt, InputRole.TEMPLATE),
                 (distal_tailing, InputRole.CATALYST)), (distal_tailed,))
    p7 = graph.input("modified_p7", "T-overhang modified P7 adapter", MaterialKind.ADAPTER,
                     state=(DSB_P7, DSB_P7_BOTTOM))
    p7_ligase = graph.input("p7_ligase", "modified-P7 adapter ligase", MaterialKind.CATALYST)
    dual = graph.product("dual", "modified-P5/P7 break fragments", MaterialKind.MOLECULAR,
                         state=distal_tailed.state)
    graph.react("ligate_p7", "Ligate modified P7 only at the distal end",
                StepKind.LIGATION,
                ((distal_tailed, InputRole.TEMPLATE), (p7, InputRole.STRUCTURAL),
                 (p7_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "DSBCapture modified P7", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The modified P7 adapter carries the T overhang complementary to the dA-tailed sonication end."))
    reviewed = ["nuclei", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "modified_p5", "p5_ligase", "ligate_p5", "lambda_exonuclease",
                "remove_free_adapter", "sonication_condition", "sonicate", "streptavidin",
                "capture", "distal_repair", "repair_distal", "distal_tailing", "tail_distal",
                "modified_p7", "p7_ligase", "ligate_p7"]
    _finish_pcr(graph, "library", dual, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label="DSBCapture",
                forward_state=DSB_PCR_F, reverse_state=DSB_PCR_R_AD06,
                step_label="INFERRED — amplify with PCR F and an AD-indexed PCR R; the"
                           " sample-index window is measured from Supplementary Fig. 1c"
                           " artwork, not stated in the text")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
