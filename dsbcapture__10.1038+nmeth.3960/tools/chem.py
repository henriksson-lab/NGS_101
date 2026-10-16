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
        ("PCR-add P5, i7 and P7; size-select the library",
         molecule("DSBCapture sequencing library", (role("P5 and Read 1 site"), role("DNA beginning at the DSB"), role("Read 2 and Index 1 site"), role("six-base sample index"), role("P7"))),
         "The original detailed method prints PCR F and the indexed PCR R family and specifies a six-base, single-index TruSeq library."),
    ))


TITLE = 'DSBCapture — direct in-situ sequencing of DNA breaks'
NOTES = '01_dsbcapture.html'
SOURCE = 'Lensing et al., Nature Methods (2016); <a href="https://doi.org/10.1038/nmeth.3960">doi:10.1038/nmeth.3960</a>. The original <a href="https://doi.org/10.1038/protex.2016.52">Protocol Exchange method</a> prints the oligo set on pp. 8–9.'
SUMMARY = 'A biotinylated modified P5 adapter is ligated to repaired breaks inside fixed nuclei; after sonication and capture, modified P7 is installed only at the random distal end, so Read 1 begins at the DSB.'
CAVEAT = "<b>Source scope.</b> The paper depicts the oligos in raster figures. The original Protocol Exchange method, doi:10.1038/protex.2016.52, prints the modified P5/P7 and control duplexes, PCR primers, and six-base single-index design on PDF pp. 7–9. The model derives common TruSeq bases from shared constants; an explicit index read on the instrument is not specified by the detailed method."
FINAL_CAPTION = 'Published modified-P5/modified-P7 library; the proximal orientation is shown. Read 1 starts at the repaired native break. The six-base sample index is printed in the original detailed method.'
READ_LENGTHS = {'Read 1': 75, 'Index 1 (i7)': 6, 'Read 2': 75}
READ_USES = {'Index 1 (i7)': 'undetermined'}
SEQUENCING_INTRO = "The paper used 2×75-bp NextSeq runs (and 2×50-bp MiSeq for the EcoRV pilot); downstream break calling retains the proximal Read 1. The linked original detailed method specifies a six-base single index; whether and how the instrument acquired its separate index read remains unstated."
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _illumina_break_library("DSBCapture library", inferred=False, tagged=False, index_inferred=False)
def oligos():
    return (f'modified P5 top — 5′-/5Biosg/{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′', f'modified P5 bottom — 5′-/5Phos/{DSB_P5_BOTTOM}-3′', f'modified P7 top — 5′-/5Phos/{DSB_P7}-3′', f'modified P7 bottom — 5′-{DSB_P7_BOTTOM[:-1]}*{DSB_P7_BOTTOM[-1]}-3′', f'PCR F — 5′-{DSB_PCR_F}-3′', f'PCR R (AD06 example, one of eight printed) — 5′-{DSB_PCR_R_AD06}-3′', f'control modified P5 top — 5′-{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′ (original detailed method, PDF p. 8; no biotin; same phosphorylated complement)')

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
                step_label="Amplify with PCR F and a published six-base AD-indexed PCR R")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
