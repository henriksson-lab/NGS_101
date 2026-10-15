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
from topology import HairpinTopology

BLESS_PROXIMAL = "TACTACCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGGTAGTA"

BLESS_DISTAL = "CGTCGTCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGACGACG"

BLESS_PCR_PROXIMAL = "CCCTAGCGTAACTCTCGAGGTAGTA"

BLESS_PCR_DISTAL = "CTAGCGTAACTCTCGAGACGACG"

BLESS_PROXIMAL_READ_TAG = "TCGAGGTAGTA"

BLESS_DISTAL_READ_TAG = "TCGAGACGACG"

DESIGN = BreakCaptureDesign('BLESS', 'biotinylated proximal hairpin', 'distal hairpin', 'intact fixed nuclei', False)

def workflow():
    key = 'bless'
    d = DESIGN
    initial = molecule("yeast cells" if key == "i-bless" else d.immobilization,
                       (role("genomic DNA with a native DSB"),))
    prefix = ((
        "Embed, spheroplast and deproteinate cells in agarose beads",
        molecule("agarose-protected genomic DNA", (role("native DSB in deproteinated DNA"),)),
        "Agarose-bead immobilization is the defining i-BLESS change and protects yeast DNA from handling-induced shearing.",
    ),) if key == "i-bless" else ()
    return chain(initial, prefix + (
        ("Blunt and 5′-phosphorylate native break ends",
         molecule("ligation-ready native breaks", (role("blunt 5-prime-phosphorylated DSB"),)),
         "End repair makes 5′ overhangs, 3′ overhangs and blunt breaks enter the same ligation state."),
        ("Ligate the biotinylated proximal hairpin in situ",
         molecule("proximal-labelled break", (role("biotinylated proximal hairpin"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated proximal hairpin", "break-adjacent genomic DNA", "T4 DNA ligase"),)),
         "The proximal barcode marks the end that existed in the cell before extraction."),
        ("Fragment DNA and capture proximal-labelled molecules on streptavidin",
         molecule("captured break fragment", (role("biotinylated proximal hairpin"), role("break-to-random genomic fragment"), role("fragmentation-created free end"))),
         "Fragmentation creates the distal end only after the native break has already been labelled."),
        ("Blunt the free end and ligate the distal hairpin",
         molecule("dual-hairpin break fragment", (role("proximal hairpin"), role("break-to-random genomic fragment"), role("distal hairpin")),
                  junctions=(("proximal hairpin", "break-to-random genomic fragment", "native-break ligation"),
                             ("break-to-random genomic fragment", "distal hairpin", "distal ligation"))),
         "Different proximal and distal barcodes preserve which end is the native DSB."),
        ("Linearize with I-SceI, PCR-amplify, then cut with XhoI",
         molecule("barcoded BLESS insert", (role("proximal 11-nt tag"), role("break-to-random genomic DNA"), role("distal 11-nt tag"))),
         "I-SceI releases the captured molecule from its hairpins; XhoI leaves the short end-identifying tags."),
        ("INFERRED — install the sequencing-kit ends",
         molecule("Illumina-compatible break library", (role("Illumina left end", inferred=True), role("proximal or distal tag"), role("genomic insert"), role("Illumina right end", inferred=True))),
         "The paper names TruSeq DNA Sample Preparation Kit v2 (Methods, Next-generation sequencing) but does not print its complete adapter construct."),
    ))


TITLE = 'BLESS — Breaks Labeling, Enrichment on Streptavidin and Sequencing'
NOTES = '01_bless.html'
SOURCE = 'Crosetto et al., Nature Methods (2013); <a href="https://doi.org/10.1038/nmeth.2408">doi:10.1038/nmeth.2408</a>.'
SUMMARY = 'A barcoded biotin hairpin labels repaired DSBs in fixed nuclei before fragmentation; streptavidin capture and a differently barcoded distal hairpin preserve the native-break end.'
CAVEAT = "The proximal/distal tags are exact. The paper names TruSeq DNA Sample Preparation Kit v2 (Methods, Next-generation sequencing) but does not print its complete adapter geometry or selected index member; those outer regions are marked inferred."
FINAL_CAPTION = 'One insert orientation is shown. Proximal and distal 11-nt tags identify which genomic end was the native break; inferred outer sequencing-kit regions are shaded.'
READ_LENGTHS = {'Read 1': 75, 'Index 1 (i7)': 6, 'Read 2': 75}
SEQUENCING_INTRO = 'BLESS used paired-end Illumina sequencing; either end can carry the proximal tag, which is selected computationally.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _illumina_break_library('BLESS library', inferred=True, tagged=True, proximal_tag=BLESS_PROXIMAL_READ_TAG, distal_tag=BLESS_DISTAL_READ_TAG, index_length=6)
def oligos():
    return (f'proximal hairpin — 5′-/5Phos/{BLESS_PROXIMAL[:45]}/iBiodT/{BLESS_PROXIMAL[46:]}-3′', f'distal hairpin — 5′-/5Phos/{BLESS_DISTAL}-3′', f'proximal PCR primer — 5′-{BLESS_PCR_PROXIMAL}-3′', f'distal PCR primer — 5′-{BLESS_PCR_DISTAL}-3′')

def protocol_graph():
    key = 'bless'
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    sample = graph.input("sample", "yeast cells" if key == "i-bless" else "intact fixed nuclei",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("DSB-bearing sample", (role("genomic DNA with native DSBs"),)),
                         primary=True)
    current = sample
    reviewed = ["sample"]
    if key == "i-bless":
        deproteinization = graph.input("deproteinization_condition",
                                       "agarose immobilization, spheroplasting and deproteinization condition",
                                       MaterialKind.CONDITION)
        protected = graph.product("protected_dna", "agarose-protected deproteinated yeast DNA",
                                  MaterialKind.MOLECULAR, state=sample.state)
        graph.react("immobilize", "Immobilize and deproteinate cells in agarose beads",
                    StepKind.DEPROTEINIZATION,
                    ((sample, InputRole.TEMPLATE),
                     (deproteinization, InputRole.CONDITION)), (protected,))
        current = protected
        reviewed.extend(("deproteinization_condition", "immobilize"))
    repair = graph.input("proximal_repair", "native-break blunting and phosphorylation enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_breaks", "blunt 5-prime-phosphorylated native breaks",
                             MaterialKind.MOLECULAR, state=current.state)
    graph.react("repair_breaks", "Blunt and 5-prime-phosphorylate native breaks",
                StepKind.END_REPAIR,
                ((current, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    proximal = graph.input("proximal_adapter", "biotinylated proximal hairpin",
                           MaterialKind.ADAPTER)
    proximal_ligase = graph.input("proximal_ligase", "proximal-hairpin DNA ligase",
                                  MaterialKind.CATALYST)
    labelled = graph.product("proximal_labelled", "proximal-hairpin-labelled native breaks",
                             MaterialKind.MOLECULAR,
                             state=molecule("proximal-labelled break", (
                                 role("biotinylated proximal hairpin"),
                                 role("break-adjacent genomic DNA")),
                                 junctions=(("biotinylated proximal hairpin",
                                             "break-adjacent genomic DNA",
                                             "proximal ligation"),)))
    graph.react("ligate_proximal", "Ligate the proximal hairpin in situ",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (proximal, InputRole.STRUCTURAL),
                 (proximal_ligase, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "BLESS proximal hairpin", kind="blunt",
                    statement="The repaired native DSB is blunt and 5′ phosphorylated before proximal-hairpin ligation."))
    fragmentation = graph.input("fragmentation_condition", "genomic DNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "proximal-labelled break-to-random fragments",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment", "Fragment DNA after proximal break labelling",
                StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragments,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "captured proximal-labelled break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled genomic fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture proximal-labelled molecules on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the biotinylated native-break adapter")
    distal_repair = graph.input("distal_repair", "distal-end blunting enzyme mix",
                                MaterialKind.CATALYST)
    distal_ready = graph.product("distal_ready", "blunt captured distal fragment ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the fragmentation-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_ready,))
    distal = graph.input("distal_adapter", "distinct distal hairpin",
                         MaterialKind.ADAPTER)
    distal_ligase = graph.input("distal_ligase", "distal-hairpin DNA ligase",
                                MaterialKind.CATALYST)
    dual = graph.product("dual", "dual-hairpin break fragments", MaterialKind.MOLECULAR,
                         state=molecule("dual-hairpin break fragment", (
                             role("proximal hairpin"), role("break-to-random genomic DNA"),
                             role("distal hairpin")),
                             junctions=(("proximal hairpin", "break-to-random genomic DNA",
                                         "native-break ligation"),
                                        ("break-to-random genomic DNA", "distal hairpin",
                                         "distal ligation"))))
    graph.react("ligate_distal", "Ligate the distinct distal hairpin",
                StepKind.LIGATION,
                ((distal_ready, InputRole.TEMPLATE), (distal, InputRole.STRUCTURAL),
                 (distal_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "BLESS distal hairpin", kind="blunt",
                    statement="The fragmentation-created distal end is blunted before the distinct distal hairpin is ligated."))
    iscei = graph.input("iscei", "I-SceI hairpin-release nuclease", MaterialKind.CATALYST)
    released = graph.product("released", "linearized dual-tag break inserts",
                             MaterialKind.MOLECULAR, state=dual.state)
    graph.react("linearize", "Linearize the captured hairpin product with I-SceI",
                StepKind.DIGESTION,
                ((dual, InputRole.TEMPLATE), (iscei, InputRole.CATALYST)), (released,))
    xhoi = graph.input("xhoi", "XhoI tag-trimming nuclease", MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "break inserts retaining proximal/distal tags",
                            MaterialKind.MOLECULAR, state=released.state)
    graph.react("trim_tags", "Cut with XhoI while retaining end-identifying tags",
                StepKind.DIGESTION,
                ((released, InputRole.TEMPLATE), (xhoi, InputRole.CATALYST)), (trimmed,))
    reviewed.extend(("proximal_repair", "repair_breaks", "proximal_adapter",
                     "proximal_ligase", "ligate_proximal", "fragmentation_condition",
                     "fragment", "streptavidin", "capture", "distal_repair",
                     "repair_distal", "distal_adapter", "distal_ligase",
                     "ligate_distal", "iscei", "linearize", "xhoi", "trim_tags"))
    _finish_pcr(graph, "library", trimmed, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label=TITLE)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


STEM_NT = 42
PROXIMAL_TOPOLOGY = HairpinTopology(
    "BLESS L1 proximal hairpin", BLESS_PROXIMAL, 0,
    len(BLESS_PROXIMAL) - STEM_NT, STEM_NT,
    colour="var(--c-r1)")
DISTAL_TOPOLOGY = HairpinTopology(
    "BLESS L3 distal hairpin", BLESS_DISTAL, 0,
    len(BLESS_DISTAL) - STEM_NT, STEM_NT,
    colour="var(--c-r2)")
