"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role

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
         distal_tailed_scene(),
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
FINAL_CAPTION = 'AD06 example. Published modified-P5/modified-P7 library; the proximal orientation is shown. Read 1 starts at the repaired native break. The six-base sample index is printed in the original detailed method.'
READ_LENGTHS = {'Read 1': 75, 'Read 2': 75}
READ_USES = {'Index 1 (i7)': 'undetermined'}
SEQUENCING_INTRO = "The paper used 2×75-bp NextSeq runs (and 2×50-bp MiSeq for the EcoRV pilot); downstream break calling retains the proximal Read 1. The linked original detailed method specifies a six-base single index; whether and how the instrument acquired its separate index read remains unstated."
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
# Ends from the printed oligo list (Protocol Exchange, pp. 8-9): the modified P5 top strand,
# the dA-tailed distal end ligated to the modified P7 top strand (its Index 1/Read 2 arm),
# and the P7 end copied from the 5' 24 bases shared by every printed PCR R primer.
DSB_P7_END = revcomp(DSB_PCR_R_AD06[:len(il.P7)])
def dsbcapture_library():
    # Original Protocol Exchange pp8–9 explicitly supplies this six-base AD06 example.
    lib=Construct([Segment("P5 + Read1 arm",DSB_P5,"p5"),
        role("DNA beginning at repaired DSB"),Segment("ligation junction","A"),
        Segment("Index1 + Read2 arm",DSB_P7,"r2"),
        Segment("6-base i7 reverse complement",DSB_AD006,"cbc",
                feature=feature("sample_index_i7","sample_index","whitelist",whitelist="published DSBCapture AD primer set")),
        Segment("P7 reverse complement",DSB_P7_END,"p7")],name="DSBCapture AD06 library example")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems:raise ValueError("DSBCapture library: "+"; ".join(problems))
    return lib,primers


def pre_index_template():
    return Construct([Segment("modified P5",DSB_P5,"p5"),
        role("DNA beginning at repaired DSB"),Segment("distal dA junction","A"),
        Segment("modified P7",DSB_P7,"r2")],name="modified P5/P7 adapter-ligated PCR input")


def distal_tailed_scene():
    # Only the sonication-created distal end is prepared here. The proximal
    # adapter is already covalently installed; it must not gain a second dA end.
    body=[Segment("modified P5",DSB_P5,"p5"),Segment("genomic insert","X"*30,placeholder=True)]
    bottom=complement_segments(body)
    sc=Scene();sc.strand("top",[*body,Segment("distal3′dA","A")])
    sc.anneal("bottom",bottom,to="top",pair=(bottom[-1].name,body[0].name))
    sc.overhang("top","distal3′dA",polarity="3-prime",owner="substrate",text="distal3′dA only")
    return sc

def p5_adapter_scene():
    sc=Scene();sc.strand("P5top",[Segment("P5 paired core",DSB_P5[:-1]),Segment("P5terminalT",DSB_P5[-1])],mod5="biotin")
    sc.anneal("P5bottom",[Segment("P5 complement",DSB_P5_BOTTOM)],to="P5top",pair=("P5 complement","P5 paired core"),mod5="p")
    sc.ligation_end("P5top","P5terminalT",kind="cohesive",polarity="3-prime",owner="adapter",text="3′T; terminal phosphorothioate")
    return sc


def p7_adapter_scene():
    sc=Scene();sc.strand("P7top",[Segment("modified P7",DSB_P7)],mod5="p")
    sc.anneal("P7bottom",[Segment("P7 complement",DSB_P7_BOTTOM[:-1]),Segment("P7terminalT",DSB_P7_BOTTOM[-1])],
        to="P7top",pair=("P7 complement","modified P7"))
    sc.ligation_end("P7bottom","P7terminalT",kind="cohesive",polarity="3-prime",owner="adapter",text="3′T; terminal phosphorothioate")
    return sc


def native_tailed_scene():
    body=[Segment("native-break-adjacent genomic DNA","X"*30,placeholder=True)]
    bottom=complement_segments(body)
    sc=Scene();sc.strand("top",[*body,Segment("native3′dA","A")])
    sc.anneal("bottom",bottom,to="top",pair=(bottom[-1].name,body[0].name))
    sc.overhang("top","native3′dA",polarity="3-prime",owner="substrate",text="repaired native-break3′dA")
    return sc


FINAL_LIBRARY,SEQ_PRIMERS=dsbcapture_library()
def oligos():
    return (f'modified P5 top — 5′-/5Biosg/{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′', f'modified P5 bottom — 5′-/5Phos/{DSB_P5_BOTTOM}-3′', f'modified P7 top — 5′-/5Phos/{DSB_P7}-3′', f'modified P7 bottom — 5′-{DSB_P7_BOTTOM[:-1]}*{DSB_P7_BOTTOM[-1]}-3′', f'PCR F — 5′-{DSB_PCR_F}-3′', f'PCR R (AD06 example, one of eight printed) — 5′-{DSB_PCR_R_AD06}-3′', f'control modified P5 top — 5′-{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′ (original detailed method, PDF p. 8; no biotin; same phosphorylated complement)')

PROTEX = ("Original Protocol Exchange method doi:10.1038/protex.2016.52 "
          "(supplement_0.pdf), steps 42–131 and oligonucleotide list")
METHODS = "Lensing et al. 2016 Online Methods, DSBCapture and Annealing of oligonucleotides"


def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    nuclei = graph.input("nuclei", "2 % formaldehyde-fixed, lysed, SDS-permeabilized and briefly proteinase K-treated nuclei containing native DNA breaks",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("fixed nuclei", (role("genomic DNA with native DSBs"),)),
                         primary=True, source=METHODS)
    repair = graph.input("repair_enzyme", "NEB Quick Blunting Kit Blunting Enzyme Mix + dNTPs, 25 °C 45 min in situ",
                         MaterialKind.CATALYST, source=PROTEX)
    repaired = graph.product("repaired", "blunt repaired native breaks",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("repair", "Blunt native break ends in situ", StepKind.END_REPAIR,
                ((nuclei, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,),
                source=PROTEX)
    tailing = graph.input("tailing_enzyme", "Klenow Fragment 3′–5′ exo− + dATP, 37 °C 45 min in situ",
                          MaterialKind.CATALYST, source=PROTEX)
    tailed = graph.product("tailed", "dA-tailed repaired native breaks",
                           MaterialKind.MOLECULAR, state=native_tailed_scene())
    graph.react("tail", "Add one 3-prime dA to native break ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,),
                source=PROTEX)
    p5 = graph.input("modified_p5", "biotinylated T-overhang modified P5 adapter (annealed 10 µM in T4 ligase buffer, 95 °C 10 min, slow cool)",
                     MaterialKind.ADAPTER, state=p5_adapter_scene(), source=METHODS)
    p5_ligase = graph.input("p5_ligase", "T4 DNA ligase, 16 °C 15–20 h in situ (a no-ligase control is run in parallel)",
                            MaterialKind.CATALYST, source=PROTEX)
    labelled = graph.product("labelled", "modified-P5-labelled native breaks",
                             MaterialKind.MOLECULAR, state=Construct([Segment("modified P5",DSB_P5,"p5"),role("native-break-adjacent genomic DNA")],name="P5-labelled native break"))
    graph.react("ligate_p5", "Ligate modified P5 at the native break",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p5_ligase, InputRole.CATALYST)), (labelled,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",p5_adapter_scene(),native_tailed_scene(),
                    statement="The modified P5 adapter carries a single 3′ T complementary to the dA-tailed native break."),
                source=PROTEX)
    exonuclease = graph.input("lambda_exonuclease", "lambda exonuclease, 50 U, 37 °C 30 min after W&B washes",
                              MaterialKind.CATALYST, source=PROTEX)
    cleaned = graph.product("cleaned", "covalently P5-labelled break molecules",
                            MaterialKind.MOLECULAR, state=labelled.state)
    adapter_digest = graph.product("free_adapter_digest", "excess unligated modified-P5 adapter digested by lambda exonuclease",
                                   MaterialKind.WASTE)
    graph.react("remove_free_adapter", "Digest excess free modified-P5 adapter",
                StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (cleaned, adapter_digest), source=PROTEX)
    graph.discard(adapter_digest, "nuclei are pelleted and washed in NEBuffer 2 + Triton X-100 after the digestion (steps 73–75)")
    extraction = graph.input("extraction", "proteinase K 200 µg/ml, 55 °C 30 min then 65 °C 30 min; potassium acetate/isopropanol precipitation, 70 % ethanol washes",
                             MaterialKind.CONDITION, source=PROTEX)
    extracted = graph.product("extracted", "purified P5-labelled genomic DNA",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    extraction_waste = graph.product("extraction_supernatant", "precipitation supernatant and ethanol washes",
                                     MaterialKind.WASTE)
    graph.react("extract", "Lyse nuclei, reverse cross-links and precipitate genomic DNA",
                StepKind.EXTRACTION,
                ((cleaned, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (extracted, extraction_waste), source=PROTEX)
    graph.discard(extraction_waste, "supernatants removed after isopropanol precipitation and ethanol washes (steps 89–98)")
    shear = graph.input("sonication_condition", "Bioruptor HIGH, 30 s on/30 s off, ~25 cycles at 4 °C; target 200–500 bp",
                        MaterialKind.CONDITION, source=PROTEX)
    fragments = graph.product("fragments", "P5-labelled break-to-sonication fragments",
                              MaterialKind.MOLECULAR, state=extracted.state)
    graph.react("sonicate", "Sonicate extracted break-labelled DNA", StepKind.FRAGMENTATION,
                ((extracted, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,),
                source=PROTEX)
    streptavidin = graph.input("streptavidin", "Dynabeads MyOne Streptavidin C1, 5 µl per 50 µg DNA, 4 °C 45 min",
                               MaterialKind.SELECTION, source=PROTEX)
    captured = graph.product("captured", "captured modified-P5 break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled sonication fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated break fragments",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected), source=PROTEX)
    graph.discard(rejected, "supernatant removed on the magnet and beads washed 3× in W&B buffer; it does not retain the biotinylated native-break end")
    distal_repair = graph.input("distal_repair", "Quick Blunting Enzyme Mix + dNTPs on beads, 25 °C 45 min",
                                MaterialKind.CATALYST, source=PROTEX)
    distal_blunt = graph.product("distal_blunt", "blunt sonication-created distal ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the sonication-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_blunt,), source=PROTEX)
    distal_tailing = graph.input("distal_tailing", "Klenow Fragment 3′–5′ exo− + dATP on beads, 37 °C 45 min",
                                 MaterialKind.CATALYST, source=PROTEX)
    distal_tailed = graph.product("distal_tailed", "dA-tailed distal break fragments",
                                  MaterialKind.MOLECULAR, state=distal_tailed_scene())
    graph.react("tail_distal", "Add one 3-prime dA to the distal end", StepKind.TAILING,
                ((distal_blunt, InputRole.TEMPLATE),
                 (distal_tailing, InputRole.CATALYST)), (distal_tailed,), source=PROTEX)
    p7 = graph.input("modified_p7", "T-overhang modified P7 adapter (annealed like modified P5)", MaterialKind.ADAPTER,
                     state=p7_adapter_scene(), source=METHODS)
    p7_ligase = graph.input("p7_ligase", "T4 DNA ligase on beads, 16 °C 15–20 h", MaterialKind.CATALYST,
                            source=PROTEX)
    dual = graph.product("dual", "modified-P5/P7 break fragments", MaterialKind.MOLECULAR,
                         state=pre_index_template())
    graph.react("ligate_p7", "Ligate modified P7 only at the distal end",
                StepKind.LIGATION,
                ((distal_tailed, InputRole.TEMPLATE), (p7, InputRole.STRUCTURAL),
                 (p7_ligase, InputRole.CATALYST)), (dual,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",p7_adapter_scene(),distal_tailed_scene(),
                    statement="The modified P7 adapter carries the T overhang complementary to the dA-tailed sonication end."),
                source=PROTEX)
    forward = graph.input("library_forward", "DSBCapture PCR F primer", MaterialKind.OLIGO,
                          state=DSB_PCR_F, source=PROTEX)
    reverse = graph.input("library_reverse", "DSBCapture PCR R primer (AD06 shown; one of eight AD-indexed primers)",
                          MaterialKind.OLIGO, state=DSB_PCR_R_AD06, source=PROTEX)
    polymerase = graph.input("library_polymerase", "NEBNext High-Fidelity 2× PCR Master Mix, 15 cycles, 65 °C annealing, on beads",
                             MaterialKind.CATALYST, source=PROTEX)
    amplified = graph.product("amplified", "amplified DSBCapture library in PCR supernatant",
                              MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    spent_beads = graph.product("spent_beads", "streptavidin beads after PCR", MaterialKind.WASTE)
    graph.react("library_pcr", "Amplify on beads with PCR F and a published six-base AD-indexed PCR R",
                StepKind.PCR,
                ((dual, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified, spent_beads),
                geometry=PCRGeometry.exact(dual.state,PCRPrimer("PCR F",DSB_PCR_F),
                    PCRPrimer("PCR R AD06",DSB_PCR_R_AD06),
                    statement="Original Protocol Exchange oligo list: PCR F44nt and AD06 PCR R bind the actual modifiedP5/P7 precursor before PCR adds the index/P7 tail."),
                source=PROTEX)
    graph.discard(spent_beads, "beads pelleted; the PCR supernatant is recovered and the beads discarded (step 127)")
    sizer = graph.input("size_selection", "MinElute cleanup then BluePippin 1.5 % cassette, 250–1200 bp (or 2 % gel excision + MinElute gel extraction)",
                        MaterialKind.SELECTION, source=PROTEX)
    final = graph.product("library_library", "DSBCapture sequencing library", MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    off_size = graph.product("off_size", "PCR products outside 250–1200 bp", MaterialKind.WASTE)
    graph.react("size_select", "Purify and size-select 250–1200 bp", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (sizer, InputRole.SELECTION)), (final, off_size),
                source=PROTEX)
    graph.discard(off_size, "excluded by BluePippin or gel band excision")
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Read 2"),read_lengths={"Read 1":75,"Read 2":75},
        read_uses={"Read 1":"include","Read 2":"include"},
        limitation="Representative2×75 NextSeq500 branch. The2×50 MiSeq EcoRV pilot is a separate run variant; six-base index content does not establish acquired index cycles."))
    reviewed = ["nuclei", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "modified_p5", "p5_ligase", "ligate_p5", "lambda_exonuclease",
                "remove_free_adapter", "extraction", "extract", "sonication_condition",
                "sonicate", "streptavidin", "capture", "distal_repair", "repair_distal",
                "distal_tailing", "tail_distal", "modified_p7", "p7_ligase", "ligate_p7",
                "library_forward", "library_reverse", "library_polymerase", "library_pcr",
                "size_selection", "size_select", "library_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
