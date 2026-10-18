"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from dsb_capture import BreakCaptureDesign
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
        ("Fragment DNA with HaeIII and capture proximal-labelled molecules on streptavidin",
         molecule("captured break fragment", (role("biotinylated proximal hairpin"), role("break-to-random genomic fragment"), role("fragmentation-created free end"))),
         "HaeIII fragmentation creates the distal end only after the native break has already been labelled."),
        ("Ligate the distal hairpin to the HaeIII end of captured fragments",
         molecule("dual-hairpin break fragment", (role("proximal hairpin"), role("break-to-random genomic fragment"), role("distal hairpin")),
                  junctions=(("proximal hairpin", "break-to-random genomic fragment", "native-break ligation"),
                             ("break-to-random genomic fragment", "distal hairpin", "distal ligation"))),
         "Different proximal and distal barcodes preserve which end is the native DSB."),
        ("Release with I-SceI, PCR-amplify with P1/P3, then cut with XhoI",
         molecule("barcoded BLESS insert", (role("proximal 11-nt tag"), role("break-to-random genomic DNA"), role("distal 11-nt tag"))),
         "I-SceI releases the captured molecule from its hairpins; XhoI leaves the short end-identifying tags."),
        ("INFERRED — install the sequencing-kit ends",
         molecule("Illumina-compatible break library", (role("Illumina left end", inferred=True), role("proximal or distal tag"), role("genomic insert"), role("Illumina right end", inferred=True))),
         "The paper names TruSeq DNA Sample Preparation Kit v2, used without fragmentation or size selection (Methods, Next-generation sequencing), but does not print its complete adapter construct."),
    ))


TITLE = 'BLESS — Breaks Labeling, Enrichment on Streptavidin and Sequencing'
NOTES = '01_bless.html'
SOURCE = 'Crosetto et al., Nature Methods (2013); <a href="https://doi.org/10.1038/nmeth.2408">doi:10.1038/nmeth.2408</a>.'
SUMMARY = 'A barcoded biotin hairpin labels repaired DSBs in fixed nuclei before fragmentation; streptavidin capture and a differently barcoded distal hairpin preserve the native-break end.'
CAVEAT = "The proximal/distal tags are exact. The paper names TruSeq DNA Sample Preparation Kit v2 (Methods, Next-generation sequencing); its exact supplied adapter geometry and selected index member/length remain unassigned in this model, so those outer regions are inferred."
FINAL_CAPTION = 'One insert orientation is shown. Proximal and distal 11-nt tags identify which genomic end was the native break; inferred outer sequencing-kit regions are shaded.'
# Supplementary Table 1 (published Springer supplement) lists GAII single-end runs for most
# BLESS samples and HiSeq 2000 paired-end runs for N1/N2; no read length or index read is
# stated, so no cycle counts are declared here.
SEQUENCING_INTRO = ('Supplementary Table 1 lists Illumina GAII single-end runs for most BLESS samples and HiSeq 2000 '
                    'paired-end runs (TruSeq PE Cluster Kit v3, TruSeq SBS Kit v3-HS) for the neocarzinostatin samples; '
                    'read lengths and any index read are not stated, so no cycle counts are drawn. Either end can carry the '
                    'proximal tag; reads are scanned for the exact proximal or distal barcode (Methods, Illumina data analysis).')
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
# INFERRED outer ends: the BLESS sources name the TruSeq DNA Sample Preparation Kit v2
# but print no adapter, index primer or junction; canonical TruSeq single-index bases
# (Illumina Adapter Sequences) and the kit's dA/T ligation junction are drawn as inferred.
def bless_library():
    # TruSeq v2 is source-named, but this pass does not establish its supplied
    # index length/member or end-processing details. Preserve that distinction.
    lib=Construct([Segment("P5 + Read1 arm",il.TRUSEQ_P5_FULL,"p5",inferred=True),
        Segment("proximal read tag",BLESS_PROXIMAL_READ_TAG,"cbc",feature=feature("proximal_break_end","feature_barcode","fixed")),
        role("break-to-random genomic DNA"),
        Segment("distal read tag",BLESS_DISTAL_READ_TAG,"cbc",feature=feature("distal_fragment_end","feature_barcode","fixed")),
        Segment("ligation junction","A",inferred=True),
        Segment("Index1 + Read2 arm",il.INDEX1_PRIMER,"r2",inferred=True),
        role("unassigned kit i7 index",tag="cbc",inferred=True,
             molecular_feature=feature("sample_index_i7","sample_index","unknown")),
        Segment("P7 reverse complement",il.P7_RC,"p7",inferred=True)],name="BLESS conditional TruSeq v2 library")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems:raise ValueError("BLESS library: "+"; ".join(problems))
    return lib,primers

FINAL_LIBRARY,SEQ_PRIMERS=bless_library()
# Supplementary Table 5 marks the biotinylated base by cyan highlighting only (46th base,
# in the T7 loop); it prints no 5'-phosphate on L1/L3, and the Methods do not say whether
# the annealed linkers were phosphorylated, so no /5Phos/ is displayed.
BIOTIN_INDEX = 45


def oligos():
    return (f'L1 proximal hairpin — 5′-{BLESS_PROXIMAL[:BIOTIN_INDEX]}/iBiodT/{BLESS_PROXIMAL[BIOTIN_INDEX + 1:]}-3′ (5′-phosphate not stated)',
            f'L3 distal hairpin — 5′-{BLESS_DISTAL}-3′ (5′-phosphate not stated)',
            f'P1 proximal PCR primer — 5′-{BLESS_PCR_PROXIMAL}-3′',
            f'P3 distal PCR primer — 5′-{BLESS_PCR_DISTAL}-3′')


def _primer(name, bases):
    return (Segment(name, bases),)

def hairpin_reagent(sequence, *, proximal):
    # One covalent oligo; folding is declared by the paired HairpinTopology below
    # and shown in the public molecular-detail view, not a fictitious second oligo.
    sc=Scene()
    sc.strand("hairpin oligo",[Segment("first42nt stem",sequence[:STEM_NT]),
        Segment("hairpin loop",sequence[STEM_NT:-STEM_NT],note="internal biotin-dT in L1 only" if proximal else ""),
        Segment("opposite42nt stem",sequence[-STEM_NT:])],label="L1 proximal" if proximal else "L3 distal")
    return sc


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
                           MaterialKind.ADAPTER,state=hairpin_reagent(BLESS_PROXIMAL,proximal=True),
                           chemical_state=PROXIMAL_TOPOLOGY)
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
    fragmentation = graph.input("fragmentation_enzyme", "HaeIII (6 U per million fixed cells, 18-20 h, 16 °C)",
                                MaterialKind.CATALYST)
    fragments = graph.product("fragments", "proximal-labelled break-to-HaeIII-site fragments",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment", "Extract gDNA and fragment it with HaeIII after proximal break labelling",
                StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (fragmentation, InputRole.CATALYST)),
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
    distal = graph.input("distal_adapter", "distinct distal hairpin",
                         MaterialKind.ADAPTER,state=hairpin_reagent(BLESS_DISTAL,proximal=False),
                         chemical_state=DISTAL_TOPOLOGY)
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
    graph.react("ligate_distal", "Ligate the distinct distal hairpin on the beads",
                StepKind.LIGATION,
                ((captured, InputRole.TEMPLATE), (distal, InputRole.STRUCTURAL),
                 (distal_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "BLESS distal hairpin", kind="blunt",
                    statement="The distal hairpin is ligated directly to the HaeIII-created end of captured fragments; the Methods describe no separate distal blunting step, and HaeIII (GG^CC) cuts blunt."))
    iscei = graph.input("iscei", "I-SceI hairpin-release nuclease", MaterialKind.CATALYST)
    released = graph.product("released", "linearized dual-tag break inserts",
                             MaterialKind.MOLECULAR, state=dual.state)
    graph.react("linearize", "Linearize the captured hairpin product with I-SceI",
                StepKind.DIGESTION,
                ((dual, InputRole.TEMPLATE), (iscei, InputRole.CATALYST)), (released,))
    p1 = graph.input("bless_p1", "P1 proximal PCR primer (Supplementary Table 5)",
                     MaterialKind.OLIGO, state=_primer("P1", BLESS_PCR_PROXIMAL))
    p3 = graph.input("bless_p3", "P3 distal PCR primer (Supplementary Table 5)",
                     MaterialKind.OLIGO, state=_primer("P3", BLESS_PCR_DISTAL))
    phusion = graph.input("phusion", "Phusion polymerase (18 cycles, Ta 55 °C)",
                          MaterialKind.CATALYST)
    amplified = graph.product("amplified", "P1/P3-amplified BLESS fragments (gel-purified)",
                              MaterialKind.MOLECULAR, state=released.state)
    graph.react("bless_pcr", "Amplify released fragments with linker-specific primers P1/P3",
                StepKind.PCR,
                ((released, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
                 (p3, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    released.state, reason="The released template is drawn as roles; the exact I-SceI-cut linker remnant on each end is not modelled base by base."))
    xhoi = graph.input("xhoi", "XhoI tag-trimming nuclease", MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "break inserts retaining proximal/distal tags",
                            MaterialKind.MOLECULAR, state=released.state)
    graph.react("trim_tags", "Cut with XhoI while retaining end-identifying tags, then gel-purify",
                StepKind.DIGESTION,
                ((amplified, InputRole.TEMPLATE), (xhoi, InputRole.CATALYST)), (trimmed,))
    kit = graph.input("truseq_kit", "TruSeq DNA Sample Preparation Kit v2 reagents",
                      MaterialKind.CONDITION)
    library = graph.product("library_library", TITLE + " sequencing library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True)
    graph.react("truseq_prep", "INFERRED — TruSeq DNA Sample Preparation Kit v2 library preparation without fragmentation or size selection",
                StepKind.MAGIC, ((trimmed, InputRole.TEMPLATE), (kit, InputRole.CONDITION)),
                (library,),
                reason="The paper names TruSeq DNA Sample Preparation Kit v2 but the kit original (end processing, adapter ligation, PCR primers, index member) remains uninspected; the outer adapter shell is drawn inferred.")
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS),
                                           acquired_roles=("Read 1", "Read 2"),
        limitation="Representative paired-end HiSeq2000 N1/N2 branch from Supplementary Table1; most other BLESS libraries use GAII single-end sequencing. No index acquisition/length or read cycles are inferred from this kit shell."))
    reviewed.extend(("proximal_repair", "repair_breaks", "proximal_adapter",
                     "proximal_ligase", "ligate_proximal", "fragmentation_enzyme",
                     "fragment", "streptavidin", "capture", "distal_adapter", "distal_ligase",
                     "ligate_distal", "iscei", "linearize", "bless_p1", "bless_p3", "phusion",
                     "bless_pcr", "xhoi", "trim_tags", "truseq_kit", "truseq_prep",
                     "library_library"))
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
