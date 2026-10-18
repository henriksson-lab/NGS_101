"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, Workflow, feature
from circular import CircularizationGeometry
import illumina as il
from pcr import PCRGeometry
from protocol_graph import DetailView, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role
from rna_end_readout import StopReadout

PSEUDO_ADAPTER = Segment("3-prime adapter; 5-prime chemistry unresolved", "TGGAATTCTCGGGTGCCAAGG", "r2")
PSEUDO_BLOCK = Segment("adapter terminal dideoxyC", "C", note="/3ddC/: no 3′ OH")

PSEUDO_FWD = Segment("PCR forward", il.P5[:20], "p5")

PSEUDO_REV = Segment("indexed PCR reverse", il.P7 + "XXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA", "p7", placeholder=True)

# The printed reverse primer split into its tail, barcode window and RT-primer-matching end.
PSEUDO_REV_STATE = Construct((
    Segment("P7 tail", PSEUDO_REV.top[:24], "p7"),
    Segment("XXXXXX sample barcode", PSEUDO_REV.top[24:30], "cbc", placeholder=True,
            feature=feature("sample_index", "sample_index", "unknown")),
    Segment("RT-primer-matching 3-prime end", PSEUDO_REV.top[30:], "r2"),
), name="barcoded PCR reverse primer")


TITLE = 'Pseudo-seq — transcriptome-wide pseudouridine mapping'
NOTES = '01_pseudo_seq.html'
SOURCE = 'Carlile et al., Nature (2014); <a href="https://doi.org/10.1038/nature13802">doi:10.1038/nature13802</a>.'
SUMMARY = 'Persistent Ψ–CMC adducts stop reverse transcriptase; enriched cDNA termini relative to a mock-treated control locate pseudouridines.'
CAVEAT = ("<b>Published adapter discrepancy.</b> Carlile et al., Methods, ‘Pseudo-seq Library Preparation’ (PMC4224642, paragraph on 3′ ligation) call the adapter pre-adenylated, but print its IDT sequence as <code>/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/</code>. A 5′ phosphate and a 5′ adenylate are different states. The same paragraph runs the ligation with T4 RNA ligase (NEB M0204) in 1X buffer without ATP, which is consistent with an adenylated donor but does not settle the ordered modification. This page follows the printed base sequence and modification string while depicting ligation as described; whether the ordered or reacted adapter was adenylated remains unresolved.")
FINAL_LIBRARY = Construct((role('P5'), role('stop-defined RNA cDNA'), role('sample index', tag='cbc', length_bp=6, molecular_feature=feature('sample_index', 'sample_index', 'unknown')), role('P7')), name='Pseudo-seq library')
FINAL_CAPTION = 'PCR product from circularized, stop-defined cDNA.'
SEQ_PRIMERS = ()
READ_LENGTHS = None
SEQUENCING_INTRO = ''
SEQUENCING_ENDING = 'The defining experiment sequenced the amplified library on an Illumina HiSeq 2000.'
SEQUENCING_UNAVAILABLE = 'Library oligos are source printed; this model has not reconstructed the spacer-containing circular template or sequencing-primer landing. Read acquisition remains unassigned.'

def oligos():
    return ('3′ adapter — 5′-/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/-3′', 'RT primer — 5′-/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA-3′', 'PCR forward — 5′-AATGATACGGCGACCACCGA-3′', 'indexed PCR reverse — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′')

def workflow() -> Workflow:
    readout = StopReadout("CMC", "pseudouridine", "parallel mock −CMC")
    return chain(molecule("fragmented RNA", (role("RNA containing U and Ψ"),)), (  # ZnCl2 (yeast) or Zn acetate (HeLa) fragmentation precedes CMC
        ("Derivatize RNA with CMC and reverse U/G adducts",
         molecule("CMC-treated RNA", (role("RNA"), role("persistent Ψ–CMC adduct"))),
         f"{readout.control} distinguishes Ψ-dependent stops from background."),
        ("Repair ends, size-select, and ligate the reported 3′ adapter",
         molecule("adapter-ligated RNA", (role("CMC-treated RNA fragment"), PSEUDO_ADAPTER, PSEUDO_BLOCK),
                  junctions=(("CMC-treated RNA fragment", "3-prime adapter; 5-prime chemistry unresolved", "RNA ligation"),)),
         "CIP/PNK treatment precedes 8% urea-PAGE selection of 80–100, 100–120 and 120–140-nt RNA fragments; T4 RNA ligase 1 joins the adapter in ATP-free buffer."),
        ("Reverse-transcribe until AMV RT stops one nucleotide 3′ of Ψ–CMC",
         molecule("truncated cDNA", (*rt_primer_parts(), role("cDNA ending at Ψ-dependent stop"))),
         "Truncated cDNAs are gel size-selected; the 5′ end of each read reports the position of the block."),
        ("Circularize cDNA and amplify with indexed primers",
         molecule("Pseudo-seq library", tuple(FINAL_LIBRARY)),
         "CircLigase closes the single-stranded cDNA before Illumina PCR."),
    ))



PSEUDO_SOURCE = "Carlile2014 PMC4224642 Methods Pseudo-seq Library Preparation, lines54–64"
RT_HEAD = "GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT"
RT_SITE = "GCCTTGGCACCCGAGAATTCCA"


def rt_primer_parts():
    return [Segment("RT 5′ DNA component", RT_HEAD), role("first /iSp18/ non-nucleotide spacer"),
            Segment("internal DNA", "CACTCA"), role("second /iSp18/ non-nucleotide spacer"),
            Segment("RT 3′ DNA binding component", RT_SITE)]

def rt_primer_state():
    # /iSp18/ are non-nucleotide spacers, not zero-length DNA junctions.
    sc = Scene()
    sc.strand("primer", rt_primer_parts(), label="5′-phosphorylated RT primer", mod5="phosphate")
    return sc


def rt_priming_scene():
    sc = Scene()
    sc.strand("RNA", [role("CMC-treated RNA fragment"),
              Segment("full printed adapter including ddC base", PSEUDO_ADAPTER.top + PSEUDO_BLOCK.top)],
              label="reported ligated RNA", mod3="ddC (blocked)")
    sc.anneal("primer", rt_primer_parts(), to="RNA",
              pair=("RT 3′ DNA binding component", "full printed adapter including ddC base"),
              label="spacer-containing RT primer", mod5="phosphate")
    sc.arrow("primer", "AMV-RT extends toward the RNA stop")
    return sc


def _select(graph, identifier, template, label, excluded_label):
    criterion = graph.input(identifier + "_gel", label, MaterialKind.SELECTION, source=PSEUDO_SOURCE)
    selected = graph.product(identifier + "_product", label + " fraction", MaterialKind.MOLECULAR, state=template.state)
    excluded = graph.product(identifier + "_excluded", excluded_label, MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.SELECTION,
                ((template, InputRole.TEMPLATE), (criterion, InputRole.SELECTION)), (selected, excluded), changes_molecule=False)
    graph.discard(excluded, "physically excluded during the source gel-band recovery")
    return selected


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("input", "yeast total/poly(A)+ or HeLa poly(A)+ RNA (source alternatives)",
                      MaterialKind.BIOLOGICAL, primary=True,
                      state=molecule("source RNA", (role("RNA with U and Ψ"),)))
    zinc = graph.input("fragment_condition", "Yeast: 10 mM ZnCl₂, 94 °C, 5 min total/55 s poly(A)+; "
                       "HeLa: 10 mM zinc acetate, 60 °C 10 min (choose the matching source sample)", MaterialKind.CONDITION)
    fragments = graph.product("fragment_product", "zinc-fragmented RNA", MaterialKind.MOLECULAR,
                              state=molecule("RNA fragments", (role("fragmented RNA with U and Ψ"),)))
    graph.react("fragment", "Fragment RNA using the sample-matched zinc condition", StepKind.FRAGMENTATION,
                ((rna, InputRole.TEMPLATE), (zinc, InputRole.CONDITION)), (fragments,))
    cmc = graph.input("cmc_condition", "0.2/0.4 M CMC in BEU, 40 °C 30 min; parallel BEU mock without CMC",
                      MaterialKind.CONDITION)
    modified = graph.product("cmc_product", "CMC-modified RNA (with parallel mock control)", MaterialKind.MOLECULAR,
                             state=molecule("CMC-treated RNA", (role("RNA with CMC adducts on U/G/Ψ"),)))
    graph.react("cmc", "Derivatize RNA with CMC, retaining a matching mock comparison", StepKind.CHEMICAL_CONVERSION,
                ((fragments, InputRole.TEMPLATE), (cmc, InputRole.CONDITION)), (modified,))
    alkali = graph.input("reverse_cmc_condition", "50 mM Na₂CO₃ pH10.4, 2 mM EDTA, 50 °C 2 h; precipitation",
                         MaterialKind.CONDITION)
    persistent = graph.product("reverse_cmc_product", "RNA retaining persistent Ψ–CMC", MaterialKind.MOLECULAR,
                               state=molecule("reversed CMC RNA", (role("persistent Ψ–CMC with U/G adducts reversed"),)))
    graph.react("reverse_cmc", "Reverse U/G adducts after CMC treatment", StepKind.CHEMICAL_CONVERSION,
                ((modified, InputRole.TEMPLATE), (alkali, InputRole.CONDITION)), (persistent,))
    cip = graph.input("repair_enzyme", "CIP M0290 and T4 PNK M0201", MaterialKind.CATALYST)
    repaired = graph.product("repair_product", "dephosphorylated RNA fragments", MaterialKind.MOLECULAR,
                             state=molecule("RNA acceptor", (role("CMC-treated RNA fragment"), role("3′ OH acceptor end"))))
    graph.react("repair", "Dephosphorylate RNA ends with CIP and PNK", StepKind.DEPHOSPHORYLATION,
                ((persistent, InputRole.TEMPLATE), (cip, InputRole.CATALYST)), (repaired,))
    selected = _select(graph, "rna_size", repaired, "8% urea-TBE PAGE: source 80–100, 100–120 or 120–140-nt fraction",
                       "RNA outside the selected size fraction")
    adapter_state = Scene()
    adapter_state.strand("adapter", [PSEUDO_ADAPTER, PSEUDO_BLOCK], label="printed adapter", mod5="phosphate (printed; rApp stated in prose)", mod3="ddC")
    adapter = graph.input("ligate_adapter", "Printed /5Phos/ 3′ adapter; prose calls it preadenylated (unresolved)",
                          MaterialKind.ADAPTER, state=adapter_state, source=PSEUDO_SOURCE)
    ligase = graph.input("ligate_ligase", "T4 RNA ligase M0204; ATP-free buffer; 22 °C 2.5 h", MaterialKind.CATALYST)
    tagged_state = molecule("reported adapter-ligated RNA", (role("CMC-treated RNA fragment"), PSEUDO_ADAPTER, PSEUDO_BLOCK),
                  junctions=(("CMC-treated RNA fragment", PSEUDO_ADAPTER.name, "reported RNA–adapter junction"),))
    tagged = graph.product("ligate_product", "reported adapter-ligated RNA (donor activation unresolved)",
                           MaterialKind.MOLECULAR, state=tagged_state)
    graph.react("ligate", "Ligate the reported 3′ adapter (published activation conflict)", StepKind.LIGATION,
                ((selected, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)),
                (tagged,), geometry=AdapterLigationGeometry.unavailable("Methods call the donor preadenylated but print /5Phos/; "
                    "ATP-free ligation does not resolve which donor activation was actually used. Exact end state stays unresolved."))
    primer = graph.input("rt_primer", "5′-phosphorylated RT primer with two internal /iSp18/ spacers", MaterialKind.OLIGO,
                          state=rt_primer_state(), source=PSEUDO_SOURCE)
    enzyme = graph.input("rt_enzyme", "AMV-RT M5108, 42 °C 1 h, source RT buffer/dNTP/MgCl₂", MaterialKind.CATALYST)
    cdna_state = molecule("stop-defined single-strand cDNA", (*rt_primer_parts(),
                           role("cDNA extension ending at a Ψ-dependent stop")))
    cdna = graph.product("rt_product", "truncated Pseudo-seq cDNA", MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("rt", "Reverse-transcribe with AMV-RT", StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL), (enzyme, InputRole.CATALYST)), (cdna,))
    short = _select(graph, "cdna_size", cdna, "8% urea-TBE PAGE: recover truncated cDNA",
                    "cDNA and primer outside the recovered band")
    ligase_c = graph.input("circle_enzyme", "CircLigase CL4115K", MaterialKind.CATALYST)
    circle = graph.product("circle_product", "circularized stop cDNA (spacer/junction geometry not assembled)",
                            MaterialKind.MOLECULAR, state=molecule("circle boundary", (role("circularized spacer-containing cDNA"),)))
    graph.react("circle", "Circularize selected single-strand cDNA", StepKind.CIRCULARIZATION,
                ((short, InputRole.TEMPLATE), (ligase_c, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.unavailable("5′ phosphate and primer spacers are source known; "
                    "this model has not reconstructed the circular template and amplification path through its DNA components."))
    forward = graph.input("library_forward", "source PCR forward", MaterialKind.OLIGO,
                           state=Construct((PSEUDO_FWD,)), source=PSEUDO_SOURCE)
    reverse = graph.input("library_reverse", "source indexed PCR reverse (six-base sample index)", MaterialKind.OLIGO,
                           state=PSEUDO_REV_STATE, source=PSEUDO_SOURCE)
    pol = graph.input("library_polymerase", "Phusion M0530", MaterialKind.CATALYST)
    amplified = graph.product("library_amplified", "Pseudo-seq PCR products", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify circularized cDNA with Phusion and printed indexed primers", StepKind.PCR,
                ((circle, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (amplified,), geometry=PCRGeometry.unavailable(circle.state,
                    reason="Exact PCR primers are supplied, but the spacer-containing circular template and primer landing "
                           "have not been reconstructed; the final product cannot establish its own input geometry."))
    purified = _select(graph, "library_size", amplified, "Gel purify PCR products before sequencing",
                        "PCR products outside the selected library band")
    library = graph.product("library_library", "Pseudo-seq Illumina library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    precip = graph.input("library_recovery_condition", "ethanol precipitation of gel-purified PCR products", MaterialKind.CONDITION)
    graph.react("library_recovery", "Precipitate and recover the gel-purified library", StepKind.EXTRACTION,
                ((purified, InputRole.TEMPLATE), (precip, InputRole.CONDITION)), (library,), changes_molecule=False)
    graph.sequence(library, SequencingPlan(acquisition="unspecified", limitation=SEQUENCING_ENDING +
                    " Study-specific read length, read type and primer landing remain unassigned."))
    graph.declare_reviewed_scope(("input", "fragment_condition", "fragment", "cmc_condition", "cmc", "reverse_cmc_condition", "reverse_cmc",
        "repair_enzyme", "repair", "rna_size_gel", "rna_size", "ligate_adapter", "ligate_ligase", "ligate", "rt_primer", "rt_enzyme", "rt",
        "cdna_size_gel", "cdna_size", "circle_enzyme", "circle", "library_forward", "library_reverse", "library_polymerase", "library_pcr",
        "library_size_gel", "library_size", "library_recovery_condition", "library_recovery", "library_library"), PSEUDO_SOURCE)
    graph.annotate(DetailView("rt", "Printed RT-primer landing",
        "Computed pairing includes the adapter's terminal ddC nucleotide. The two /iSp18/ "
        "regions are non-nucleotide spacers; the donor activation discrepancy remains unresolved.",
        state=rt_priming_scene()))
    return graph.require_publishable()
