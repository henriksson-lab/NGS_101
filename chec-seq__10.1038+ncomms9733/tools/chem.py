"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, da_tailed, molecule, role
from protein_dna_mapping import henikoff_single_index_library, _finish_henikoff_library

CHEC_UNIVERSAL_ADAPTER = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

CHEC_INDEXED_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

CHEC_P5_PCR = "AATGATACGGCGACCACCGAG"

CHEC_P7_PCR = "CAAGCAGAAGACGGCATACGAG"

HENIKOFF_DA_LIGATION = reported_adapter_end(
    "phosphorothioate-protected TruSeq adapter", kind="cohesive",
    polarity="3-prime", terminal_sequence="T",
    statement=("The cited Henikoff low-input procedure publishes the protected TruSeq "
               "adapter and its terminal dT/dA ligation compatibility."))


def chec_library():
    return henikoff_single_index_library(
        "ChEC-seq library", (role("DNA released by two nearby protein-tethered MNase cuts"),))


def workflow():
    return chain(molecule("living budding yeast expressing a chromatin-protein–MNase fusion", (
        role("chromatin"), role("DNA-bound protein fused to 3xFLAG-MNase"))), (
        ("Wash cells and permeabilize for 5 minutes with 0.1% digitonin",
         molecule("permeabilized cells with inactive tethered MNase", (
             role("accessible intracellular chromatin"), role("protein-tethered MNase"))),
         "Low nuclear calcium keeps MNase largely inactive before the timed reaction. A matched soluble nuclear MNase strain measures accessibility and enzyme-sequence background."),
        ("Add 2 mM calcium and remove timed digestion aliquots",
         molecule("chromatin cut near occupied protein sites", (
             role("distal chromatin"), role("MNase cleavage end"), role("protein-protected binding region"),
             role("MNase cleavage end"), role("distal chromatin"))),
         "Calcium activates the fusion in situ. Closely spaced seconds-to-minutes samples distinguish rapid tether-proximal cleavage from slower mass-action background."),
        ("Quench each aliquot with EDTA/EGTA, SDS and high salt",
         molecule("stopped ChEC digest", (
             role("large background genomic DNA"), role("short two-cut protein-site fragment"))),
         "Chelators remove the essential calcium and SDS disrupts cells and protein complexes, freezing each kinetic endpoint."),
        ("Extract nucleic acid, digest RNA and retain the small-DNA supernatant",
         molecule("size-selected ChEC DNA", (role("short fragment bounded by two MNase cleavage ends"),)),
         "At a 2.5:1 AMPure-bead:sample ratio the large genomic DNA binds; the unbound fraction is recovered, enriching short protected fragments before library construction."),
        ("End-repair, phosphorylate and add one 3-prime dA",
         molecule("A-tailed ChEC fragment", (
             role("repaired MNase boundary"), role("protein-site insert"), role("repaired MNase boundary"), role("3-prime dA"))),
         "T4 DNA polymerase, Klenow and T4 PNK make ligatable blunt phosphorylated ends; Klenow exo-minus then creates the dA overhang."),
        ("Ligate the six-base indexed TruSeq adapter and PCR-enrich",
         molecule("single-index paired-end ChEC-seq library", (
             role("P5 and Read 1 arm"), role("two-cut ChEC insert"), role("Read 2 arm, six-base i7 and P7"))),
         "The low-input paired-end protocol preserves short inserts and uses KAPA polymerase with the published P5 and P7 enrichment primers."),
    ))


def oligos():
    return (
        "TruSeq Universal Adapter — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        "TruSeq Indexed Adapter — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′",
        f"P5 library-enrichment primer — 5′-{CHEC_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{CHEC_P7_PCR}-3′",
    )


FINAL_LIBRARY, SEQ_PRIMERS = chec_library()
TITLE = 'ChEC-seq — calcium-triggered chromatin endogenous cleavage sequencing'
NOTES = '01_chec_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/ncomms9733">Zentner et al., Nature Communications (2015)</a>; library oligos from the cited Henikoff laboratory low-input paired-end protocol.'
SUMMARY = 'A protein of interest is genetically fused to MNase. Digitonin admits calcium into living yeast, activating cleavage beside occupied binding sites; short two-cut fragments are selected and converted into a paired-end Illumina library.'
CAVEAT = 'ChEC-seq reports calcium-triggered cleavage by a mobile protein–MNase fusion, not an antibody-captured footprint. A free nuclear MNase control and early kinetic samples are essential because slower cleavage contains accessibility, enzyme-bias and transient-collision signal.'
FINAL_CAPTION = 'The selected insert is bounded by two MNase-created ends near a protein-binding event. End repair and dA-tailing erase the original overhang geometry but retain both genomic boundary coordinates; the historical construct has one six-base i7 index and no i5 index.'
READ_LENGTHS = {'Read 1': 25, 'Index 1 (i7)': 6, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining data set used 25-base paired-end reads. Standard TruSeq Read 1 and Read 2 primers approach the short ChEC insert from opposite ends, and the six-cycle Index 1 read identifies the sample.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = ("Zentner et al., Nature Communications (2015); doi:10.1038/ncomms9733; "
              "library chemistry: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ChEC-seq", authority=source)
    cells = graph.input(
        "fusion_cells", "living yeast expressing the chromatin-protein–MNase fusion",
        MaterialKind.BIOLOGICAL,
        state=molecule("living budding yeast expressing a protein–MNase fusion", (
            role("chromatin"), role("DNA-bound protein fused to 3xFLAG-MNase"))),
        primary=True)
    digitonin = graph.input("digitonin", "0.1% digitonin permeabilization condition",
                            MaterialKind.CONDITION)
    chromatin = graph.product(
        "permeabilized_chromatin", "accessible intracellular chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("permeabilized cells", (role("accessible intracellular chromatin"),)))
    tethered_mnase = graph.product("tethered_mnase", "inactive protein-tethered MNase",
                                   MaterialKind.CATALYST)
    graph.react("permeabilize", "Permeabilize cells with digitonin",
                StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE), (digitonin, InputRole.CONDITION)),
                (chromatin, tethered_mnase))
    calcium = graph.input("calcium", "2 mM calcium activation condition",
                          MaterialKind.CONDITION)
    cut = graph.product(
        "cut_chromatin", "chromatin cut near occupied protein sites",
        MaterialKind.MOLECULAR,
        state=molecule("chromatin cut near occupied protein sites", (
            role("distal chromatin"), role("MNase cleavage end"),
            role("protein-protected binding region"), role("MNase cleavage end"),
            role("distal chromatin"))))
    graph.react("timed_digest", "Activate tethered MNase and remove timed aliquots",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (tethered_mnase, InputRole.CATALYST),
                 (calcium, InputRole.CONDITION)), (cut,))
    quench = graph.input("quench_reagents", "EDTA/EGTA, SDS and high-salt quench",
                         MaterialKind.CONDITION)
    stopped = graph.product(
        "stopped_digest", "stopped ChEC digest", MaterialKind.MOLECULAR,
        state=molecule("stopped ChEC digest", (
            role("large background genomic DNA"), role("short two-cut protein-site fragment"))))
    graph.react("quench", "Quench each timed digestion aliquot", StepKind.QUENCHING,
                ((cut, InputRole.TEMPLATE), (quench, InputRole.CONDITION)), (stopped,))
    extraction = graph.input("extraction_condition", "nucleic-acid extraction and RNase condition",
                             MaterialKind.CONDITION)
    extracted = graph.product("extracted_dna", "extracted ChEC DNA",
                              MaterialKind.MOLECULAR, state=stopped.state)
    graph.react("extract", "Extract nucleic acid and digest RNA", StepKind.EXTRACTION,
                ((stopped, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (extracted,))
    beads = graph.input("ampure_selection", "2.5:1 AMPure size-selection rule",
                        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_fragments", "short ChEC DNA fragments", MaterialKind.MOLECULAR,
        state=molecule("size-selected ChEC DNA", (
            role("short fragment bounded by two MNase cleavage ends"),)))
    large = graph.product("large_genomic_dna", "large bead-bound genomic DNA",
                          MaterialKind.WASTE)
    graph.react("size_select", "Retain the small-DNA supernatant", StepKind.SELECTION,
                ((extracted, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, large))
    graph.discard(large, "binds the AMPure beads while the short ChEC fragments remain in supernatant")
    reviewed = [
        "fusion_cells", "digitonin", "permeabilize", "calcium", "timed_digest",
        "quench_reagents", "quench", "extraction_condition", "extract",
        "ampure_selection", "size_select",
    ]
    _finish_henikoff_library(
        graph, "chec", selected, "ChEC", final, primers, source, reviewed,
        {"Read 1": 25, "Index 1 (i7)": 6, "Read 2": 25},
        HENIKOFF_DA_LIGATION, CHEC_P5_PCR, CHEC_P7_PCR)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
