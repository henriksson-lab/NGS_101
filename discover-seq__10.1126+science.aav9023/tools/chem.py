"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import nebnext as nb
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from chemdraw import Construct, Segment, feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, closed_circle, da_tailed, molecule, role
from editing_outcomes import _finish_pcr

DISCOVER_AR001 = "ATCACG"


def discover_library() -> tuple[Construct, tuple]:
    """Historical single-index TruSeq ChIP product with a six-cycle index read."""
    return truseq_library([
        role("MRE11-bound genomic DNA from a Cas9 cut to a random sonication end"),
    ], "DISCOVER-Seq TruSeq ChIP library", dual_index=False, i7_nt=6)


def workflow():
    return chain(molecule("living cells or tissue receiving a CRISPR nuclease", (role("chromosomal on-target and off-target sites"),)), (
        ("Let MRE11 bind newly made double-strand breaks",
         molecule("edited chromatin with repair complexes", (role("Cas9 cleavage end"), role("MRE11 at the DNA end"), role("surrounding chromatin"))),
         "The sampling time is optimized for the delivery system; the detailed protocol suggests approximately 8–16 hours after RNP delivery."),
        ("Formaldehyde-crosslink cells, lyse nuclei and sonicate chromatin",
         molecule("150–500-bp crosslinked chromatin fragments", (role("random sonication end"), role("MRE11-bound DNA reaching a nuclease cut"), role("Cas9 cleavage end"))),
         "Crosslinking preserves the transient MRE11–DNA complex while mechanical shearing makes ChIP-sized fragments."),
        ("Immunoprecipitate MRE11 on Protein A magnetic beads",
         molecule("MRE11-enriched chromatin", (role("antibody–MRE11 complex"), role("DNA adjacent to a cellular DSB"), role("Cas9 cleavage boundary"))),
         "An unedited or non-targeting sample measures background; ChIP enrichment is checked by qPCR before library construction."),
        ("Reverse crosslinks and purify the ChIP DNA",
         molecule("purified break-adjacent DNA", (role("random sonication end"), role("MRE11 footprint region"), role("native nuclease-cleavage end"))),
         "Proteinase K and RNase A remove the captured protein and RNA while retaining the DNA endpoints used by BLENDER."),
        ("End-repair and add one 3-prime dA",
         molecule("A-tailed ChIP fragments", (role("end-repaired random end"), role("break-adjacent insert"), role("end-repaired cleavage end"), role("3-prime dA"))),
         "The historical TruSeq ChIP kit blunts overhangs and then adenylates the repaired fragments."),
        ("Ligate a six-base indexed TruSeq ChIP adapter",
         molecule("adapter-ligated ChIP library", (role("TruSeq Read 1 side"), role("break-adjacent genomic insert"), role("TruSeq Read 2 side with six-base i7 index"))),
         "The discontinued kit supplies single-index adapters; Table 4/5 instructs recording six index bases in the sample sheet."),
        ("PCR-enrich, remove adapter dimers and pool",
         molecule("single-index DISCOVER-Seq library", (role("P5 and Read 1 arm"), role("MRE11-enriched insert"), role("Read 2 arm, i7 and P7"))),
         "Amplification completes the flow-cell ends; the detailed protocol adds a 0.8× bead cleanup when adapter dimers remain."),
    ))


def oligos():
    return (
        f"TruSeq ChIP AR001 sample index — 5′-{DISCOVER_AR001}-3′",
        f"TruSeq universal PCR primer — 5′-{il.TRUSEQ_P5_FULL}-3′",
        f"TruSeq AR001 indexed P7 primer model — 5′-{il.P7}{DISCOVER_AR001}{il.TRUSEQ_READ2}-3′",
    )


TITLE = 'DISCOVER-Seq — MRE11 ChIP-seq at CRISPR cleavage sites'
NOTES = '01_discover_seq.html'
SOURCE = '<a href="https://doi.org/10.1126/science.aav9023">Wienert et al., Science (2019)</a>; stepwise implementation from the authors’ <a href="https://doi.org/10.1038/s41596-020-0309-5">Nature Protocols procedure</a> and the named TruSeq ChIP kit guide.'
SUMMARY = 'MRE11 recruited to newly cut chromosomal DNA is fixed in place and immunoprecipitated; ordinary paired-end ChIP-seq then preserves the sharp genomic endpoints that BLENDER recognizes as nuclease cuts.'
CAVEAT = 'DISCOVER-Seq does not use a custom capture oligo or UMI. Its defining implementation delegates library construction to the now-discontinued single-index TruSeq ChIP kit; index identities are public, but the kit did not publish every supplied oligo as an ordering sequence.'
FINAL_CAPTION = 'Either genomic orientation can be captured. One library strand is shown from a random sonication boundary through MRE11-enriched DNA to the cellular nuclease-cleavage boundary; the six-base i7 sample index comes from the historical TruSeq ChIP adapter.'
READ_LENGTHS = {'Read 1': 100, 'Index 1 (i7)': 6, 'Read 2': 100}
SEQUENCING_INTRO = 'The detailed protocol specifies 100-bp paired-end sequencing. Standard TruSeq Read 1, Index 1 and Read 2 primers bind the final library; paired fragment endpoints, not a barcode read, provide the cut-site signal.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = discover_library()

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    chromatin = graph.input("chromatin", "edited cells or tissue with nuclease-created breaks",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("edited chromatin", (role("chromosomal nuclease cuts"),)),
                            primary=True)
    mre11 = graph.input("mre11_antibody", "MRE11-binding antibody probe",
                        MaterialKind.SELECTION)
    bound = graph.product("mre11_bound", "MRE11-bound break chromatin",
                          MaterialKind.MOLECULAR, state=chromatin.state)
    unbound = graph.product("mre11_unbound", "chromatin without captured MRE11-bound ends",
                            MaterialKind.WASTE)
    graph.react("bind_mre11", "Recognize MRE11 at newly made DNA breaks", StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (mre11, InputRole.SELECTION)),
                (bound, unbound), changes_molecule=False)
    graph.discard(unbound, "does not carry the MRE11-bound break-end signal")
    formaldehyde = graph.input("formaldehyde", "formaldehyde cross-linking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "cross-linked MRE11–break chromatin",
                          MaterialKind.MOLECULAR, state=bound.state)
    graph.react("crosslink", "Cross-link MRE11–DNA complexes", StepKind.CROSS_LINKING,
                ((bound, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    shear = graph.input("sonication", "chromatin sonication condition", MaterialKind.CONDITION)
    fragments = graph.product("fragments", "150–500-bp cross-linked chromatin fragments",
                              MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("sonicate", "Sonicate fixed chromatin", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    beads = graph.input("protein_a_beads", "Protein A magnetic beads", MaterialKind.SELECTION)
    captured = graph.product("captured", "MRE11-enriched chromatin fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("chip_flowthrough", "MRE11-depleted ChIP flow-through",
                             MaterialKind.WASTE)
    graph.react("chip", "Immunoprecipitate MRE11-bound chromatin", StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "is not retained by MRE11 immunoprecipitation")
    reversal = graph.input("crosslink_reversal", "cross-link reversal and protein-removal condition",
                           MaterialKind.CONDITION)
    dna = graph.product("chip_dna", "purified break-adjacent ChIP DNA",
                        MaterialKind.MOLECULAR, state=captured.state)
    graph.react("purify_dna", "Reverse cross-links and remove protein",
                StepKind.DEPROTEINIZATION,
                ((captured, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (dna,))
    repair = graph.input("repair_enzyme", "TruSeq ChIP end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired ChIP fragments",
                             MaterialKind.MOLECULAR, state=dna.state)
    graph.react("repair", "End-repair ChIP DNA", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed ChIP fragments",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("tail", "Add one 3-prime dA", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("truseq_adapter", "six-base indexed TruSeq ChIP adapter",
                          MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "TruSeq adapter ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated", "adapter-ligated ChIP fragments",
                            MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate", "Ligate the indexed TruSeq ChIP adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "TruSeq ChIP adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The named TruSeq ChIP workflow ligates its indexed adapter to dA-tailed fragments."))
    reviewed = ["chromatin", "mre11_antibody", "bind_mre11", "formaldehyde", "crosslink",
                "sonication", "sonicate", "protein_a_beads", "chip", "crosslink_reversal",
                "purify_dna", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "truseq_adapter", "adapter_ligase", "ligate"]
    _finish_pcr(graph, "library", ligated, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label="DISCOVER-Seq library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
