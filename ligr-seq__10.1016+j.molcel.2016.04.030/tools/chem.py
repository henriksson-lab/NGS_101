"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from interaction_capture import ReversibleDuplexCrosslink
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_rt_pcr
from rna_library_operations import molecule

DEFINITION = Definition(
    "LIGR-seq — ligation of interacting RNA", "10.1016/j.molcel.2016.04.030",
    "Sharma et al., Molecular Cell (2016)",
    "AMT-crosslinked RNA duplexes survive RNase R enrichment and are proximity-ligated before cDNA library construction.",
    "truseq", ("5-prime random barcode", "RNA A-derived cDNA", "RNA B-derived cDNA"),
    "The study analyzed single-end, strand-specific reads after trimming a 5-prime random barcode. The cited SMARTer kit is described as Illumina-compatible; this page infers standard sequencing-primer endpoints, while the exact kit priming junctions and selected index remain unverified.",
    "RNA",
    "The article supplement names a SMARTer Stranded RNA-Seq kit generically and describes reaction roles, without giving a catalogue number, kit revision, selected index, or RT/template-switching primer bases. The 2013 Clontech file is only a product components list. An online transcript of the 2014 Clontech User Manual (revision 050714, pp. 6–8) labels its SMARTer oligo and N6 primer proprietary and prints twelve six-base index values, but its original vendor PDF is not available in the source bundle. The study's exact kit generation and selected index remain unverified."
)


from protocol_specs import psoralen_workflow


def _source_workflow():
    return psoralen_workflow('AMT', 'RNase R depletes uncross-linked linear RNA')

SPEC = get('ligr-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers, *, reagent: str,
                   selection_label: str):
    """Build a psoralen contact graph from source-specific reagent choices."""
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = ReversibleDuplexCrosslink(reagent, selection_label)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA molecules",)), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " and long-wave UV condition",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR,
                                state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    fragmentation = graph.input("fragmentation_condition", "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA", StepKind.FRAGMENTATION,
                ((crosslinked, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragments,))
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    ligase = graph.input("proximity_ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((selected, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("photoreversal", reagent + " photoreversal condition",
                           MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    reviewed = ["rna", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "selector", "select_duplexes", "proximity_ligase", "proximity_ligation",
                "photoreversal", "photoreverse"]
    finish_rt_pcr(graph, "library", reversed_chimera, final, primers, d.sequencing,
                   reviewed, label=d.title)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS, reagent='AMT psoralen', selection_label='RNase-R-resistant duplex selection')

SPEC.protocol_graph = protocol_graph
