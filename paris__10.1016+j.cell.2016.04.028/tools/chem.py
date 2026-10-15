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

DEFINITION = Definition("PARIS — psoralen analysis of RNA interactions and structures", "10.1016/j.cell.2016.04.028", "Lu et al., Cell (2016)", "AMT photocrosslinking, two-dimensional gel enrichment and proximity ligation preserve the two partners of an RNA duplex in one chimera.", "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"), "The detailed PARIS protocol publishes the complete custom P6 sequencing primer and all adapter, RT and PCR oligos. Its circular-cDNA junction still needs to be represented before primer placement can be computed rather than stated by role.", "RNA", "The exact oligos were not absent: Lu et al.'s detailed PARIS protocol prints them, including modifications and the 24-member RT-index list. The current diagram retains a role-level final construct until the circular-cDNA library geometry is reconstructed.")


from protocol_specs import psoralen_workflow
from paris_chem import oligos as paris_oligos


def _source_workflow():
    return psoralen_workflow('AMT', 'cross-linked duplexes selected above the diagonal on a 2D gel')

SPEC = get('paris', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=paris_oligos,
           sequencing_unavailable='The custom sequencing-primer bases are published; exact placement remains unresolved only because the circular-cDNA junction is not yet represented.')


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
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS, reagent='AMT psoralen', selection_label='two-dimensional-gel duplex selection')

SPEC.protocol_graph = protocol_graph
