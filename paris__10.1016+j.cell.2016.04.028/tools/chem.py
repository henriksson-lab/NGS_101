"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from interaction_capture import ReversibleDuplexCrosslink
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from rna_library_operations import molecule

DEFINITION = Definition(
    "PARIS — psoralen analysis of RNA interactions and structures", "10.1016/j.cell.2016.04.028",
    "Lu et al., Cell (2016)",
    "AMT photocrosslinking, two-dimensional gel enrichment and proximity ligation preserve the two partners of an RNA duplex in one chimera.",
    "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"),
    "The original 2016 supplement distinguishes HeLa random-hexamer RT from HEK293/mES adapter-primed RT, which inherits Flynn 2016 icSHAPE library preparation. The circular-cDNA junction still needs molecular reconstruction before primer placement can be computed.",
    "RNA",
    "The original 2016 Extended Experimental Procedures specify two library routes. The graph follows the HEK293/mES adapter-ligation route inherited from Flynn 2016; the separate HeLa random-hexamer primer is shown below. The additionally labelled 2018 oligos document a later implementation with different RT-primer random/index tracts. Exact circular-cDNA geometry is not yet reconstructed."
)


from protocol_specs import psoralen_workflow
from paris_chem import oligos as paris_oligos


def _source_workflow():
    return psoralen_workflow('AMT', 'cross-linked duplexes selected above the diagonal on a 2D gel')

SPEC = get('paris', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=lambda: (
               '<b>2016 HeLa random-hexamer RT primer</b> — 5′-/5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN-3′ (D = A/G/T; XXXX = four-base sample barcode; first N6 = duplicate tag; terminal N6 = RNA priming)',
               '<b>2016 HEK293/mES inherited RT-primer-1 example</b> — 5′-/5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGGA/iSp18/GGATCC/iSp18/TACTGAACCGC-3′ (Flynn 2016 icSHAPE reagent list; AACC is its printed example index)',
               '<b>Later implementation: Lu 2018 oligos</b> — these are version-labelled comparison reagents, not an assertion that the 2016 RT primers were identical.',
               *paris_oligos()),
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
    adapter = graph.input("rna_adapter", "2016 HEK293/mES barcoded RNA-library adapter", MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "RNA adapter ligase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_rna", "adapter-tagged PARIS RNA chimera", MaterialKind.MOLECULAR,
                           state=molecule("adapter-tagged chimera", ("RNA chimera", "RNA library adapter")))
    graph.react("ligate_adapter", "Ligate the RNA-library adapter (2016 HEK293/mES route)", StepKind.LIGATION,
                ((reversed_chimera, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The original supplement inherits Flynn 2016 step 26 onward; adapter chemistry is source-located, but this graph has not reconstructed the ligation-end Scene."))
    rt_primer = graph.input("library_rt_primer", "Flynn 2016 adapter-annealing RT-primer-1", MaterialKind.OLIGO)
    rt = graph.input("library_rt_enzyme", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "PARIS first-strand cDNA", MaterialKind.MOLECULAR,
                         state=molecule("PARIS cDNA", ("RT primer and chimera cDNA",)))
    graph.react("library_rt", "Reverse-transcribe the adapter-tagged RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    circle_ligase = graph.input("circle_ligase", "CircLigase II", MaterialKind.CATALYST)
    circle = graph.product("cdna_circle", "circularized PARIS cDNA", MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("circularize_cdna", "Circularize the cDNA before library PCR", StepKind.CIRCULARIZATION,
                ((cdna, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known("single-strand intramolecular ligation",
                    statement="Inherited Flynn 2016 icSHAPE step 58 circularizes purified cDNA with CircLigase II; exact circular junction is not yet reconstructed in this model."))
    reviewed.extend(("rna_adapter", "adapter_ligase", "ligate_adapter", "library_rt_primer",
                     "library_rt_enzyme", "library_rt", "circle_ligase", "circularize_cdna"))
    finish_pcr(graph, "library", circle, final, primers, d.sequencing, reviewed, label=d.title)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS, reagent='AMT psoralen', selection_label='two-dimensional-gel duplex selection')

SPEC.protocol_graph = protocol_graph
