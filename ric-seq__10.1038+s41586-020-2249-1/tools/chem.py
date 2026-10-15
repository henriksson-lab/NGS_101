"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from interaction_capture import FixedRNPContact
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_rt_pcr
from rna_library_operations import molecule

DEFINITION = Definition("RIC-seq — RNA in situ conformation sequencing", "10.1038/s41586-020-2249-1", "Cai et al., Nature (2020)", "Crosslinked RNA–protein assemblies are fragmented and proximity-ligated in situ to capture global spatial RNA–RNA contacts.", "custom", ("RNA A cDNA", "RNA B cDNA"), "Paired-end reads identify the two RNA fragments joined at each in-situ ligation junction.", "RNA")


from protocol_specs import workflow
from interaction_capture import FixedRNPContact


def _source_workflow():
    contact = FixedRNPContact("formaldehyde crosslinking")
    return workflow(molecule("fixed cells", ("spatially proximal RNAs",)), (
        ("Cross-link RNA–protein assemblies in situ and fragment RNA",
         contact.complex(),
         "Formaldehyde preserves spatial contacts, including protein-mediated ones."),
        ("Label RNA ends and proximity-ligate them in situ",
         contact.chimera(),
         "The ligation junction records a spatial RNA–RNA contact."),
        ("Reverse crosslinks, enrich chimeras and make a paired-end library",
         molecule("RIC-seq library", ("left adapter", "RNA A cDNA", "RNA B cDNA", "right adapter"), joined=True,
                  junctions=(("left adapter", "RNA A cDNA", "adapter ligation"),
                             ("RNA A cDNA", "RNA B cDNA", "copied proximity-ligation junction"),
                             ("RNA B cDNA", "right adapter", "adapter ligation"))),
         "Paired sequence identifies both sides of the RNA contact."),
    ))

SPEC = get('ric-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = FixedRNPContact("formaldehyde crosslinking")
    cells = graph.input("cells", "cells containing spatially proximal RNA–protein assemblies",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("fixed-cell input", ("spatially proximal RNAs",)), primary=True)
    formaldehyde = graph.input("crosslinker", "formaldehyde crosslinking condition",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "crosslinked RNA–protein neighbourhoods",
                          MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("crosslink", "Crosslink RNA–protein assemblies in situ",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    fragmentation = graph.input("fragmentation_condition", "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "fragmented fixed RNA neighbourhoods",
                              MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("fragment", "Fragment RNA in the fixed assemblies", StepKind.FRAGMENTATION,
                ((fixed, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)), (fragments,))
    ligase = graph.input("proximity_ligase", "RNA proximity ligase", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "in-situ proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("ligate", "Proximity-ligate labelled RNA ends in situ",
                StepKind.PROXIMITY_LIGATION,
                ((fragments, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("reversal", "crosslink-reversal and protein-removal condition",
                           MaterialKind.CONDITION)
    purified = graph.product("purified", "de-crosslinked RNA chimeras",
                             MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("reverse_crosslinks", "Reverse crosslinks and recover RNA chimeras",
                StepKind.DEPROTEINIZATION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    reviewed = ["cells", "crosslinker", "crosslink", "fragmentation_condition", "fragment",
                "proximity_ligase", "ligate", "reversal", "reverse_crosslinks"]
    finish_rt_pcr(graph, "library", purified, final, primers, d.sequencing, reviewed,
                   label="RIC-seq chimera")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
