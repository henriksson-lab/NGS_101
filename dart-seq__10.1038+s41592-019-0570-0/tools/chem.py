"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_rt_pcr
from rna_library_operations import molecule

DEFINITION = Definition("DART-seq — antibody-free m6A mapping", "10.1038/s41592-019-0570-0", "Meyer, Nature Methods (2019)", "APOBEC1 fused to the YTH m6A-binding domain writes nearby C-to-U marks that become identifiable mismatches in an RNA-seq library.", "truseq", ("cDNA carrying m6A-adjacent C-to-T mark",), "Standard RNA-seq reads report C-to-T conversions adjacent to YTH-bound m6A sites.", "RNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("cellular RNA", ("RNA containing m6A",)), (
        ("Express APOBEC1 fused to the YTH m6A-binding domain",
         molecule("targeted RNA", ("APOBEC1–YTH", "m6A-adjacent cytidine")),
         "YTH recruits cytidine deamination near m6A."),
        ("Convert nearby cytidine to uridine in the cell",
         molecule("marked RNA", ("RNA", "C-to-U mark near m6A")),
         "The modification is read as a C-to-T mismatch after reverse transcription."),
        ("Prepare an ordinary RNA-seq library",
         molecule("DART-seq library", ("left adapter", "cDNA carrying C-to-T mark", "right adapter"), joined=True,
                  junctions=(("left adapter", "cDNA carrying C-to-T mark", "adapter ligation"),
                             ("cDNA carrying C-to-T mark", "right adapter", "adapter ligation"))),
         "The library reports an enzymatic mark, not direct sequencing of m6A."),
    ))

SPEC = get('dart-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "cellular RNA containing m6A sites", MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA containing m6A",)), primary=True)
    fusion = graph.input("fusion", "APOBEC1–YTH m6A-reader fusion", MaterialKind.SELECTION)
    targeted = graph.product("targeted", "RNA bound by APOBEC1–YTH near m6A",
                             MaterialKind.MOLECULAR, state=rna.state)
    unbound = graph.product("unbound", "RNA not selected by the reader fusion", MaterialKind.WASTE)
    graph.react("bind", "Recruit APOBEC1 to m6A-proximal RNA",
                StepKind.BINDING,
                ((rna, InputRole.TEMPLATE), (fusion, InputRole.SELECTION)),
                (targeted, unbound))
    graph.discard(unbound, "does not receive the targeted deamination mark")
    deamination = graph.input("deamination_condition", "APOBEC cytidine-deamination condition",
                              MaterialKind.CONDITION)
    marked = graph.product("marked", "RNA carrying m6A-adjacent C-to-U marks",
                           MaterialKind.MOLECULAR,
                           state=molecule("marked RNA", ("RNA", "C-to-U mark near m6A")))
    graph.react("deaminate", "Convert nearby cytidine to uridine",
                StepKind.CHEMICAL_CONVERSION,
                ((targeted, InputRole.TEMPLATE), (deamination, InputRole.CONDITION)),
                (marked,))
    reviewed = ["rna", "fusion", "bind", "deamination_condition", "deaminate"]
    finish_rt_pcr(graph, "library", marked, final, primers, d.sequencing, reviewed,
                   label="DART-seq marked RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
