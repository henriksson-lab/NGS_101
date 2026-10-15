"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from rna_library_operations import molecule

DEFINITION = Definition("Nano3P-seq — Nanopore 3′-end capture", "10.1038/s41592-022-01714-w", "Begik et al., Nature Methods (2023)", "Template-switching reverse transcription captures native RNA 3′ ends for PCR-free Nanopore cDNA sequencing without RNA-adapter ligation.", "nanopore", ("motor adapter", "full-length cDNA from native 3-prime end"), "Nanopore traverses the cDNA from the captured native RNA 3′ end, retaining tail length and composition information.", "RNA", "The defining paper prints all three custom Nano3P-seq oligos and their modifications. The SQK-DCS109 Adapter Mix motor-adapter sequence remains proprietary to Oxford Nanopore and is marked inferred.")


from protocol_specs import workflow
from nano3p_chem import final_library as nano3p_final_library
from nano3p_chem import oligos as nano3p_oligos


def _source_workflow():
    return workflow(molecule("native RNA", ("RNA body", "native 3-prime end")), (
        ("Prime at the native RNA 3′ end", molecule("RNA–primer complex", ("RNA", "3-prime RT primer")),
         "Capture does not require a poly(A) tail or RNA-adapter ligation."),
        ("Template-switch during reverse transcription",
         molecule("full cDNA", ("template-switch handle", "RNA-derived cDNA"), joined=True),
         "The RT places a sequencing-compatible handle at the opposite cDNA end."),
        ("Attach the Nanopore sequencing adapter without PCR",
         molecule("Nanopore cDNA molecule", ("motor adapter", "full-length cDNA"), joined=True,
                  junctions=(("motor adapter", "full-length cDNA", "adapter ligation"),)),
         "One molecule reports RNA abundance and native 3′-tail composition."),
    ))

SPEC = get('nano3p-seq', definition=DEFINITION,
           workflow_fn=_source_workflow,
           final_override=nano3p_final_library(),
           oligos_fn=nano3p_oligos,
           sequencing_mechanism="nanopore_motor",
           motor_adapter_name="SQK-DCS109 AMX motor adapter")


def _local_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    rna = graph.input("rna", "native RNA with an unmodified 3-prime end",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("native RNA", ("RNA body", "native 3-prime end")),
                      primary=True)
    starter = graph.input("rt_starter", "preannealed blocked R_RNA / 5-prime-phosphorylated D_DNA starter duplex",
                          MaterialKind.OLIGO,
                          source=source + "; Methods print both oligos and modifications")
    rt = graph.input("reverse_transcriptase", "template-switching reverse transcriptase",
                     MaterialKind.CATALYST)
    full = graph.product("full_cdna", "TGIRT template-switched cDNA copied from the native RNA 3-prime end",
                         MaterialKind.MOLECULAR,
                         state=molecule("full cDNA", ("D_DNA handle", "RNA-derived cDNA")))
    graph.react("template_switch_rt", "Template-switch TGIRT from the starter duplex onto the native RNA 3-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((rna, InputRole.TEMPLATE), (starter, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (full,))
    compa = graph.input("compa", "published CompA_DNA ligation complement",
                        MaterialKind.OLIGO, source=source + "; Methods")
    ready = graph.product("ligation_ready", "CompA-annealed cDNA end",
                          MaterialKind.MOLECULAR, state=full.state)
    graph.react("anneal_compa", "Anneal CompA_DNA to make the adapter-ligation end",
                StepKind.HYBRIDIZATION,
                ((full, InputRole.TEMPLATE), (compa, InputRole.STRUCTURAL)), (ready,))
    motor = graph.input("motor_adapter", "Nanopore motor sequencing adapter",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "Nanopore adapter ligase", MaterialKind.CATALYST)
    library = graph.product("library", "PCR-free Nano3P-seq molecule", MaterialKind.LIBRARY,
                            state=final)
    graph.react("ligate_motor", "Attach the Nanopore motor adapter without PCR",
                StepKind.LIGATION,
                ((ready, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source defines the motor-adapter role; kit-revision-specific end chemistry is not asserted."))
    graph.sequence(library, SequencingPlan(limitation=d.sequencing))
    reviewed = ["rna", "rt_starter", "reverse_transcriptase", "template_switch_rt",
                "compa", "anneal_compa", "motor_adapter", "adapter_ligase",
                "ligate_motor", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY)

SPEC.protocol_graph = protocol_graph
