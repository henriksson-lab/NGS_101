"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from nanoseq_chem import nanoseq_restrictive_end_repair
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("NanoSeq — nanorate duplex sequencing", "10.1038/s41586-021-03477-4", "Abascal et al., Nature (2021)", "Restriction fragmentation and dideoxynucleotide-based restrictive end repair prevent lesion copying before duplex consensus sequencing.", "truseq", ("duplex tag A", "genomic insert", "duplex tag B"), "Paired-end Illumina reads are grouped into complementary strand families and then a duplex consensus.", "DNA")


from protocol_specs import workflow
from nanoseq_chem import nanoseq_restrictive_end_repair


def _source_workflow():
    repaired = nanoseq_restrictive_end_repair()
    return workflow(molecule("high-molecular-weight genomic DNA", ("genomic DNA",)), (
        ("Digest with a restriction enzyme instead of random shearing", repaired,
         "Defined intact ends limit lesion-generating end repair."),
        ("Perform restrictive end repair with dideoxynucleotides", repaired,
         "Nick extension is terminated before it can copy damage into both strands."),
        ("Ligate duplex adapters, amplify and recover both strand families",
         molecule("duplex-consensus library", ("duplex tag A", "insert", "duplex tag B"), joined=True,
                  junctions=(("duplex tag A", "insert", "adapter ligation"),
                             ("insert", "duplex tag B", "adapter ligation"))),
         "Complementary tags associate Watson and Crick read families."),
    ))

SPEC = get('nanoseq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("high-molecular-weight genomic DNA", ("genomic DNA",)),
                      primary=True)
    restriction = graph.input("restriction_enzyme", "restriction endonuclease",
                              MaterialKind.CATALYST)
    fragments = graph.product("fragments", "restriction-defined genomic fragments",
                              MaterialKind.MOLECULAR, state=nanoseq_restrictive_end_repair())
    graph.react("digest", "Restriction-fragment genomic DNA", StepKind.DIGESTION,
                ((dna, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)), (fragments,))
    repair = graph.input("repair_mix", "dideoxynucleotide restrictive end-repair mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "restrictively end-repaired fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("restrictive_end_repair", "Perform restrictive end repair",
                StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapters = graph.input("duplex_adapters", "complementary duplex-family adapters",
                           MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "duplex-adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "duplex-tagged genomic fragments", MaterialKind.MOLECULAR,
                            state=molecule("duplex-consensus precursor", (
                                "duplex tag A", "genomic insert", "duplex tag B"),
                                junctions=(("duplex tag A", "genomic insert", "adapter ligation"),
                                           ("genomic insert", "duplex tag B", "adapter ligation"))))
    graph.react("ligate", "Ligate duplex-family adapters", StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The source establishes complementary duplex tags; the role-level page does not assert one adapter-end sequence."))
    reviewed = ["dna", "restriction_enzyme", "digest", "repair_mix",
                "restrictive_end_repair", "duplex_adapters", "ligase", "ligate"]
    finish_pcr(graph, "library", adapted, final, primers, d.sequencing, reviewed,
                label="NanoSeq duplex")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
