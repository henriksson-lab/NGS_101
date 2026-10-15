"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from adapter_ligation import reported_adapter_end
from duplex_error import CODEC_ADAP5_1
from duplex_error import CODEC_ADAP7_1
from duplex_error import codec_linked_duplex
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("CODEC — single-duplex error correction", "10.1038/s41588-023-01376-0", "Bae et al., Nature Genetics (2023)", "An adapter quadruplex and strand-displacing extension place readable copies of both strands of one original DNA duplex into a single Illumina cluster molecule.", "custom", ("Watson copy", "CODEC linkage", "Crick copy"), "Standard TruSeq Read 1 and Read 2 primers bind within the two published bridge oligos and read outward through a 3-nt UMI, an 18-nt inline sample index and one original strand copy.", "DNA", "Supplementary Table 5 publishes all four adapter oligos in eight matched sample-index sets. Set 1 is shown as a representative; the linked post-extension product remains role-scaled where the source figure does not print one contiguous strand.")


from protocol_specs import workflow
from duplex_error import codec_linked_duplex
from duplex_error import codec_final_library
from duplex_error import codec_oligos


def _source_workflow():
    return workflow(molecule("original genomic duplex", ("Watson strand", "Crick strand")), (
        ("End-repair and ligate the preassembled CODEC adapter quadruplex",
         molecule("adapter-bound duplex", ("left CODEC adapter", "original duplex", "right CODEC adapter"), joined=True,
                  junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                             ("original duplex", "right CODEC adapter", "adapter ligation"))),
         "The four adapter oligos, modifications and matched index sets are printed in Supplementary Table 5."),
        ("Strand-displacing extension copies and physically links both complements",
         codec_linked_duplex(), "Both original strands become recoverable in one cluster molecule."),
        ("PCR and outward Illumina sequencing", codec_linked_duplex(),
         "Reads recover insert, UMI and index information for within-pair consensus."),
    ))

SPEC = get('codec', definition=DEFINITION,
           workflow_fn=_source_workflow,
           final_override=codec_final_library(),
           primers_override="standard",
           oligos_fn=codec_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    duplex = graph.input("duplex", "original genomic Watson/Crick duplex",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("original genomic duplex", ("Watson strand", "Crick strand")),
                         primary=True)
    repair = graph.input("end_repair_enzyme", "end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired genomic duplex", MaterialKind.MOLECULAR,
                             state=duplex.state)
    graph.react("end_repair", "End-repair the original duplex", StepKind.END_REPAIR,
                ((duplex, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapter = graph.input("codec_adapter", "preassembled four-oligo CODEC adapter quadruplex",
                          MaterialKind.ADAPTER,
                          source=source + "; Supplementary Table 5")
    ligase = graph.input("ligase", "CODEC adapter ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "CODEC-adapter-bound duplex", MaterialKind.MOLECULAR,
                            state=molecule("adapter-bound duplex", (
                                "left CODEC adapter", "original duplex", "right CODEC adapter"),
                                junctions=(("left CODEC adapter", "original duplex", "adapter ligation"),
                                           ("original duplex", "right CODEC adapter", "adapter ligation"))))
    graph.react("ligate", "Ligate the preassembled CODEC adapter quadruplex",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive",
                    *reported_adapter_end(
                        "LD4-adap5-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP5_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    *reported_adapter_end(
                        "LD4-adap7-1 arm", kind="cohesive", polarity="3-prime",
                        terminal_sequence="T", body_sequence=CODEC_ADAP7_1[:-1],
                        statement="Supplementary Table 5 prints the phosphorothioate-linked terminal T."
                    ).scenes,
                    statement=("The paper dA-tails the insert and Supplementary Table 5 prints a "
                               "terminal phosphorothioate-linked T on each ligating adapter arm.")))
    strand_polymerase = graph.input("strand_polymerase", "strand-displacing polymerase",
                                    MaterialKind.CATALYST)
    linked = graph.product("linked", "physically linked copies of both original complements",
                           MaterialKind.MOLECULAR, state=codec_linked_duplex())
    graph.react("link_duplex", "Copy and physically link both complements by strand displacement",
                StepKind.EXTENSION,
                ((adapted, InputRole.TEMPLATE), (strand_polymerase, InputRole.CATALYST)),
                (linked,))
    reviewed = ["duplex", "end_repair_enzyme", "end_repair", "codec_adapter", "ligase",
                "ligate", "strand_polymerase", "link_duplex"]
    finish_pcr(graph, "library", linked, final, primers, d.sequencing, reviewed,
                label="CODEC linked duplex")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
