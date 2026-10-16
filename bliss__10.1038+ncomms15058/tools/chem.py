"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from bliss_chem import adapter_a1 as bliss_adapter_a1
from chemdraw import Segment, oligo
from illumina import TRUSEQ_SMALL_RNA_RA3
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_rt_pcr
from rna_library_operations import molecule

DEFINITION = Definition("BLISS — breaks labeling in situ and sequencing", "10.1038/ncomms15058", "Yan et al., Nature Communications (2017)", "Native DSBs are blunted and ligated in situ to a T7/RA5 adapter carrying a sample barcode and UMI, then linearly amplified by IVT.", "custom", ("UMI", "sample barcode", "break-adjacent DNA"), "The genomic read begins immediately beside the captured DSB; the inline UMI distinguishes molecular break events.", "DNA")


from protocol_specs import workflow
from bliss_chem import oligos as bliss_oligos


def _oligos():
    return (*bliss_oligos(), oligo("RA3 RNA 3′ adapter (custom-synthesized TruSeq Small RNA sequence)",
                                  [Segment("RA3", TRUSEQ_SMALL_RNA_RA3)]))


def _source_workflow():
    return workflow(molecule("fixed DNA", ("native double-strand break",)), (
        ("Blunt the break ends in situ",
         molecule("blunted break", ("genomic DNA", "blunt DSB end")),
         "Only ends already present in fixed material are exposed for tagging."),
        ("Ligate a UMI, sample barcode, RA5 and T7-promoter adapter",
         molecule("tagged break", ("genomic break", "UMI", "sample barcode", "RA5", "T7 promoter"), joined=True,
                  junctions=(("genomic break", "UMI", "break-adapter ligation"),)),
         "The ligation boundary marks the nucleotide adjacent to the break."),
        ("Linearly amplify by T7 transcription",
         molecule("BLISS IVT RNA", ("RA5", "UMI", "sample barcode", "break-adjacent RNA"), joined=True),
         "The break-adjacent sequence is copied into RNA."),
        ("Ligate the RA3 adapter to the IVT RNA, reverse-transcribe and index-amplify",
         molecule("BLISS library", ("RA5", "UMI", "sample barcode", "break-adjacent DNA"), joined=True),
         "The original Protocol Exchange method names RA3, RTP, RP1 and a chosen RPI primer."),
    ))

SPEC = get('bliss', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "fixed DNA containing native double-strand breaks",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("fixed DNA", ("native double-strand break",)), primary=True)
    repair = graph.input("blunting_enzyme", "break-end blunting enzyme mix",
                         MaterialKind.CATALYST)
    blunted = graph.product("blunted", "blunted native DNA break ends", MaterialKind.MOLECULAR,
                            state=molecule("blunted break", ("genomic DNA", "blunt DSB end")))
    graph.react("blunt", "Blunt the native break ends in situ", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (blunted,))
    adapter = graph.input("break_adapter", "published UMI/sample-barcode/RA5/T7 break adapter",
                          MaterialKind.ADAPTER, state=bliss_adapter_a1(),
                          source=source + "; Supplementary Data 1 and Supplementary Fig. 1a")
    ligase = graph.input("ligase", "break-adapter DNA ligase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "adapter-tagged break ends", MaterialKind.MOLECULAR,
                           state=molecule("tagged break", (
                               "genomic break", "UMI", "sample barcode", "RA5", "T7 promoter"),
                               junctions=(("genomic break", "UMI", "break-adapter ligation"),)))
    graph.react("ligate", "Ligate the barcoded T7 adapter to break ends", StepKind.LIGATION,
                ((blunted, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.known(
                    "blunt", bliss_adapter_a1(),
                    statement=("Supplementary Fig. 1a marks ligation at the blunt, "
                               "5-prime-phosphorylated sample-barcode end; Supplementary "
                               "Data 1 prints both oligos.")))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rna = graph.product("rna", "linearly amplified break-adjacent RNA", MaterialKind.MOLECULAR,
                        state=tagged.state)
    graph.react("transcribe", "Linearly amplify tagged breaks by T7 transcription",
                StepKind.TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    ra3 = graph.input("ra3_adapter", "custom TruSeq Small RNA RA3 3′ adapter",
                      MaterialKind.ADAPTER, source=source + "; Protocol Exchange v2 pp. 3, 10")
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 2, truncated (NEB M0242L)", MaterialKind.CATALYST)
    ligated_rna = graph.product("ra3_ligated_rna", "RA3-ligated BLISS IVT RNA",
                               MaterialKind.MOLECULAR,
                               state=molecule("RA3-ligated IVT RNA",
                                              ("RA5", "UMI", "sample barcode", "break-adjacent RNA", "RA3"),
                                              joined=True))
    graph.react("ligate_ra3", "Ligate RA3 to the IVT RNA 3′ end", StepKind.LIGATION,
                ((rna, InputRole.TEMPLATE), (ra3, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (ligated_rna,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The custom RA3 bases are printed by Illumina, but the original BLISS protocol does not specify its custom-synthesis end modifications."))
    reviewed = ["dna", "blunting_enzyme", "blunt", "break_adapter", "ligase", "ligate",
                "t7_polymerase", "transcribe", "ra3_adapter", "rna_ligase", "ligate_ra3"]
    finish_rt_pcr(graph, "library", ligated_rna, final, primers, d.sequencing, reviewed,
                   label="BLISS break")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
