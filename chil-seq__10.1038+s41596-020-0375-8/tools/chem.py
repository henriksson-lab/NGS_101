"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from chil_chem import AD1_NOMX as CHIL_AD1_NOMX
from chil_chem import AD2_INDEX as CHIL_AD2_INDEX
from pcr import PCRGeometry
from pcr import PCRPrimer
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from rna_library_operations import molecule

DEFINITION = Definition("ChIL-seq — chromatin integration labeling", "10.1038/s41596-020-0375-8", "Handa et al., Nature Protocols (2020)", "An oligo-conjugated antibody targets chromatin; nearby Tn5 integration and T7 transcription linearly amplify target-proximal DNA before sequencing.", "custom", ("ChIL-probe end", "target-proximal DNA"), "Illumina sequencing starts from the ChIL-probe-derived end and reports nearby genomic sequence.", "DNA")


from protocol_specs import workflow, _payload
from chil_chem import oligos as chil_oligos
from batch_ngs import nextera_library


def _source_workflow():
    return workflow(molecule("fixed chromatin", ("antibody target", "nearby DNA")), (
        ("Bind an antibody carrying a ChIL oligonucleotide",
         molecule("targeted chromatin", ("antibody–ChIL probe", "nearby genomic DNA")),
         "The probe contains a T7 promoter and sequence used to recover nearby DNA."),
        ("Tagment nearby genomic DNA in situ",
         molecule("tagged target DNA", ("ChIL probe", "target-proximal genomic DNA", "Tn5 end")),
         "The antibody position and the transposed DNA are joined during subsequent copying."),
        ("Extend, transcribe linearly with T7, then make a cDNA library",
         molecule("ChIL-seq cDNA", ("ChIL-probe end", "target-proximal DNA copy"), joined=True),
         "T7 transcription amplifies the selected genomic neighbourhood before PCR."),
    ))

SPEC = get('chil-seq', definition=DEFINITION,
           workflow_fn=_source_workflow,
           final_override=nextera_library(_payload(DEFINITION), DEFINITION.title + " final library", dual_index=False),
           oligos_fn=chil_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    chromatin = graph.input("chromatin", "fixed chromatin containing an antibody target",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("fixed chromatin", ("antibody target", "nearby DNA")),
                            primary=True)
    antibody = graph.input("antibody", "antibody carrying the ChIL oligonucleotide",
                           MaterialKind.SELECTION)
    bound = graph.product("bound", "ChIL-probe-targeted chromatin", MaterialKind.MOLECULAR,
                          state=molecule("targeted chromatin", (
                              "antibody–ChIL probe", "nearby genomic DNA")))
    unbound = graph.product("unbound", "untargeted chromatin", MaterialKind.WASTE)
    graph.react("bind", "Bind the ChIL-oligo-conjugated antibody",
                StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (bound, unbound))
    graph.discard(unbound, "does not retain the selected chromatin target")
    transposome = graph.input("transposome", "loaded ChIL-compatible transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented target-proximal genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagged target DNA", (
                               "ChIL probe", "target-proximal genomic DNA", "Tn5 end")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment nearby genomic DNA in situ", StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain target-proximal tagged DNA")
    extension_primer = graph.input("extension_primer", "published ChIL-probe extension architecture",
                                   MaterialKind.OLIGO,
                                   source=source + "; Supplementary Table 4")
    polymerase = graph.input("extension_polymerase", "DNA polymerase", MaterialKind.CATALYST)
    extended = graph.product("extended", "ChIL probe joined to a nearby DNA copy",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("extend", "Extend from the ChIL probe into target-proximal DNA",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (extension_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (extended,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rna = graph.product("rna", "linearly amplified target-proximal RNA copy",
                        MaterialKind.MOLECULAR, state=extended.state)
    graph.react("transcribe", "Linearly amplify with T7 transcription",
                StepKind.TRANSCRIPTION,
                ((extended, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,))
    rt_primer = graph.input("library_rt_primer", "Read2+N11 first-strand synthesis primer",
                            MaterialKind.OLIGO,
                            source=source + "; Supplementary Table 4")
    rt = graph.input("library_rt_enzyme", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "ChIL target-proximal cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("ChIL cDNA", ("target-proximal cDNA",)))
    graph.react("library_rt", "Reverse-transcribe ChIL RNA with Read2+N11",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    forward = graph.input("library_forward", "Ad1_noMX PCR primer", MaterialKind.OLIGO,
                          source=source + "; Supplementary Table 4")
    reverse = graph.input("library_reverse", "Ad2_index PCR primer", MaterialKind.OLIGO,
                          source=source + "; Supplementary Table 4")
    pcr_polymerase = graph.input("library_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("library_library", "ChIL target-proximal sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify and index the ChIL-seq library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, PCRPrimer("Ad1_noMX", CHIL_AD1_NOMX),
                    PCRPrimer("Ad2_index", CHIL_AD2_INDEX),
                    statement="Supplementary Table 4 prints both complete PCR primers."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed = ["chromatin", "antibody", "bind", "transposome", "tn5", "tagment",
                "extension_primer", "extension_polymerase", "extend", "t7_polymerase",
                "transcribe", "library_rt_primer", "library_rt_enzyme", "library_rt",
                "library_forward", "library_reverse", "library_polymerase",
                "library_pcr", "library_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
