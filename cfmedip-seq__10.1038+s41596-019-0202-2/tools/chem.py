"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("cfMeDIP-seq — cell-free DNA methylome enrichment", "10.1038/s41596-019-0202-2", "Shen et al., Nature Protocols (2019)", "Adapter-ligated low-input cfDNA is supplemented with controls and filler, denatured, and enriched with a 5mC antibody without bisulfite conversion.", "truseq", ("methylated cfDNA fragment",), "Ordinary paired-end Illumina reads report the endpoints and sequence of immunoprecipitated cfDNA fragments.", "DNA")


from protocol_specs import workflow


def _source_workflow():
    return workflow(molecule("plasma cfDNA", ("naturally fragmented cfDNA",)), (
        ("End-repair and add one 3-prime dA to each end",
         molecule("dA-tailed cfDNA", ("3-prime dA", "cfDNA", "3-prime dA")),
         "The cohesive dA ends are the substrates for the following adapter ligation."),
        ("Ligate Illumina adapters to the dA-tailed cfDNA",
         molecule("adapter-ligated cfDNA", ("P5-side adapter", "cfDNA", "P7-side adapter"), joined=True,
                  junctions=(("P5-side adapter", "cfDNA", "adapter ligation"),
                             ("cfDNA", "P7-side adapter", "adapter ligation"))),
         "Adapters are installed before methylated-fragment enrichment."),
        ("Add methylated and unmethylated spike-ins plus λ filler DNA",
         molecule("low-input library mixture", ("cfDNA library", "controls", "filler DNA")),
         "Filler supports low-input immunoprecipitation; controls measure recovery."),
        ("Denature and immunoprecipitate 5-methylcytosine-containing fragments",
         molecule("5mC-enriched library", ("adapter", "methylated cfDNA fragment", "adapter"), joined=True),
         "The assay enriches methylated fragments without bisulfite conversion."),
    ))

SPEC = get('cfmedip-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("cfdna", "naturally fragmented plasma cell-free DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("plasma cfDNA", ("naturally fragmented cfDNA",)), primary=True)
    repair = graph.input("repair_mix", "cfDNA end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired cfDNA", MaterialKind.MOLECULAR,
                             state=dna.state)
    graph.react("end_repair", "End-repair cfDNA fragments", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "3-prime-dA-tailed cfDNA", MaterialKind.MOLECULAR,
                           state=molecule("dA-tailed cfDNA", ("3-prime dA", "cfDNA", "3-prime dA")))
    graph.react("tail", "Add one 3-prime dA to each fragment end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("adapters", "Illumina dT-overhang adapters", MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated cfDNA", MaterialKind.MOLECULAR,
                            state=molecule("adapter-ligated cfDNA", (
                                "P5-side adapter", "cfDNA", "P7-side adapter"),
                                junctions=(("P5-side adapter", "cfDNA", "adapter ligation"),
                                           ("cfDNA", "P7-side adapter", "adapter ligation"))))
    graph.react("ligate", "Ligate Illumina adapters to dA-tailed cfDNA", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The role-level source summary specifies dA/dT ligation without fixing one vendor adapter revision."))
    controls = graph.input("controls", "methylated/unmethylated controls and lambda filler DNA",
                           MaterialKind.MOLECULAR)
    mixture = graph.product("mixture", "low-input cfDNA library mixture",
                            MaterialKind.MOLECULAR,
                            state=molecule("library mixture", ("cfDNA library", "controls", "filler DNA")))
    graph.react("pool_controls", "Add controls and filler DNA",
                StepKind.POOL,
                ((adapted, InputRole.TEMPLATE), (controls, InputRole.TEMPLATE)), (mixture,))
    denaturing = graph.input("denaturing_condition", "single-strand denaturing condition",
                             MaterialKind.CONDITION)
    single = graph.product("single_strand", "denatured cfDNA library strands",
                           MaterialKind.MOLECULAR, state=mixture.state)
    complement = graph.product("complement", "opposite denatured strands", MaterialKind.WASTE)
    graph.react("denature", "Denature the low-input library mixture",
                StepKind.STRAND_SEPARATION,
                ((mixture, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)),
                (single, complement))
    graph.discard(complement, "represented by the retained single-stranded library population")
    antibody = graph.input("antibody", "5-methylcytosine immunoprecipitation antibody",
                           MaterialKind.SELECTION)
    enriched = graph.product("enriched", "5mC-enriched adapter-ligated cfDNA",
                             MaterialKind.MOLECULAR, state=single.state)
    depleted = graph.product("depleted", "unmethylated antibody-depleted material",
                             MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate 5mC-containing fragments",
                StepKind.AFFINITY_CAPTURE,
                ((single, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, depleted))
    graph.discard(depleted, "does not belong to the antibody-enriched methylated fraction")
    reviewed = ["cfdna", "repair_mix", "end_repair", "tailing_enzyme", "tail", "adapters",
                "ligase", "ligate", "controls", "pool_controls", "denaturing_condition",
                "denature", "antibody", "immunoprecipitate"]
    finish_pcr(graph, "library", enriched, final, primers, d.sequencing, reviewed,
                label="cfMeDIP-seq")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
