"""Local source and reaction model for cap-seq."""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Workflow
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import role, role_library
from rare_cell_and_affinity import _scene, _workflow
from types import SimpleNamespace

CAP_HUMAN_WASH_MM = 600
CAP_MOUSE_WASH_MM = 560

def cap_library() -> Construct:
    return role_library("CAP-seq library", (role("unmethylated CpG-rich insert"),))

def cap_workflow() -> Workflow:
    genomic = _scene("input genomic DNA prepared as for MAP-seq", role("genomic DNA fragment with CpG sites"))
    adapted = _scene("adapter-ligated genomic library", role("Solexa left adapter"),
                     role("genomic insert with CpG sites"), role("Solexa right adapter"))
    bound = _scene("CXXC-bound library", role("nickel-Sepharose-bound recombinant CXXC"),
                   role("adapter-ligated unmethylated CpG-rich fragment"))
    final = cap_library()
    return _workflow(genomic, (
        ("Ligate Solexa sequencing adaptors before purification, as for MAP [57]", adapted,
         "Adaptor sequences and end geometry are delegated to the predecessor and not assigned here."),
        ("Bind adapter-ligated input to recombinant CXXC on nickel Sepharose", bound,
         "CXXC affinity selects clustered unmethylated CpGs; the adapters were installed before purification."),
        (f"Wash with {CAP_HUMAN_WASH_MM} mM NaCl (human) or "
         f"{CAP_MOUSE_WASH_MM} mM NaCl (mouse)", bound,
         "The lower-stringency mouse wash retains its relatively CpG-poor islands."),
        ("Elute at 1 M NaCl; pool, concentrate and precipitate as for MAP [57]; sequence",
         _scene("CAP-seq library", *tuple(final)),
         "The defining experiment used single-end Solexa sequencing."),
    ))

SPEC = SimpleNamespace(
    TITLE="CAP-seq — CXXC affinity purification sequencing",
    NOTES="01_cap_seq.html",
    SOURCE='Defining source: <a href="https://doi.org/10.1371/journal.pgen.1001134">Illingworth et al., PLOS Genetics (2010)</a>.',
    SUMMARY="Adapter-ligated genomic DNA is selected on immobilized recombinant CXXC for clusters of unmethylated CpGs, then eluted and sequenced.",
    CAVEAT="The defining paper delegates Solexa library preparation to Skene et al. (2010), whose Supplemental Data is not in the local source bundle. Adapter and sequencing-primer bases have not been verified, so the library ends are shown by role.",
    FINAL_LIBRARY=cap_library(),
    SEQ_PRIMERS=(),
    READ_LENGTHS=None,
    FINAL_CAPTION="Solexa library enriched for fragments carrying clustered unmethylated CpGs.",
    SEQUENCING_INTRO="",
    SEQUENCING_ENDING="The defining study used single-end Solexa reads (Wellcome Trust Sanger Institute) mapped with MAQ; reads with mapping score of at least 30 were retained computationally. Read length is not stated in the inspected Methods.",
    SEQUENCING_UNAVAILABLE="Exact sequencing-primer sites remain unresolved pending the predecessor's Supplemental Data and verification of its Solexa adapters.",
    workflow=cap_workflow,
    oligos=lambda: (),
)

def _cap_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    genomic = graph.input("genomic_dna", "genomic DNA from human or mouse sperm, blood or cerebellum",
                          MaterialKind.BIOLOGICAL,
                          state=_scene("genomic DNA", role("genomic DNA with CpG sites")), primary=True)
    map_prep = graph.input("map_prep", "DNA preparation as for MAP-seq (reference 57)",
                           MaterialKind.CONDITION)
    prepared = graph.product("prepared_dna", "CAP input DNA prepared as for MAP-seq",
                             MaterialKind.MOLECULAR,
                             state=_scene("prepared input DNA", role("genomic DNA fragment; preparation delegated")))
    graph.react("prepare_dna", "Prepare input DNA as for MAP-seq", StepKind.MAGIC,
                ((genomic, InputRole.TEMPLATE), (map_prep, InputRole.CONDITION)), (prepared,),
                reason="Methods line 218 says DNA was prepared as for MAP-seq and delegates the recipe to reference 57 (Skene et al. 2010), whose Supplemental Data was not obtained; fragmentation method and fragment ends are not assigned.")
    adaptors = graph.input("adaptors", "Solexa sequencing adaptors (bases not printed in the defining article)",
                           MaterialKind.ADAPTER)
    adapted = graph.product("adapted_dna", "adapter-ligated genomic DNA library",
                            MaterialKind.MOLECULAR,
                            state=_scene("adapter-ligated genomic library",
                                         role("Solexa left adapter"),
                                         role("genomic insert with CpG sites"),
                                         role("Solexa right adapter")))
    graph.react("ligate", "Ligate Solexa sequencing adaptors before CXXC purification",
                StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adaptors, InputRole.STRUCTURAL)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Methods line 226 delegates adaptor ligation to the MAP protocol of reference 57; adaptor bases, ends and ligase are not given in the defining article."))
    cxxc = graph.input("cxxc", "DNaseI-treated recombinant CXXC bound to nickel-charged Sepharose (50 mg/ml bead volume)",
                       MaterialKind.SELECTION)
    binding = graph.input("binding_condition", "column buffer containing 0.1 M NaCl",
                          MaterialKind.CONDITION)
    bound = graph.product("bound", "CXXC-bound unmethylated CpG-rich fragments",
                          MaterialKind.MOLECULAR, state=adapted.state)
    depleted = graph.product("depleted", "CXXC-unbound genomic library", MaterialKind.WASTE)
    graph.react("bind", "Bind adapter-ligated input to immobilized CXXC in 0.1 M NaCl",
                StepKind.AFFINITY_CAPTURE,
                ((adapted, InputRole.TEMPLATE), (cxxc, InputRole.SELECTION),
                 (binding, InputRole.CONDITION)),
                (bound, depleted))
    graph.discard(depleted, "lacks sufficient clustered unmethylated CpG affinity")
    wash = graph.input("wash_condition",
                       f"{CAP_HUMAN_WASH_MM} mM NaCl (human) or {CAP_MOUSE_WASH_MM} mM NaCl (optimized mouse) wash",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed CXXC-bound CpG-rich fragments",
                           MaterialKind.MOLECULAR, state=bound.state)
    wash_waste = graph.product("wash_waste", "fragments removed by high-salt washing",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash the CXXC-bound library", StepKind.CLEANUP,
                ((bound, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste))
    graph.discard(wash_waste, "does not remain bound at the chosen species-specific stringency")
    elution = graph.input("elution_condition", "1-M NaCl elution buffer",
                          MaterialKind.CONDITION)
    eluted = graph.product("eluted", "pooled 1-M NaCl eluate fractions",
                           MaterialKind.MOLECULAR, state=washed.state)
    graph.react("elute", "Elute with 1 M NaCl and pool eluted fractions", StepKind.EXTRACTION,
                ((washed, InputRole.TEMPLATE), (elution, InputRole.CONDITION)), (eluted,))
    precipitation = graph.input("precipitation", "concentration and precipitation as for MAP (reference 57)",
                                MaterialKind.CONDITION)
    library = graph.product("library", "CXXC-enriched CAP-seq library",
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    supernatant = graph.product("precipitation_supernatant", "precipitation supernatant and salt",
                                MaterialKind.WASTE)
    graph.react("precipitate", "Concentrate and precipitate pooled eluate as for MAP",
                StepKind.CLEANUP,
                ((eluted, InputRole.TEMPLATE), (precipitation, InputRole.CONDITION)),
                (library, supernatant))
    graph.discard(supernatant, "is separated from the precipitated DNA")
    graph.sequence(library, SequencingPlan(limitation=spec.SEQUENCING_ENDING))
    graph.annotate(Highlight("precipitate", "Post-enrichment library steps are delegated",
                             "Methods line 226 says CAP Solexa libraries were prepared as previously described for MAP (reference 57). Whether any amplification or size selection followed CXXC elution is not stated in the defining article and is not drawn."))
    reviewed = ["genomic_dna", "map_prep", "prepare_dna", "adaptors", "ligate", "cxxc",
                "binding_condition", "bind", "wash_condition", "wash", "elution_condition",
                "elute", "precipitation", "precipitate", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _cap_graph(SPEC)

SPEC.protocol_graph = protocol_graph
