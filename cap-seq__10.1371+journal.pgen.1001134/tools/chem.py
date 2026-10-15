"""Local source and reaction model for cap-seq."""
from chemdraw import Construct, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import role, role_library
from rare_cell_and_affinity import _scene, _workflow
from types import SimpleNamespace

CAP_HUMAN_WASH_MM = 600
CAP_MOUSE_WASH_MM = 560

def cap_library() -> Construct:
    return role_library("CAP-seq library", (role("unmethylated CpG-rich insert"),))

def cap_workflow() -> Workflow:
    adapted = _scene("adapter-ligated genomic library", role("Solexa left adapter"),
                     role("genomic insert with CpG sites"), role("Solexa right adapter"))
    bound = _scene("CXXC-bound library", role("nickel-Sepharose-bound recombinant CXXC"),
                   role("adapter-ligated unmethylated CpG-rich fragment"))
    final = cap_library()
    return _workflow(adapted, (
        ("Bind adapter-ligated input to recombinant CXXC on nickel Sepharose", bound,
         "CXXC affinity selects clustered unmethylated CpGs; the adapters were installed before purification."),
        (f"Wash with {CAP_HUMAN_WASH_MM} mM NaCl (human) or "
         f"{CAP_MOUSE_WASH_MM} mM NaCl (mouse)", bound,
         "The lower-stringency mouse wash retains its relatively CpG-poor islands."),
        ("Elute at 1 M NaCl, pool and sequence",
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
    SEQUENCING_ENDING="The defining study used single-end Solexa reads to map CXXC-enriched fragments.",
    SEQUENCING_UNAVAILABLE="Exact sequencing-primer sites remain unresolved pending the predecessor's Supplemental Data and verification of its Solexa adapters.",
    workflow=cap_workflow,
    oligos=lambda: (),
)

def _cap_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    adapted = graph.input("adapted_dna", "adapter-ligated genomic DNA library",
                          MaterialKind.LIBRARY,
                          state=_scene("adapter-ligated genomic library",
                                       role("Solexa left adapter"),
                                       role("genomic insert with CpG sites"),
                                       role("Solexa right adapter")), primary=True)
    cxxc = graph.input("cxxc", "immobilized recombinant CXXC on nickel Sepharose",
                       MaterialKind.SELECTION)
    bound = graph.product("bound", "CXXC-bound unmethylated CpG-rich fragments",
                          MaterialKind.MOLECULAR, state=adapted.state)
    depleted = graph.product("depleted", "CXXC-unbound genomic library", MaterialKind.WASTE)
    graph.react("bind", "Bind adapter-ligated input to immobilized CXXC",
                StepKind.AFFINITY_CAPTURE,
                ((adapted, InputRole.TEMPLATE), (cxxc, InputRole.SELECTION)),
                (bound, depleted))
    graph.discard(depleted, "lacks sufficient clustered unmethylated CpG affinity")
    wash = graph.input("wash_condition", "species-adjusted high-salt wash condition",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed CXXC-bound CpG-rich fragments",
                           MaterialKind.MOLECULAR, state=bound.state)
    wash_waste = graph.product("wash_waste", "fragments removed by high-salt washing",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash the CXXC-bound library", StepKind.CLEANUP,
                ((bound, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste))
    graph.discard(wash_waste, "does not remain bound at the chosen species-specific stringency")
    elution = graph.input("elution_condition", "1-M NaCl elution condition",
                          MaterialKind.CONDITION)
    library = graph.product("library", "CXXC-enriched CAP-seq library",
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    graph.react("elute", "Elute and pool CXXC-enriched fragments", StepKind.EXTRACTION,
                ((washed, InputRole.TEMPLATE), (elution, InputRole.CONDITION)), (library,))
    graph.sequence(library, SequencingPlan(limitation=spec.SEQUENCING_ENDING))
    reviewed = ["adapted_dna", "cxxc", "bind", "wash_condition", "wash",
                "elution_condition", "elute", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _cap_graph(SPEC)

SPEC.protocol_graph = protocol_graph
