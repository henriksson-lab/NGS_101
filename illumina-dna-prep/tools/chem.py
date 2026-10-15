"""Local source and reaction model for illumina-dna-prep."""
from chemdraw import Construct, Scene, Segment, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from illumina_vendor_workflows import Definition, _duplex, _finish_pcr, _nextera_final, _single, _tagged_intermediate, _workflow, nextera_final_rows
from types import SimpleNamespace


def _dna_prep_workflow(final: Construct) -> Workflow:
    tagged = _tagged_intermediate("genomic DNA insert")
    bead = _single("bead-bound", [Segment("bead-linked transposome", "[BLT]",
                                                   placeholder=True),
                                           Segment("genomic DNA", "N" * 36,
                                                   placeholder=True)],
                   label="BLT-bound genomic DNA")
    return _workflow("genomic DNA", _duplex("genomic DNA"), (
        ("Bind genomic DNA to bead-linked transposomes", "BLT-bound DNA", bead,
         "At inputs of 100 ng or more, the finite bead-linked complexes support yield normalization."),
        ("On-bead tagmentation", "bead-bound tagged fragments", tagged,
         "BLT fragments a set number of DNA molecules and transfers partial adapters."),
        ("Post-tagmentation wash", "washed tagged fragments", tagged,
         "The solid phase allows removal of unbound DNA and reaction components."),
        ("Input-dependent limited-cycle index PCR", "dual-index sequencing library",
         Scene.duplex(list(final), label="Illumina DNA Prep library"),
         "PCR completes the platform arms; lower input requires additional cycles."),
    ))

def _dna_prep_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                      state=_duplex("genomic DNA"), primary=True)
    beads = graph.input("blt_beads", "bead-linked transposome complexes",
                        MaterialKind.SELECTION, inferred=True,
                        source=source + " Reagent-level BLT composition is proprietary.")
    bound = graph.product("bound", "genomic DNA bound to bead-linked transposomes",
                          MaterialKind.MOLECULAR,
                          state=_single("bead-bound", [
                              Segment("bead-linked transposome", "[BLT]", placeholder=True),
                              Segment("genomic DNA", "N" * 36, placeholder=True)],
                              label="BLT-bound genomic DNA"))
    unbound = graph.product("unbound", "unbound excess genomic DNA", MaterialKind.WASTE)
    graph.react("bind_blt", "Bind genomic DNA to bead-linked transposomes",
                StepKind.BINDING,
                ((dna, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (bound, unbound))
    graph.discard(unbound, "is removed by the finite bead-linked transposome capacity")
    adapter = graph.input("transposome_adapter", "BLT partial Nextera adapter",
                          MaterialKind.ADAPTER, inferred=True,
                          source=source + " Public endpoint architecture is used; BLT reagent composition is proprietary.")
    tn5 = graph.input("tn5", "bead-linked transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "bead-bound tagged genomic fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("genomic DNA insert"))
    tn5_waste = graph.product("tn5_waste", "non-library or same-end tagmentation products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment genomic DNA on bead", StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not form a heterologous amplifiable library fragment")
    wash = graph.input("wash_condition", "post-tagmentation wash condition",
                       MaterialKind.CONDITION)
    washed = graph.product("washed", "washed bead-bound tagged fragments",
                           MaterialKind.MOLECULAR, state=tagged.state)
    wash_waste = graph.product("wash_waste", "removed unbound DNA and reaction components",
                               MaterialKind.WASTE)
    graph.react("wash", "Wash bead-bound tagged fragments", StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (wash, InputRole.CONDITION)),
                (washed, wash_waste))
    graph.discard(wash_waste, "is removed during the post-tagmentation wash")
    reviewed = ["dna", "blt_beads", "bind_blt", "transposome_adapter", "tn5",
                "tagment", "wash_condition", "wash"]
    _finish_pcr(graph, "library", washed, library, primers, d.sequencing_intro,
                reviewed, label="Illumina DNA Prep")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "Illumina DNA Prep", "01_illumina_dna_prep.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/illumina_prep/1000000025416_12_idp.pdf">Illumina DNA Prep Reference Guide 1000000025416 v12</a>.',
        "Bead-linked transposomes fragment and tag genomic DNA, a wash removes unbound material, and input-dependent limited-cycle PCR completes unique-dual-index libraries.",
        "Illumina does not disclose the reagent-level bead/transposome composition. The molecular endpoint uses the public Nextera adapter and primer architecture.",
        "A representative Illumina DNA Prep fragment after bead-linked tagmentation, wash and dual-index PCR.",
        "Standard Nextera primers read both genomic-insert ends and the unique dual sample indexes.",
        "nextera", "genomic DNA insert", _dna_prep_workflow)


def protocol_graph():
    return _dna_prep_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)

FINAL_LIBRARY, SEQ_PRIMERS = _nextera_final(DEFINITION.insert_name, DEFINITION.title + " library")
d = DEFINITION
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
    CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
    workflow=lambda: _dna_prep_workflow(FINAL_LIBRARY),
    final_rows=lambda: nextera_final_rows(FINAL_LIBRARY),
    protocol_graph=protocol_graph,
)
globals().update(vars(SPEC))
