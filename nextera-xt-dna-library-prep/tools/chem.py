"""Local source and reaction model for nextera-xt-dna-library-prep."""
from chemdraw import Construct, Scene, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from illumina_vendor_workflows import Definition, _duplex, _finish_pcr, _nextera_final, _tagged_intermediate, _workflow, nextera_final_rows
from types import SimpleNamespace


def _nextera_xt_workflow(final: Construct) -> Workflow:
    tagged = _tagged_intermediate("small-genome / amplicon insert")
    return _workflow("1 ng double-stranded DNA", _duplex("small-genome or amplicon DNA"), (
        ("Solution-transposome tagmentation", "tagged fragments", tagged,
         "The engineered transposome fragments DNA and transfers partial adapters in one reaction."),
        ("Neutralize transposase", "neutralized tagged fragments", tagged,
         "Neutralization stops further tagmentation before amplification."),
        ("Limited-cycle PCR with two index primers", "dual-index sequencing library",
         Scene.duplex(list(final), label="Nextera XT library"),
         "PCR both amplifies the insert and completes P5, P7, i5 and i7 sequences."),
    ))

def _nextera_xt_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "1-ng double-stranded small-genome or amplicon DNA",
                      MaterialKind.BIOLOGICAL,
                      state=_duplex("small-genome or amplicon DNA"), primary=True)
    transposome = graph.input("transposome", "solution transposome carrying partial Nextera adapters",
                              MaterialKind.ADAPTER, inferred=True,
                              source=source + " Commercial transposome composition is proprietary.")
    tn5 = graph.input("tn5", "engineered transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "heterologously tagged DNA fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("small-genome / amplicon insert"))
    unusable = graph.product("unusable_tn5_products", "same-end or non-library tagmentation products",
                             MaterialKind.WASTE)
    graph.react("tagment", "Fragment and tag DNA with the solution transposome",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, unusable))
    graph.discard(unusable, "cannot bridge the heterologous P5/P7 library arms")
    neutralizer = graph.input("neutralizer", "transposase-neutralization condition",
                              MaterialKind.CONDITION)
    neutral = graph.product("neutral", "neutralized tagged fragments",
                            MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("neutralize", "Neutralize transposase before PCR", StepKind.QUENCHING,
                ((tagged, InputRole.TEMPLATE), (neutralizer, InputRole.CONDITION)),
                (neutral,))
    reviewed = ["dna", "transposome", "tn5", "tagment", "neutralizer", "neutralize"]
    _finish_pcr(graph, "library", neutral, library, primers, d.sequencing_intro,
                reviewed, label="Nextera XT")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "Nextera XT DNA Library Prep", "01_nextera_xt.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nextera-xt/15031942_07_nextera-xt-product-documentation.pdf">Nextera XT DNA Library Prep Product Documentation 15031942 v07</a>.',
        "A solution transposome fragments and tags 1 ng DNA; neutralization and limited-cycle PCR complete a dual-index Nextera library.",
        "The commercial transposome composition is proprietary. Public Nextera mosaic-end, primer-entry and sequencing-primer sequences are imported from the central canonical model.",
        "A representative heterologous S5/S7 product; same-end tagmentation products cannot bridge P5 to P7 and are not the sequenced library.",
        "Standard Nextera primers read the insert from both mosaic-end boundaries and read i5/i7 sample indexes.",
        "nextera", "small-genome / amplicon insert", _nextera_xt_workflow)


def protocol_graph():
    return _nextera_xt_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)

FINAL_LIBRARY, SEQ_PRIMERS = _nextera_final(DEFINITION.insert_name, DEFINITION.title + " library")
d = DEFINITION
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
    CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
    workflow=lambda: _nextera_xt_workflow(FINAL_LIBRARY),
    final_rows=lambda: nextera_final_rows(FINAL_LIBRARY),
    protocol_graph=protocol_graph,
)
globals().update(vars(SPEC))
