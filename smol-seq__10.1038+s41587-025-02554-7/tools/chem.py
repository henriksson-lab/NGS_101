from types import SimpleNamespace
from chemdraw import feature
from protocol_model import chain, role_library

from readout_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, _authority, _unavailable, molecule, role
DEFINITION = Definition("smol-seq — small-molecule sequencing",
        "10.1038/s41587-025-02554-7", "Tan and Fraser, Nature Biotechnology (2026)",
        "A structure-switching aptamer releases an analyte-specific DNA barcode when it binds its small-molecule target; barcode sequencing converts metabolite abundance into counts.",
        "The sequenced feature barcode identifies the aptamer—and therefore the analyte—not the small molecule itself.",
        "01_smol_seq.html")


def _smol_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    sensor = graph.input(
        "aptamer_sensor", "structure-switching aptamer sensor",
        MaterialKind.MOLECULAR,
        state=molecule("structure-switching aptamer sensor", (
            role("target-binding aptamer"), ANALYTE_BARCODE, role("quenched capture strand"))),
        primary=True)
    analyte = graph.input("small_molecule_analyte", "cognate small-molecule analyte",
                          MaterialKind.MOLECULAR)
    bound = graph.product(
        "bound_sensor", "analyte-bound switched aptamer", MaterialKind.MOLECULAR,
        state=molecule("analyte-bound aptamer", (role("small-molecule analyte"),
            role("folded aptamer"), ANALYTE_BARCODE)))
    unbound = graph.product("unbound_sensor", "unbound aptamer sensor", MaterialKind.WASTE)
    graph.react(
        "analyte_binding", "Bind cognate analyte and switch aptamer structure",
        StepKind.BINDING,
        ((sensor, InputRole.TEMPLATE), (analyte, InputRole.ANALYTE)), (bound, unbound))
    graph.discard(unbound, "does not release an analyte-responsive reporter")

    release = graph.input("release_condition", "reporter-release condition",
                          MaterialKind.CONDITION)
    reporter = graph.product(
        "released_reporter", "released analyte-identity reporter", MaterialKind.MOLECULAR,
        state=molecule("released reporter", (ANALYTE_BARCODE,)))
    spent = graph.product("spent_sensor", "spent aptamer/capture-strand complex",
                          MaterialKind.WASTE)
    graph.react(
        "reporter_release", "Release the analyte-specific DNA reporter",
        StepKind.STRAND_SEPARATION,
        ((bound, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (reporter, spent))
    graph.discard(spent, "the released DNA identity reporter is the sequenced analyte proxy")

    forward = graph.input("index_forward", "reporter index-PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("index_reverse", "reporter index-PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "reporter PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "indexed smol-seq reporter library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(
        "index_pcr", "PCR-index released reporter molecules", StepKind.PCR,
        ((reporter, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
         (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
        (library,), geometry=_unavailable(
            final, "The aptamer panel determines the reporter family; no one reporter sequence is universal."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "aptamer_sensor", "small_molecule_analyte", "analyte_binding",
        "release_condition", "reporter_release", "index_forward", "index_reverse",
        "pcr_polymerase", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'smol-seq'
    return _smol_graph(d, final)

ANALYTE_BARCODE = role("aptamer identity barcode", tag="cbc",
                       molecular_feature=feature("analyte_identity", "feature_barcode", "unknown"))
SAMPLE_INDEX = role("sample index", tag="cbc", length_bp=8,
                    molecular_feature=feature("sample_index", "sample_index", "unknown"))

def _smol():
    return chain(molecule("structure-switching aptamer sensor", (
        role("target-binding aptamer"), ANALYTE_BARCODE, role("quenched capture strand"))), (
        ("Bind the cognate small molecule and switch aptamer structure",
         molecule("analyte-bound aptamer", (role("small-molecule analyte"),
                                             role("folded aptamer"), ANALYTE_BARCODE)),
         "Target binding destabilizes the duplex that retained the DNA identity barcode."),
        ("Release the analyte-specific DNA barcode",
         molecule("released reporter", (ANALYTE_BARCODE,)),
         "Reporter abundance is a proxy for the amount of the metabolite recognized by that aptamer."),
        ("PCR-index the released reporters and sequence them",
         molecule("smol-seq library", (role("P5 end"), SAMPLE_INDEX,
                                        ANALYTE_BARCODE, role("P7 end"))),
         "Multiplexed barcode counts convert chemically diverse analytes into one DNA-sequencing readout."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library("smol-seq final library", (ANALYTE_BARCODE,), index=SAMPLE_INDEX)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("The defining source establishes adapter and identifier roles, but the page does not promote panel-specific or supplementary oligo sequences into universal protocol constants."),
    workflow=_smol, oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
