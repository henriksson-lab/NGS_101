from types import SimpleNamespace
from chemdraw import feature
from protocol_model import chain, role_library

from readout_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, _authority, _unavailable, molecule, role
DEFINITION = Definition("TIRTL-seq — throughput-intensive rapid TCR library sequencing",
        "10.1038/s41592-025-02907-9", "Pogorelyy et al., Nature Methods (2026)",
        "Miniaturized plate RT, multiplex V/C amplification, plate barcodes and well-specific UDIs generate deep alpha/beta TCR repertoires whose chain pairs are inferred across replicate wells.",
        "TCR insert reads are associated with the PCR-I plate barcode and PCR-II well UDI pair.",
        "01_tirtl_seq.html")


def _tirtl_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "replicate_well_cells", "T cells distributed across replicate wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("cells distributed across replicate wells", (role("T-cell RNA"),)),
        primary=True)
    lysis = graph.input("lysis_condition", "miniaturized-well lysis condition",
                        MaterialKind.CONDITION)
    rt_primers = graph.input("rt_primer_pool", "TCR reverse-transcription primer pool",
                             MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(
        "tcr_cdna", "TCR alpha/beta cDNA",
        MaterialKind.MOLECULAR,
        state=molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))))
    graph.react(
        "reverse_transcription", "Lyse and reverse-transcribe each replicate well",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cells, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
         (rt_primers, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)), (cdna,))

    vc_forward = graph.input("vc_forward_pool", "Nextera-tailed V-segment primer pool",
                             MaterialKind.OLIGO)
    vc_reverse = graph.input("vc_reverse_pool", "plate-barcoded C-segment primer pool",
                             MaterialKind.OLIGO)
    pcr1_polymerase = graph.input("pcr1_polymerase", "first-PCR polymerase",
                                  MaterialKind.CATALYST)
    pcr1 = graph.product(
        "pcr1_product", "plate-barcoded TCR V-to-C amplicons",
        MaterialKind.MOLECULAR,
        state=molecule("PCR-I products", (role("partial Nextera end"),
            role("TCR V(D)J insert"), PLATE_BARCODE, role("partial Nextera end"))))
    graph.react(
        "pcr1", "Multiplex-PCR TCR V-to-C segments",
        StepKind.PCR,
        ((cdna, InputRole.TEMPLATE), (vc_forward, InputRole.STRUCTURAL),
         (vc_reverse, InputRole.STRUCTURAL), (pcr1_polymerase, InputRole.CATALYST)),
        (pcr1,), geometry=_unavailable(
            pcr1.state,
            "V/C primer identities depend on the receptor segment and the source does not define one universal amplicon."))

    pcr2_forward = graph.input("pcr2_forward", "full-adapter PCR-II forward primer",
                               MaterialKind.OLIGO)
    pcr2_reverse = graph.input("pcr2_reverse", "well-UDI PCR-II reverse primer",
                               MaterialKind.OLIGO)
    pcr2_polymerase = graph.input("pcr2_polymerase", "second-PCR polymerase",
                                  MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "well-indexed TIRTL-seq library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(
        "pcr2", "Add full adapters and a well-specific UDI pair",
        StepKind.PCR,
        ((pcr1, InputRole.TEMPLATE), (pcr2_forward, InputRole.STRUCTURAL),
         (pcr2_reverse, InputRole.STRUCTURAL), (pcr2_polymerase, InputRole.CATALYST)),
        (library,), geometry=_unavailable(
            final, "The page records the UDI and adapter roles but not a universal fixed primer pair."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    graph.declare_reviewed_scope((
        "replicate_well_cells", "lysis_condition", "rt_primer_pool",
        "reverse_transcriptase", "reverse_transcription", "vc_forward_pool",
        "vc_reverse_pool", "pcr1_polymerase", "pcr1", "pcr2_forward",
        "pcr2_reverse", "pcr2_polymerase", "pcr2", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'tirtl-seq'
    return _tirtl_graph(d, final)

WELL_UDI = role("well unique dual index", tag="cbc",
                molecular_feature=feature("well_udi", "sample_index", "unknown"))
PLATE_BARCODE = role("plate barcode", tag="cbc",
                     molecular_feature=feature("plate_barcode", "sample_index", "unknown"))

def _tirtl():
    return chain(molecule("cells distributed across replicate wells", (
        role("T-cell RNA"),)), (
        ("Lyse and reverse-transcribe in each miniaturized well",
         molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))),
         "Alpha/beta pairing is reconstructed from shared occurrence patterns across replicate wells."),
        ("Multiplex-PCR V-to-C segments with partial Nextera-tailed primers",
         molecule("PCR-I products", (role("partial Nextera end"), role("TCR V(D)J insert"),
                                      PLATE_BARCODE, role("partial Nextera end"))),
         "The C-segment reverse primer carries the plate barcode; alpha and beta primer pools are amplified together."),
        ("PCR-II adds full adapters and a well-specific UDI pair",
         molecule("TIRTL-seq library", (role("P5 / Read 1 end"), WELL_UDI,
                                         role("TCR amplicon"), PLATE_BARCODE,
                                         role("Read 2 / P7 end"))),
         "The second PCR makes 384 well-addressable libraries that can be pooled for Illumina sequencing."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library("TIRTL-seq final library", (WELL_UDI, role("TCR V(D)J amplicon"), PLATE_BARCODE))
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("The defining source establishes adapter and identifier roles, but the page does not promote panel-specific or supplementary oligo sequences into universal protocol constants."),
    workflow=_tirtl, oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
