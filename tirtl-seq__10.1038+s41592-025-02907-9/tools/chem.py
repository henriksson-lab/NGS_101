from types import SimpleNamespace
from chemdraw import feature
from protocol_model import chain

from readout_protocols import Construct, Definition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, _authority, _unavailable, molecule, role
DEFINITION = Definition("TIRTL-seq — throughput-intensive rapid TCR library sequencing",
        "10.1038/s41592-025-02907-9", "Pogorelyy et al., Nature Methods (2026)",
        "Miniaturized plate RT, multiplex V/C amplification, plate barcodes and well-specific UDIs generate deep alpha/beta TCR repertoires whose chain pairs are inferred across replicate wells.",
        "TCR insert reads are associated with the PCR-I plate barcode and PCR-II well UDI pair.",
        "01_tirtl_seq.html")


METHODS = "Methods, TIRTL-seq library preparation (PMC12791010 JATS)"
SUPP1 = "Supplementary Data 1 (41592_2025_2907_MOESM3_ESM.xls)"


def _tirtl_graph(d: Definition, final: Construct) -> ProtocolGraph:
    source = _authority(d)
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "replicate_well_cells", "T cells dispensed across 384 replicate wells",
        MaterialKind.BIOLOGICAL,
        state=molecule("cells distributed across replicate wells", (role("T-cell RNA"),)),
        primary=True, source=f"{source}; {METHODS}")
    lysis = graph.input("lysis_condition",
                        "0.1% Triton X-100 lysis/RT mix under Vapor-Lock overlay",
                        MaterialKind.CONDITION, source=f"{source}; {METHODS}")
    rt_primers = graph.input("rt_primer_pool",
                             "random hexamer and oligo(dT) RT primers (Thermo; sequences not printed)",
                             MaterialKind.OLIGO, source=f"{source}; {METHODS}")
    reverse_transcriptase = graph.input("reverse_transcriptase",
                                       "Maxima H Minus reverse transcriptase",
                                       MaterialKind.CATALYST, source=f"{source}; {METHODS}")
    cdna = graph.product(
        "tcr_cdna", "first-strand cDNA including TCR alpha/beta",
        MaterialKind.MOLECULAR,
        state=molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))))
    graph.react(
        "reverse_transcription",
        "Lyse and reverse-transcribe in each well (42 C 5 min, 25 C 10 min, 50 C 60 min, 94 C 5 min)",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cells, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
         (rt_primers, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST)), (cdna,),
        source=f"{source}; {METHODS}")

    vc_forward = graph.input("vc_forward_pool",
                             "88 V-segment PCR-I primers with the Nextera Read 2 (P7-side) arm",
                             MaterialKind.OLIGO, source=f"{source}; {SUPP1}, sheet Vmix_PCRI_forward")
    vc_reverse = graph.input("vc_reverse_pool",
                             "TRAC/TRBC PCR-I primers with the Nextera Read 1 (P5-side) arm and a plate barcode",
                             MaterialKind.OLIGO, source=f"{source}; {SUPP1}, sheet C_PCRI_reverse")
    pcr1_polymerase = graph.input("pcr1_polymerase", "KAPA2G Fast Multiplex Mix",
                                  MaterialKind.CATALYST, source=f"{source}; {METHODS}")
    pcr1 = graph.product(
        "pcr1_product", "plate-barcoded TCR C-to-V amplicons",
        MaterialKind.MOLECULAR,
        state=molecule("PCR-I products", PCR1_PARTS))
    graph.react(
        "pcr1", "Multiplex PCR I in the RT well (20 cycles, 59 C anneal)",
        StepKind.PCR,
        ((cdna, InputRole.TEMPLATE), (vc_forward, InputRole.STRUCTURAL),
         (vc_reverse, InputRole.STRUCTURAL), (pcr1_polymerase, InputRole.CATALYST)),
        (pcr1,), geometry=_unavailable(
            pcr1.state,
            "V/C primer identities depend on the receptor segment; the model does not select one panel row as a universal amplicon."),
        source=f"{source}; {METHODS}; {SUPP1}")

    pcr2_forward = graph.input("pcr2_forward",
                               "well-specific i5 PCR-II primer (P5 + 10-nt i5 + Read 1 arm)",
                               MaterialKind.OLIGO, source=f"{source}; {SUPP1}, sheet i5_PCRII_forward")
    pcr2_reverse = graph.input("pcr2_reverse",
                               "well-specific i7 PCR-II primer (P7 + 10-nt i7 + Read 2 arm)",
                               MaterialKind.OLIGO, source=f"{source}; {SUPP1}, sheet i7_PCRII_reverse")
    pcr2_polymerase = graph.input("pcr2_polymerase", "Q5 Hot Start High-Fidelity DNA polymerase",
                                  MaterialKind.CATALYST, source=f"{source}; {METHODS}")
    indexed = graph.product("well_indexed_products", "well-indexed TIRTL-seq amplicons",
                            MaterialKind.MOLECULAR, state=molecule("PCR-II products", FINAL_PARTS))
    graph.react(
        "pcr2", "Indexing PCR II on ~1:10 diluted PCR-I product (15 cycles)",
        StepKind.PCR,
        ((pcr1, InputRole.TEMPLATE), (pcr2_forward, InputRole.STRUCTURAL),
         (pcr2_reverse, InputRole.STRUCTURAL), (pcr2_polymerase, InputRole.CATALYST)),
        (indexed,), geometry=_unavailable(
            final, "The page records the published i5/i7 well-UDI primer roles but not one selected well pair."),
        source=f"{source}; {METHODS}; {SUPP1}")
    beads = graph.input("ampure_1x", "1X AMPure XP beads", MaterialKind.SELECTION,
                        source=f"{source}; {METHODS}")
    library = graph.product("sequencing_library", "purified TIRTL-seq library",
                            MaterialKind.LIBRARY, state=final)
    waste = graph.product("ampure_waste", "unbound primers and short products",
                          MaterialKind.WASTE)
    graph.react("bead_cleanup", "Pool the plate by inverted centrifugation and purify with 1X AMPure XP",
                StepKind.CLEANUP, ((indexed, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (library, waste), changes_molecule=False, source=f"{source}; {METHODS}")
    graph.discard(waste, "AMPure XP supernatant")
    graph.sequence(library, SequencingPlan(limitation=SEQUENCING_LIMITATION))
    graph.declare_reviewed_scope((
        "replicate_well_cells", "lysis_condition", "rt_primer_pool",
        "reverse_transcriptase", "reverse_transcription", "vc_forward_pool",
        "vc_reverse_pool", "pcr1_polymerase", "pcr1", "pcr2_forward",
        "pcr2_reverse", "pcr2_polymerase", "pcr2", "ampure_1x",
        "bead_cleanup", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'tirtl-seq'
    return _tirtl_graph(d, final)

WELL_I5 = role("well i5 index (10 nt)", tag="cbc", molecular_feature=feature(
    "well_udi_i5", "sample_index", "unknown", group="well_udi", part="i5"))
WELL_I7 = role("well i7 index (10 nt)", tag="cbc", molecular_feature=feature(
    "well_udi_i7", "sample_index", "unknown", group="well_udi", part="i7"))
PLATE_BARCODE = role("plate barcode", tag="cbc",
                     molecular_feature=feature("plate_barcode", "sample_index", "unknown"))
# Orientation follows Supplementary Data 1: the C-segment primers carry the Nextera
# Read 1 (TCGTCGGCAGCGTC...) arm and the plate barcode, the V-segment primers the
# Nextera Read 2 (GTCTCGTGGGCTCGG...) arm; i5 primers end in the Read 1 arm, i7 in Read 2.
PCR1_PARTS = (role("Nextera Read 1 arm"), PLATE_BARCODE, role("C-segment primer site"),
              role("TCR CDR3 / V(D)J region (antisense)"), role("V-segment primer site"),
              role("Nextera Read 2 arm"))
FINAL_PARTS = (role("P5"), WELL_I5) + PCR1_PARTS + (WELL_I7, role("P7"))
SEQUENCING_LIMITATION = (
    "The Methods state NovaSeq sequencing at 300,000 reads per well for the 384-well protocol "
    "but do not print read lengths or index-read cycles. The protocols.io 96-well v2 protocol "
    "(dx.doi.org/10.17504/protocols.io.kqdg3k6xpv25/v2, step 50) prefers 2x150 nt for the manual "
    "96-well variant; no acquired cycle counts are assigned here.")

def _tirtl():
    return chain(molecule("cells distributed across replicate wells", (
        role("T-cell RNA"),)), (
        ("Lyse and reverse-transcribe in each well (Maxima H Minus, random hexamer + oligo(dT))",
         molecule("TCR cDNA", (role("TCR-alpha cDNA"), role("TCR-beta cDNA"))),
         "Alpha/beta pairing is reconstructed from shared occurrence patterns across replicate wells."),
        ("Multiplex PCR I with V-segment (Read 2 arm) and plate-barcoded C-segment (Read 1 arm) primers",
         molecule("PCR-I products", PCR1_PARTS),
         "KAPA2G Fast Multiplex, 20 cycles, in the RT well; alpha and beta primer pools are amplified together."),
        ("Dilute ~1:10 and run indexing PCR II with a well-specific i5/i7 UDI pair",
         molecule("PCR-II products", FINAL_PARTS),
         "Q5 Hot Start, 15 cycles; each of 384 wells receives its own i5 and i7 index."),
        ("Pool the plate and purify with 1X AMPure XP",
         molecule("TIRTL-seq library", FINAL_PARTS),
         "The pooled library is sequenced on Illumina NovaSeq."),
    ))

d = DEFINITION
FINAL_LIBRARY = Construct(FINAL_PARTS, name="TIRTL-seq final library")
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint + " " + SEQUENCING_LIMITATION,
    SEQUENCING_UNAVAILABLE=("Supplementary Data 1 publishes the full V/C primer and well-specific UDI families. The examples above identify source-known arms; one exact library depends on the chosen TCR segment, plate and well."),
    workflow=_tirtl, oligos=lambda: (
        "<p>Published Supplementary Data 1, representative PCR I primers (5′→3′):</p>",
        "<p>P7_TRAV01-1: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTCATTCCTTAGTCGCTCTGATAGTTATGGTTA</p>",
        "<p>TIRTL_alpha_P5_ID01: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<cbc>TCCTTCG</cbc>AGCTGGTACACGGCAGGGTC</p>",
        "<p>The marked window is positions 34&ndash;40, the bases that differ among the 12 alpha (and 12 beta) plate-ID rows. The Methods describe 6- or 8-nt plate barcodes; the supplement does not mark the barcode boundary.</p>",
        "<p>Published Supplementary Data 1, example well A1 i5 PCR II primer: AATGATACGGCGACCACCGAGATCTACAC<cbc>TTGGTACGCG</cbc>TCGTCGGCAGCGTC</p>",
        "<p>Published Supplementary Data 1, example well A1 i7 PCR II primer: CAAGCAGAAGACGGCATACGAGAT<cbc>ATAGGCGCTC</cbc>GTCTCGTGGGCTCGG</p>",
    ),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
