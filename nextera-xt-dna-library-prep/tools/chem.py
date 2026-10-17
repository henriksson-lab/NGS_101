"""Local source and reaction model for nextera-xt-dna-library-prep."""
from chemdraw import Construct, Scene, Segment, Workflow
from end_combinations import LoadedEnd, LoadedEndPool
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from illumina_vendor_workflows import Definition, _duplex, _nextera_final, _sequence, _tagged_intermediate, _workflow, nextera_final_rows
import illumina as il
import nextera as nx

ADAPTER_GUIDE = ("Illumina Adapter Sequences 1000000002694 v22, 'Sequences for Illumina Prep and "
                 "Nextera Kits' (extracted text lines 266-283) and Nextera DNA Indexes N7xx/S5xx "
                 "tables (lines 13318-13345, 13505-13521)")
# Adapter guide v22 lines 266-271: "Nextera Transposase Adapters ... used for Nextera
# tagmentation", Read 1 5' TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG and Read 2 5'
# GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG.  Lines 277-283: PCR Primers, Index 1 Read
# 5' CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG, Index 2 Read 5'
# AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC.  N7xx/S5xx indexes are 8 bases.
INDEX_NT = 8
I7_PRIMER = il.P7 + "N" * INDEX_NT + nx.S7
I5_PRIMER = il.P5 + "N" * INDEX_NT + nx.S5
from types import SimpleNamespace


def _nextera_xt_workflow(final: Construct) -> Workflow:
    tagged = _tagged_intermediate("small-genome / amplicon insert")
    return _workflow("1 ng double-stranded DNA", _duplex("small-genome or amplicon DNA"), (
        ("Solution-transposome tagmentation", "tagged fragments", tagged,
         "TD buffer + 1 ng DNA (5 µl) + 5 µl ATM; TAG program 55 °C for 5 min. The engineered transposome fragments DNA and transfers partial adapters in one reaction."),
        ("Neutralize transposase", "neutralized tagged fragments", tagged,
         "Add 5 µl NT immediately at 10 °C (the transposome is still active), then 5 min at room temperature."),
        ("Limited-cycle PCR with two index primers", "dual-index sequencing library",
         Scene.duplex(list(final), label="Nextera XT library"),
         "5 µl each i7/i5 index adapter + 15 µl NPM; 72 °C 3 min, 95 °C 30 s, 12 cycles of 95/55/72 °C, 72 °C 5 min. PCR amplifies the insert and adds P5, P7, i5 and i7 sequences."),
        ("Single-sided bead cleanup", "purified double-stranded library",
         Scene.duplex(list(final), label="Nextera XT library"),
         "50 µl PCR + 30 µl IPB (3:2; amplicon inputs adjust the ratio); fragments < 200 bp, residual adapter and primer dimers are discarded."),
    ))

def _nextera_xt_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "1-ng double-stranded small-genome or amplicon DNA",
                      MaterialKind.BIOLOGICAL,
                      state=_duplex("small-genome or amplicon DNA"), primary=True)
    transposome = graph.input(
        "transposome", "ATM solution transposome carrying partial Nextera adapters",
        MaterialKind.ADAPTER, inferred=True,
        chemical_state=LoadedEndPool(
            (LoadedEnd("Read 1 transposase adapter", nx.S5 + nx.ME),
             LoadedEnd("Read 2 transposase adapter", nx.S7 + nx.ME)), nx.ME),
        source=ADAPTER_GUIDE + " prints the two Nextera transposase adapters; the XT guide "
        "does not list the ATM loaded strands, and the non-transferred strand is not printed.")
    tn5 = graph.input("tn5", "engineered transposase (in ATM)", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "heterologously tagged DNA fragments",
                           MaterialKind.MOLECULAR,
                           state=_tagged_intermediate("small-genome / amplicon insert"))
    unusable = graph.product("unusable_tn5_products", "same-end or non-library tagmentation products",
                             MaterialKind.WASTE)
    graph.react("tagment", "Tagment 1 ng DNA with ATM in TD buffer, 55 °C 5 min",
                StepKind.TAGMENTATION,
                ((dna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, unusable))
    graph.discard(unusable, "cannot bridge the heterologous P5/P7 library arms")
    neutralizer = graph.input("neutralizer", "transposase-neutralization condition",
                              MaterialKind.CONDITION)
    neutral = graph.product("neutral", "neutralized tagged fragments",
                            MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("neutralize", "Add NT immediately at 10 °C; 5 min room temperature", StepKind.QUENCHING,
                ((tagged, InputRole.TEMPLATE), (neutralizer, InputRole.CONDITION)),
                (neutral,))
    reviewed = ["dna", "transposome", "tn5", "tagment", "neutralizer", "neutralize"]
    i7_primer = PCRPrimer("Index 1 (i7) adapter", I7_PRIMER, anneal_nt=len(nx.S7))
    i5_primer = PCRPrimer("Index 2 (i5) adapter", I5_PRIMER, anneal_nt=len(nx.S5))
    i5 = graph.input("library_forward", "Index 2 (i5) adapter, S5xx (8-nt index; kit/well chosen per sample)",
                     MaterialKind.OLIGO, source=ADAPTER_GUIDE,
                     state=Construct([Segment("Index 2 (i5) adapter", I5_PRIMER, placeholder=True)]))
    i7 = graph.input("library_reverse", "Index 1 (i7) adapter, N7xx (8-nt index; kit/well chosen per sample)",
                     MaterialKind.OLIGO, source=ADAPTER_GUIDE,
                     state=Construct([Segment("Index 1 (i7) adapter", I7_PRIMER, placeholder=True)]))
    npm = graph.input("library_polymerase", "NPM (Nextera PCR Master Mix)", MaterialKind.CATALYST)
    amplified = graph.product("library_amplified", "amplified dual-indexed library",
                              MaterialKind.MOLECULAR, state=library)
    graph.react("library_pcr", "NXT PCR: 72 °C 3 min, then 12 cycles adding i7/i5 adapters",
                StepKind.PCR,
                ((neutral, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (npm, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.exact(
                    library, i5_primer, i7_primer,
                    statement="Adapter Sequences v22 prints the Nextera Index 1/Index 2 PCR primers with 8-nt N7xx/S5xx index windows."))
    ipb = graph.input("ipb", "IPB (Illumina Purification Beads), 30 µl per 50 µl PCR",
                      MaterialKind.SELECTION)
    product = graph.product("library_library", "Nextera XT sequencing library",
                            MaterialKind.LIBRARY, state=library)
    short = graph.product("short_fragments", "fragments < 200 bp, residual adapter and primer dimers",
                          MaterialKind.WASTE)
    graph.react("cleanup", "Single-sided IPB bead cleanup, two 80% EtOH washes, elute in RSB",
                StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE), (ipb, InputRole.SELECTION)),
                (product, short), changes_molecule=False)
    graph.discard(short, "left in the bead supernatant (technical note 770-2015-015 p. 5)")
    _sequence(graph, product, library, primers, d.sequencing_intro)
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase", "library_pcr",
                     "library_amplified", "ipb", "cleanup", "library_library", "short_fragments"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "Nextera XT DNA Library Prep", "01_nextera_xt.html",
        'Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nextera-xt/15031942_07_nextera-xt-product-documentation.pdf">Nextera XT DNA Library Prep Product Documentation 15031942 v07</a>.',
        "A solution transposome fragments and tags 1 ng DNA; neutralization and limited-cycle PCR complete a dual-index Nextera library.",
        "The XT guide describes partial-adapter transfer but does not list the ATM loaded strands; the Read 1/Read 2 transposase adapters drawn here are the generic Nextera tagmentation adapters printed in Illumina Adapter Sequences v22, and the non-transferred strand is not printed in the inspected documents. The tagmented intermediate is drawn as a closed duplex; the guide's NXT PCR program opens with 72 °C for 3 min but does not describe the transposition gaps or their filling. Public Nextera mosaic-end, primer-entry and sequencing-primer sequences are imported from the central canonical model. The drawing shows the double-stranded library before optional bead normalization.",
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
