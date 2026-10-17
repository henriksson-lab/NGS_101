"""Local source and reaction model for raisin-rna-seq."""
from chemdraw import Construct, Scene, Segment, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from batch_ngs import nextera_library, seg
from protocol_model import role
from rare_cell_and_affinity import _scene, _workflow
from types import SimpleNamespace

TRANSCRIPT = role("transcript-derived insert")
RAISIN_CST = (("CHAPS", "0.49% w/v"), ("NaCl", "146 mM"),
              ("CaCl2", "1 mM"), ("MgCl2", "21 mM"),
              ("Tris pH 8.0", "10 mM"))
RAISIN_RT_PRIMER = "AAGCAGTGGTATCAACGCAGAGTACT(30)VN"
RAISIN_RT_HANDLE = RAISIN_RT_PRIMER.removesuffix("T(30)VN")
RAISIN_LOCKED_TSO = "AAGCAGTGGTATCAACGCAGAGTACrGrG+G"
RAISIN_PCR_PRIMER = "AAGCAGTGGTATCAACGCAGAGT"
RAISIN_OLIGO_SOURCE = (
    "Drokhlyansky et al., Cell 2020, PMC8358727 Key Resources Table, Oligonucleotides"
)

def raisin_particle() -> Scene:
    sc = Scene()
    sc.strand("nucleus", [role("intact nucleus")], label="RAISIN")
    sc.strand("ribosome", [role("rough ER / attached ribosomes"),
                            role("ribosome-bound poly(A) RNA")], label="retained exterior RNA")
    sc.note("ribosome", "CHAPS/salt extraction retains nucleus-associated ribosomes and mRNA")
    return sc

def smart_cdna_scene(*, switched: bool = True) -> Scene:
    sc = Scene()
    parts = [Segment("SMART RT handle", RAISIN_RT_HANDLE, "ad"),
             Segment("oligo(dT)30", "T" * 30),
             Segment("anchored VN", "VN", placeholder=True), TRANSCRIPT]
    if switched:
        parts.append(role("copied locked TSO"))
    sc.strand("cDNA", parts, label="full-length cDNA" if switched else "first-strand cDNA")
    sc.note("cDNA", "source-printed RT primer; template-switch product shown by role"
            if switched else "source-printed oligo(dT) RT primer")
    return sc

def raisin_library() -> tuple[Construct, tuple]:
    return nextera_library([seg("RAISIN transcript-derived insert", "X" * 36,
                                placeholder=True)], "RAISIN RNA-seq library")

def raisin_workflow() -> Workflow:
    final, _ = raisin_library()
    return _workflow(_scene("frozen tissue", role("frozen tissue")), (
        ("Extract with CST (RAISIN) buffer", raisin_particle(),
         "CHAPS with salt, Ca²⁺ and Mg²⁺ preserves intact nuclei with attached ribosomes."),
        ("Sort one RAISIN per well, lyse and purify RNA", _scene(
            "single-RAISIN RNA", role("nuclear RNA"), role("ribosome-bound poly(A) RNA")),
         "The plate implementation sorts into TCL/β-mercaptoethanol and uses RNAClean XP."),
        ("Oligo-dT reverse-transcribe with a locked TSO, then amplify", smart_cdna_scene(),
         "Maxima H Minus reverse transcriptase is followed by 21 cycles of KAPA HiFi PCR; cDNA is purified twice with AMPure XP."),
        ("Tagment amplified cDNA with Nextera XT and add sample indexes",
         Scene.duplex(list(final), label="RAISIN RNA-seq library"),
         "The custom barcode-adapter sequences were available only on request. The cited Shekhar et al. procedure used 0.075 ng cDNA and ¼ Nextera XT volumes."),
    ))

final, primers = raisin_library()
SPEC = SimpleNamespace(
    TITLE="RAISIN RNA-seq — intact nuclei with ribosome-bound RNA",
    NOTES="01_raisin_rna_seq.html",
    SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
    SUMMARY="CST extraction preserves an intact nucleus together with attached ribosomes and their poly(A) RNA; the plate implementation uses SMART-seq2 followed by Nextera XT.",
    CAVEAT="The paper prints the RT, locked TSO and cDNA PCR oligos in its Key Resources Table. It withholds the custom Nextera XT barcode-adapter sequences, so their individual indexes remain symbolic. No read configuration is printed beyond NextSeq 500 with a 75-cycle high-output v2 kit.",
    FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers,
    READ_LENGTHS=None,
    FINAL_CAPTION="A Nextera XT library made from plate-amplified RAISIN cDNA.",
    SEQUENCING_INTRO="Standard Nextera primers bind the completed library as computed below.",
    SEQUENCING_ENDING="",
    SEQUENCING_UNAVAILABLE="",
    workflow=raisin_workflow,
    oligos=lambda: (
                f"5′-{RAISIN_RT_PRIMER}-3′ — SMART oligo(dT) reverse-transcription primer",
                f"5′-{RAISIN_LOCKED_TSO}-3′ — locked SMART template-switch oligo (rG, +G as printed)",
                f"5′-{RAISIN_PCR_PRIMER}-3′ — SMART cDNA amplification primer"),
)

def _raisin_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    tissue = graph.input("tissue", "frozen tissue", MaterialKind.BIOLOGICAL,
                         state=_scene("frozen tissue", role("frozen tissue")), primary=True)
    cst = graph.input("cst_buffer", "CST RAISIN extraction buffer",
                      MaterialKind.CONDITION)
    particles = graph.product("raisin_particles", "intact nuclei retaining attached ribosomes",
                              MaterialKind.MOLECULAR, state=raisin_particle())
    graph.react("extract", "Extract RAISIN particles with CST buffer",
                StepKind.EXTRACTION,
                ((tissue, InputRole.TEMPLATE), (cst, InputRole.CONDITION)), (particles,))
    sorter = graph.input("sort_rule", "single-RAISIN fluorescence-sorting rule",
                         MaterialKind.SELECTION)
    selected = graph.product("single_raisin", "one sorted RAISIN per well",
                             MaterialKind.MOLECULAR, state=particles.state)
    excluded = graph.product("excluded_particles", "unsorted or multiplet particles",
                             MaterialKind.WASTE)
    graph.react("sort", "Sort one RAISIN into each well", StepKind.SAMPLING,
                ((particles, InputRole.TEMPLATE), (sorter, InputRole.SELECTION)),
                (selected, excluded), changes_molecule=False)
    graph.discard(excluded, "does not meet the one-particle-per-well selection")
    rnaclean = graph.input("rnaclean_beads", "Agencourt RNAClean XP beads", MaterialKind.SELECTION)
    purified = graph.product("purified_rna", "bead-purified single-RAISIN RNA",
                             MaterialKind.MOLECULAR, state=selected.state)
    lysate_waste = graph.product("lysate_supernatant", "TCL lysate supernatant and bead washes",
                                 MaterialKind.WASTE)
    graph.react("rna_cleanup", "Purify well RNA with RNAClean XP beads", StepKind.CLEANUP,
                ((selected, InputRole.TEMPLATE), (rnaclean, InputRole.SELECTION)),
                (purified, lysate_waste))
    graph.discard(lysate_waste, "unbound lysate removed before reverse transcription")
    rt_primer = graph.input("rt_primer", "oligo(dT) SMART reverse-transcription primer",
                            MaterialKind.OLIGO, state=RAISIN_RT_PRIMER,
                            source=RAISIN_OLIGO_SOURCE)
    rtase = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase (EP0753)",
                        MaterialKind.CATALYST)
    first = graph.product("first_strand", "RAISIN transcript first-strand cDNA",
                          MaterialKind.MOLECULAR, state=smart_cdna_scene(switched=False))
    graph.react("reverse_transcribe", "Oligo-dT reverse-transcribe RAISIN RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((purified, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first,))
    tso = graph.input("locked_tso", "locked SMART template-switch oligo",
                      MaterialKind.OLIGO, state=RAISIN_LOCKED_TSO,
                      source=RAISIN_OLIGO_SOURCE)
    switched = graph.product("switched", "template-switched full-length RAISIN cDNA",
                             MaterialKind.MOLECULAR, state=smart_cdna_scene())
    graph.react("template_switch", "Template-switch at the transcript 5-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (switched,))
    pre_f = graph.input("preamp_forward", "SMART cDNA preamplification forward primer",
                        MaterialKind.OLIGO, state=RAISIN_PCR_PRIMER,
                        source=RAISIN_OLIGO_SOURCE)
    pre_r = graph.input("preamp_reverse", "SMART cDNA preamplification reverse primer",
                        MaterialKind.OLIGO, state=RAISIN_PCR_PRIMER,
                        source=RAISIN_OLIGO_SOURCE)
    pre_poly = graph.input("preamp_polymerase", "KAPA HiFi HotStart ReadyMix",
                           MaterialKind.CATALYST)
    amplified = graph.product("amplified", "21-cycle amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=switched.state)
    graph.react("preamp", "Amplify full-length cDNA for 21 cycles", StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (pre_f, InputRole.STRUCTURAL),
                 (pre_r, InputRole.STRUCTURAL), (pre_poly, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="The paper prints one SMART cDNA amplification primer, but the exact template binding positions cannot be placed on the symbolic transcript-derived intermediate.",
                    symmetric_species=True))
    ampure = graph.input("ampure_beads", "Agencourt AMPure XP beads (two purifications)",
                         MaterialKind.SELECTION)
    clean_cdna = graph.product("purified_cdna", "twice AMPure-purified amplified cDNA",
                               MaterialKind.MOLECULAR, state=switched.state)
    pcr_waste = graph.product("preamp_supernatant", "PCR supernatant, primers and bead washes",
                              MaterialKind.WASTE)
    graph.react("cdna_cleanup", "Purify amplified cDNA twice with AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (clean_cdna, pcr_waste))
    graph.discard(pcr_waste, "unbound reaction components removed by bead cleanup")
    transposome = graph.input("transposome", "Nextera XT loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented amplified RAISIN cDNA",
                           MaterialKind.MOLECULAR, state=amplified.state)
    waste = graph.product("tn5_waste", "same-end or non-library tagmentation products",
                          MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified RAISIN cDNA", StepKind.TAGMENTATION,
                ((clean_cdna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not form a heterologous amplifiable Nextera library")
    forward = graph.input("library_forward", "custom Nextera XT barcode adapter primer (i5 side; sequence on request)",
                          MaterialKind.OLIGO)
    reverse = graph.input("library_reverse", "custom Nextera XT barcode adapter primer (i7 side; sequence on request)",
                          MaterialKind.OLIGO)
    polymerase = graph.input("library_polymerase", "Nextera XT PCR master mix", MaterialKind.CATALYST)
    library = graph.product("library_library", "RAISIN RNA-seq sequencing library",
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    graph.react("library_pcr", "Amplify and index tagmented cDNA with custom barcode adapters",
                StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    spec.FINAL_LIBRARY, reason="The custom barcode-adapter sequences are stated to be available upon request (PMC8358727 Methods); the indexing primer pair is not printed."))
    graph.sequence(library, SequencingPlan(library=spec.FINAL_LIBRARY, primers=tuple(spec.SEQ_PRIMERS)))
    reviewed = ["tissue", "cst_buffer", "extract", "sort_rule", "sort", "rnaclean_beads",
                "rna_cleanup", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "locked_tso", "template_switch",
                "preamp_forward", "preamp_reverse", "preamp_polymerase", "preamp",
                "ampure_beads", "cdna_cleanup",
                "transposome", "tn5", "tagment", "library_forward", "library_reverse",
                "library_polymerase", "library_pcr", "library_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _raisin_graph(SPEC)

SPEC.protocol_graph = protocol_graph
