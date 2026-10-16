"""CAGEscan source-owned oligos, handling decisions and ordered reactions."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature
import illumina as il
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy,
                            StepKind, Terminal)
from protocol_model import molecule, role
import seqprimers as sp

CAGE_OLIGOS = (
    "5′-TAGTCGAACTGAAGGTCTCCAGCAXnrGrGrG-3′ — CAGEscan template-switch primer (Xn is a 3- or 6-base DNA barcode; terminal rGrGrG is RNA)",
    "5′-TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCTN6-3′ — CAGEscan RT primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGC-3′ — CAGEscan semi-suppressive PCR forward primer",
    "5′-TGACGTCGTCTAGTCGAACTGAAGGTCTCCGAACC-3′ — CAGEscan semi-suppressive PCR reverse primer",
    "5′-AATGATACGGCGACCACCGAGATCTACACTAGTCGAACTGAAGG-3′ — CAGEscan forward PCR primer",
    "5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse PCR primer",
    "5′-TAGTCGAACTGAAGGTCTCCAGCA-3′ — CAGEscan forward paired-end sequencing primer",
    "5′-CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′ — CAGEscan reverse paired-end sequencing primer",
)
ORDERED_OLIGOS = tuple(row.split("5′-", 1)[1].split("-3′", 1)[0] for row in CAGE_OLIGOS)
FINAL_FORWARD = il.P5 + ORDERED_OLIGOS[4][len(il.P5):]
FINAL_REVERSE = il.P7 + ORDERED_OLIGOS[5][len(il.P7):]
FINAL_POOL_UL = 30
FINAL_QC_UL = 5
FINAL_REMAINDER_UL = FINAL_POOL_UL - FINAL_QC_UL


def oligo_state(name, sequence):
    return Construct([Segment(name, sequence)], name=name)


def rt_oligo():
    return Construct([Segment("source RT handle", ORDERED_OLIGOS[1][:-2]),
                      Segment("random-priming N6", "N" * 6, placeholder=True, length_bp=6)], name="CAGEscan RT primer")


def tso_oligo():
    return Construct([Segment("source TSO handle", ORDERED_OLIGOS[0].split("Xn")[0]),
                      role("source DNA barcode, three or six positions", molecular_feature=feature("cagescan_sample_barcode", "sample_index", "unknown")),
                      Segment("RNA rGrGrG tail", "rGrGrG", placeholder=True, length_bp=3, note="Three terminal riboguanosines as ordered")], name="CAGEscan TSO")


def first_strand_pool():
    sc = Scene()
    sc.strand("cap_linked", [role("preferentially cap-linked template-switched cDNA")], label="cap-linked cDNA class")
    sc.strand("other_cdna", [role("other first-strand cDNA")], label="other first-strand class")
    sc.strand("primer_artifacts", [role("primer-derived short templates")], label="primer-artifact class")
    sc.note("primer_artifacts", "Conceptual classes in the same pool; not covalently joined or physically separated.")
    return sc

FIRST_POOL_UL = 30
FIRST_PILOT_UL = 1.2
FIRST_PRODUCTION_UL = 3
FIRST_REMAINDER_UL = FIRST_POOL_UL - FIRST_PILOT_UL - FIRST_PRODUCTION_UL

TITLE = "CAGEscan — paired-end cap-to-transcript sequencing"
NOTES = "01_cagescan.html"
SOURCE = '<a href="https://doi.org/10.1038/nmeth.1470">Plessy et al., Nature Methods (2010)</a>; original Online Methods pp.16–18 and Supplementary Table 8.'
SUMMARY = ("Random-primed reverse transcription and preferential template switching at capped RNA "
           "are followed by RNAClean, pilot qPCR, semi-suppressive PCR, AMPure, a second pilot qPCR, "
           "adapter-replacement PCR and a second AMPure cleanup. CAGEscan does not use the separate nanoCAGE EcoP15I route.")
CAVEAT = ("All eight source oligos, including custom paired-end sequencing primers, are published below. "
          "This model leaves their complete molecular assembly and initial PCR binding geometry unassembled; "
          "that is a model question, not evidence of unpublished primers. The barcode is an inline three- "
          "or six-base DNA sequence in the template-switch oligo, not a generic eight-base instrument index.")
FINAL_LIBRARY = Construct([role("unassembled TSO-derived library end"),
                           role("inline TSO barcode, source three or six bases", molecular_feature=feature("cagescan_sample_barcode", "sample_index", "unknown")),
                           role("cap-linked transcript cDNA"), role("random-priming-defined downstream end"),
                           role("unassembled opposite library end")], name="CAGEscan library roles")
SEQ_PRIMERS = ()
READ_LENGTHS = None
ACQUIRED_ROLES = ("Read 1", "Read 2")
READ_SPECS = (sp.read_spec("Read 1", None, "TSO/cap-proximal end and inline sample barcode; exact assembled boundary under review", "extract_then_include"),
              sp.read_spec("Read 2", None, "random-primer-defined downstream transcript end", "include"))
FINAL_CAPTION = "Source-established library roles only; custom primer sequences are published, but this diagram does not yet assemble exact primer sites."
SEQUENCING_INTRO = ""
SEQUENCING_ENDING = "The original uses Genome Analyzer II with custom forward and reverse paired-end sequencing primers (Table 8). Exact cycle counts and study barcode members remain unassigned here; no separate index read is invented."
SEQUENCING_UNAVAILABLE = "Exact end assembly and primer landings remain model questions despite the source-positive custom sequencing primers."


def oligos():
    return CAGE_OLIGOS


def _state(label, *parts):
    return molecule(label, tuple(role(part) for part in parts))


def _cleanup(graph, prefix, template, reagent):
    beads = graph.input(prefix + "_beads", reagent, MaterialKind.SELECTION)
    kept = graph.product(prefix + "_kept", "cleaned DNA fraction", MaterialKind.MOLECULAR, state=template.state)
    removed = graph.product(prefix + "_removed", "cleanup removed fraction", MaterialKind.WASTE)
    graph.react(prefix, "Clean DNA using " + reagent, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, removed))
    graph.discard(removed, "is separated from the source DNA fraction during bead cleanup")
    return kept


def _pcr(graph, prefix, template, forward, reverse, enzyme, label, *, qc=False, state=None):
    catalyst = graph.input(prefix + "_enzyme", enzyme, MaterialKind.CATALYST)
    out = graph.product(prefix + "_product", "pilot qPCR products" if qc else label,
                        MaterialKind.QC if qc else MaterialKind.MOLECULAR, state=state or template.state)
    graph.react(prefix, label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (catalyst, InputRole.CATALYST)), (out,),
                geometry=PCRGeometry.unavailable(Construct([role("unassembled initial PCR-template pool")], name=prefix + " input geometry"),
                    reason="Source primers are known; the initial template and exact annealing geometry are not assembled in this model."))
    if qc:
        graph.mark_terminal(Terminal(out.id, Disposition.QC, "pilot reaction informs source cycle selection"))
    return out


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "RNA containing capped and uncapped transcripts", MaterialKind.BIOLOGICAL,
                      state=_state("source RNA", "capped and uncapped RNA"), primary=True)
    rt_primer = graph.input("rt_primer", "Table 8 CAGEscan N6 RT primer", MaterialKind.OLIGO, state=rt_oligo())
    rtase = graph.input("rtase", "source SuperScript III or PrimeScript reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("first", "random-primed first-strand cDNA", MaterialKind.MOLECULAR,
                          state=_state("first strand", "RT-primer handle", "RNA-derived cDNA"))
    graph.react("reverse_transcribe", "Random-prime reverse transcription", StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (first,))
    tso = graph.input("tso", "Table 8 TSO: inline three/six-base DNA barcode and terminal rGrGrG", MaterialKind.OLIGO, state=tso_oligo())
    switched = graph.product("switched", "template-switched first-strand pool", MaterialKind.MOLECULAR,
                             state=first_strand_pool())
    graph.react("template_switch", "Template-switch in the source reverse-transcription reaction", StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (rtase, InputRole.CATALYST)), (switched,))
    clean_first = _cleanup(graph, "rna_clean", switched, "Agencourt RNAClean XP")
    sampling = graph.input("first_sampling_rule", "source 30-µL first-strand pool: 1.2-µL pilot and 3-µL production aliquots", MaterialKind.CONDITION)
    pilot1 = graph.product("pilot1", "1.2-µL first-strand pilot aliquot", MaterialKind.MOLECULAR, state=clean_first.state)
    production1 = graph.product("production1", "3-µL first-strand production aliquot", MaterialKind.MOLECULAR, state=clean_first.state)
    remainder1 = graph.product("remainder1", f"remaining first-strand pool ({FIRST_REMAINDER_UL:g} µL in this declared split)", MaterialKind.MOLECULAR, state=clean_first.state)
    graph.react("sample_first", "Take first-strand pilot and production aliquots", StepKind.SAMPLING,
                ((clean_first, InputRole.TEMPLATE), (sampling, InputRole.CONDITION)), (pilot1, production1, remainder1), changes_molecule=False)
    graph.mark_terminal(Terminal(remainder1.id, Disposition.UNRESOLVED, "Subsequent handling of the rest of the named 30-µL pool is not assigned; no storage or discard is invented."))
    semi_f = graph.input("semi_forward", "Table 8 semi-suppressive forward primer", MaterialKind.OLIGO, state=oligo_state("semi-suppressive forward", ORDERED_OLIGOS[2]))
    semi_r = graph.input("semi_reverse", "Table 8 semi-suppressive reverse primer", MaterialKind.OLIGO, state=oligo_state("semi-suppressive reverse", ORDERED_OLIGOS[3]))
    _pcr(graph, "pilot_semi", pilot1, semi_f, semi_r, "SYBR Premix Ex Taq", "Pilot qPCR for semi-suppressive cycle selection", qc=True)
    semi = _pcr(graph, "semi_pcr", production1, semi_f, semi_r, "TaKaRa Ex Taq HS", "Semi-suppressive second-strand synthesis/amplification")
    cleaned_semi = _cleanup(graph, "ampure_first", semi, "Agencourt AMPure XP")
    sampling2 = graph.input("second_sampling_rule", "source 8-ng pilot and remaining 20-ng production split", MaterialKind.CONDITION)
    pilot2 = graph.product("pilot2", "8-ng adapter-PCR pilot aliquot", MaterialKind.MOLECULAR, state=cleaned_semi.state)
    production2 = graph.product("production2", "remaining 20-ng adapter-PCR production aliquot", MaterialKind.MOLECULAR, state=cleaned_semi.state)
    graph.react("sample_second", "Take adapter-PCR pilot and production aliquots", StepKind.SAMPLING,
                ((cleaned_semi, InputRole.TEMPLATE), (sampling2, InputRole.CONDITION)), (pilot2, production2), changes_molecule=False)
    forward = graph.input("forward", "Table 8 CAGEscan final forward PCR primer", MaterialKind.OLIGO, state=oligo_state("final forward", FINAL_FORWARD))
    reverse = graph.input("reverse", "Table 8 CAGEscan final reverse PCR primer", MaterialKind.OLIGO, state=oligo_state("final reverse", FINAL_REVERSE))
    _pcr(graph, "pilot_final", pilot2, forward, reverse, "SYBR Premix Ex Taq", "Pilot qPCR for adapter-replacement amplification", qc=True)
    adapted = _pcr(graph, "adapter_pcr", production2, forward, reverse, "TaKaRa Ex Taq HS", "Replace semi-suppressive ends by paired-end library PCR", state=FINAL_LIBRARY)
    cleaned = _cleanup(graph, "ampure_final", adapted, "Agencourt AMPure XP")
    final_sampling = graph.input("final_sampling_rule", "source 5-of-30-µL gel aliquot", MaterialKind.CONDITION)
    qc = graph.product("gel_qc", "5-µL library gel-check aliquot", MaterialKind.QC, state=FINAL_LIBRARY)
    library = graph.product("library", f"remaining CAGEscan sequencing library ({FINAL_REMAINDER_UL:g} µL before concentration adjustment)", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("sample_final", "Take the source agarose-gel check aliquot", StepKind.SAMPLING,
                ((cleaned, InputRole.TEMPLATE), (final_sampling, InputRole.CONDITION)), (qc, library), changes_molecule=False)
    graph.mark_terminal(Terminal(qc.id, Disposition.QC, "source 1% agarose-gel quality check; not preparative gel selection"))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=ACQUIRED_ROLES, limitation=SEQUENCING_UNAVAILABLE))
    graph.annotate(SourceDiscrepancy("semi_pcr", "Apparent primer concentration unit issue in the source",
                    (SourceClaim("Online Methods printed p.17, CAGEscan cDNA amplification", "100 mM semi-suppressive primers"),
                     SourceClaim("Online Methods printed p.17, adapter-PCR pilot", "200 mM final PCR primers; nearby dNTPs are printed as 200 µM")),
                    "Keep the source-specific operations and oligos without assigning a corrected primer concentration.",
                    "Both source text and the original PDF visibly print mM for primers, distinct from µM for dNTPs.",
                    "The intended primer concentrations are not resolved by this model."))
    graph.annotate(Highlight("semi_pcr", "Preferential amplification is not physical artifact removal",
                            "The source describes shorter/homoadapter templates as less efficiently amplified. They are not declared a discarded physical fraction at PCR; cleanup fates are represented separately."))
    graph.annotate(Highlight("library", "Published custom primers; unfinished molecular assembly",
                            "Both sequencing primers and all RT/PCR oligos are source-positive. The RT and semi-reverse sequences differ upstream of a shorter exact 3′ match; this is retained without claiming a proven oligo typo. Exact initial PCR geometry remains a model question.", slot=AnnotationSlot.SEQUENCING))
    targets = ["rna", "rt_primer", "rtase", "reverse_transcribe", "tso", "template_switch", "first_sampling_rule", "sample_first", "remainder1", "semi_forward", "semi_reverse", "second_sampling_rule", "sample_second", "forward", "reverse", "final_sampling_rule", "sample_final", "gel_qc", "library"]
    for prefix in ("rna_clean", "ampure_first", "ampure_final"):
        targets.extend((prefix, prefix + "_beads"))
    for prefix in ("pilot_semi", "semi_pcr", "pilot_final", "adapter_pcr"):
        targets.extend((prefix, prefix + "_enzyme"))
    graph.declare_reviewed_scope(targets, SOURCE)
    return graph.require_publishable()

SPEC = SimpleNamespace(**{name: globals()[name] for name in (
    "TITLE", "NOTES", "SOURCE", "SUMMARY", "CAVEAT", "FINAL_LIBRARY", "SEQ_PRIMERS", "READ_LENGTHS", "ACQUIRED_ROLES", "READ_SPECS",
    "FINAL_CAPTION", "SEQUENCING_INTRO", "SEQUENCING_ENDING", "SEQUENCING_UNAVAILABLE", "oligos", "protocol_graph")})
