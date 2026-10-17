"""SMART-seq2 molecular declarations, reaction graph, and public page.

Protocol pages that explicitly inherit SMART-seq2 should construct the branch here
instead of restating its oligos or calling them unavailable.
"""
from __future__ import annotations

from chemdraw import (Construct, Scene, Segment, annotation_rows, feature,
                      panel, revcomp)
from illumina import P5, P7, P7_RC
import nextera as nx
from end_combinations import (EndCombinationGeometry, EndProductGroup,
                              LoadedEnd, LoadedEndPool)
from pcr import PCRGeometry, PCRPrimer
import rt
import seqprimers as sp
from page import head, info
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole, MaterialKind, MolecularRelation,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind, Terminal)
from protocol_views import (OligoEntry, OligoSection, PrimerBindingSection,
                            WorkflowSection, assemble_page)

DT_LEN = 30
DT_ANCHOR = "VN"
DT_LINKER = "AC"
TSO_LINKER = "ACAT"
OLIGO_DT = rt.oligo_dt(DT_LEN, DT_ANCHOR, rt.SMART_HANDLE + DT_LINKER)
TSO = rt.SMART_HANDLE + TSO_LINKER + rt.TSO_G_TAIL_LNA
ISPCR = rt.SMART_HANDLE
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
DOI = "10.1038/nmeth.2639"
PAGE_FILE = "smart-seq2.html"
TEMPLATE_SWITCH_CAPTION = (
    "INFERRED — Template-switch at the transcript 5-prime end; the TSO "
    "linker's final AT is the Nature Protocols form, not the defining paper's")


def defining_tso_linker() -> list[Segment]:
    """The 2013 AC bases and later AT extension, with their distinct evidence."""
    if not TSO_LINKER.startswith(DT_LINKER):
        raise ValueError("the later TSO linker must extend the published AC linker")
    return [Segment("AC", DT_LINKER),
            Segment("AT", TSO_LINKER[len(DT_LINKER):], inferred=True)]


def drawn_tso() -> list[Segment]:
    """The 2013 source boundary as segments for the template-switch state."""
    return [Segment("H", rt.SMART_HANDLE, "tso"),
            *defining_tso_linker(), Segment("rGrG", "GGG", "tso")]


def drawn_preamp_template(body_nt: int, poly_a_nt: int) -> Construct:
    """Source-marked cDNA state used by the Smart-seq2 reaction graph."""
    if body_nt < 1 or poly_a_nt != DT_LEN:
        raise ValueError("drawn cDNA needs a body and the published dT30 tail")
    link = defining_tso_linker()
    return Construct([
        Segment("H5", rt.SMART_HANDLE, "tso"),
        Segment("TSO linker", link[0].top),
        Segment("TSO linker AT", link[1].top, inferred=True),
        Segment("GGG", "GGG", "tso"),
        Segment("body", "X" * body_nt, placeholder=True),
        Segment("B", "B", placeholder=True),
        Segment("polyA", "A" * poly_a_nt),
        Segment("GT", revcomp(DT_LINKER)),
        Segment("H3", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="SMART-seq full-length cDNA before ISPCR")


def final_library() -> Construct:
    """Heterologous Nextera fragment after index PCR; no internal UMI."""
    return Construct([
        Segment("Illumina P5", P5, "p5"),
        Segment("i5", "N" * 8, placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("s5", nx.S5, "s5"), Segment("ME", nx.ME, "me"),
        Segment("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        Segment("ME ", nx.ME_RC, "me"), Segment("s7", nx.S7_RC, "s7"),
        Segment("i7", "N" * 8, placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("Illumina P7", P7_RC, "p7"),
    ], name="SMART-seq2 library")


def protocol_graph() -> ProtocolGraph:
    # The published predecessor supplies molecular drawing scenes used by this method.
    import sys
    from pathlib import Path
    family = Path(__file__).resolve().parents[2] / "smart-seq__10.1038+nbt.2282" / "tools"
    sys.path.insert(0, str(family))
    from smartseq_states import BODY, POLYA, gap_filled, mrna, ss2_rt, tagmented
    protocol, doi = "SMART-seq2", DOI
    source = f"Defining source: doi:{doi}"
    graph = ProtocolGraph(protocol, authority=source)
    initial = Scene(); initial.strand("mRNA", mrna())
    mrna = graph.input("mrna", "single-cell poly(A) RNA", MaterialKind.MOLECULAR,
                       state=initial, primary=True)

    priming, first_scene, switch_scene = ss2_rt(protocol)
    rt_label = ("anchored oligo-dT30VN primer, 1 ul of 10 uM (Biomers; printed in Online "
                "Methods), denatured with 1 ul 10 mM dNTPs at 72 degrees C for 3 min")
    tso_label = ("SMART-seq2 LNA template-switch oligo rGrG+G, 1 ul of 10 uM "
                 "(Supplementary Table 4)")
    preamp_template = drawn_preamp_template(BODY, POLYA)

    rt_primer = graph.input("rt_primer", rt_label, MaterialKind.OLIGO,
                            state=priming)
    # Online Methods "Reverse transcription" (deferred-followup-03) with Supplementary
    # Tables 3-5: SSRTII (Invitrogen 18064-014), 1 M betaine, 0.9 ul of 100 mM MgCl2
    # added (Table 3 column: 12 mM), no purification after RT.
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "SuperScript II reverse transcriptase (MMLV-derived)",
        MaterialKind.CATALYST)
    rt_conditions = graph.input(
        "rt_conditions",
        ("SuperScript II first-strand buffer, 0.25 ul 100 mM DTT, 2 ul 5 M betaine, 0.9 ul of 100 mM "
         "MgCl2 added (Supplementary Table 3: 12 mM); 42 degrees C 90 min, then 10 x "
         "(50 degrees C 2 min, 42 degrees C 2 min), 70 degrees C 15 min"),
        MaterialKind.CONDITION)
    first = graph.product("first_strand", "first-strand cDNA with untemplated CCC",
                          MaterialKind.MOLECULAR, state=first_scene)
    graph.react("reverse_transcription",
                "Reverse-transcribe from the anchored poly(dT) primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST),
                 (rt_conditions, InputRole.CONDITION)), (first,))

    tso = graph.input("tso", tso_label, MaterialKind.OLIGO, state=switch_scene)
    switched = graph.product("switched", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR, state=preamp_template)
    switch_label = TEMPLATE_SWITCH_CAPTION
    graph.react("template_switch", switch_label,
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))

    preamp_polymerase = graph.input(
        "preamp_polymerase",
        ("KAPA HiFi HotStart ReadyMix added directly to the RT; 98 degrees C 3 min, 18 cycles "
         "for cells (98 degrees C 15 s, 67 degrees C 20 s, 72 degrees C 6 min), 72 degrees C 5 min"),
        MaterialKind.CATALYST)
    ispcr = graph.input("ispcr", "ISPCR primer (1 ul of 10 uM; Online Methods) at both cDNA ends",
                        MaterialKind.OLIGO, state=ISPCR)
    preamp_geometry = PCRGeometry.symmetric(
        preamp_template, PCRPrimer("ISPCR", rt.SMART_HANDLE),
        statement=("The source-defined SMART handle occurs in opposite orientations "
                   "at the two cDNA ends, so one ISPCR reagent supplies both inward sites."))
    preamp_inputs = ((switched, InputRole.TEMPLATE),
                     (ispcr, InputRole.STRUCTURAL),
                     (preamp_polymerase, InputRole.CATALYST))
    primer_reviewed = ["ispcr"]
    preamplified = graph.product("preamplified", "amplified full-length cDNA",
                                 MaterialKind.MOLECULAR, state=preamp_template)
    graph.react("preamp", ("PCR-amplify full-length cDNA without bead purification "
                           "after RT"), StepKind.PCR,
                preamp_inputs, (preamplified,), geometry=preamp_geometry)
    preamp_beads = graph.input("preamp_beads", "AMPure XP beads 1:1; elution in 15 ul EB",
                               MaterialKind.SELECTION)
    amplified = graph.product("amplified", "AMPure-purified full-length cDNA",
                              MaterialKind.MOLECULAR, state=preamp_template)
    preamp_waste = graph.product("preamp_supernatant",
                                 "bead supernatant: primers, TSO and short products",
                                 MaterialKind.WASTE)
    graph.react("preamp_cleanup", "Purify pre-amplified cDNA with 1:1 AMPure XP",
                StepKind.CLEANUP,
                ((preamplified, InputRole.TEMPLATE), (preamp_beads, InputRole.SELECTION)),
                (amplified, preamp_waste))
    graph.discard(preamp_waste, "bead supernatant removed (Online Methods, PCR preamplification)")

    index_f_primer = PCRPrimer("N5xx", nx.n5xx_primer("N" * 8), anneal_nt=len(nx.S5))
    index_r_primer = PCRPrimer("N7xx", nx.n7xx_primer("N" * 8), anneal_nt=len(nx.S7))
    transposome = graph.input("transposome",
                              "s5/s7-loaded transposome, Nextera DNA sample preparation kit "
                              "(FC-121-1031)",
                              MaterialKind.ADAPTER, chemical_state=LoadedEndPool(
                                  (LoadedEnd("s5", nx.ADAPTOR_S5),
                                   LoadedEnd("s7", nx.ADAPTOR_S7)), nx.ME))
    tn5 = graph.input("tn5", ("Tagment DNA Enzyme 5 ul with 25 ul 2x Tagment DNA Buffer, "
                              "5 ng cDNA in 50 ul, 55 degrees C 5 min"), MaterialKind.CATALYST)
    tagmented_product = graph.product(
        "tagmented", "heterologous s5/s7-tagmented cDNA fragments",
        MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    # Legacy end-class validation uses WASTE for non-selected end pairs;
    # the terminal below independently declares physical recovery unresolved.
    tagment_waste = graph.product(
        "tagment_waste", "candidate s5/s5, s7/s7 and other non-library Tn5 products",
        MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified cDNA with Nextera s5/s7 adapters",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagmented_product, tagment_waste),
                geometry=EndCombinationGeometry(
                    (index_f_primer, index_r_primer),
                    (EndProductGroup("tagmented", (("s5", "s7"),)),
                     EndProductGroup("tagment_waste",
                                     (("s5", "s5"), ("s7", "s7")))),
                    "index_pcr"),
                relations=(MolecularRelation.END_COMBINATIONS,))
    graph.mark_terminal(Terminal(tagment_waste.id, Disposition.UNRESOLVED,
        "Candidate classes lack the intended heterologous P5/P7 library ends. "
        "The two-primer end-class check does not establish failure to amplify or "
        "physical removal by the column or later AMPure cleanup."))

    zymo = graph.input("zymo", "DNA Clean & Concentrator-5 column; elution in 20 ul RSB",
                       MaterialKind.SELECTION)
    tagmented_clean = graph.product("tagmented_clean", "column-purified tagmented cDNA",
                                    MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    zymo_waste = graph.product("zymo_flowthrough", "column flow-through and washes",
                               MaterialKind.WASTE)
    graph.react("tagment_cleanup", "Purify tagmented cDNA on a Zymo DCC-5 column",
                StepKind.CLEANUP,
                ((tagmented_product, InputRole.TEMPLATE), (zymo, InputRole.SELECTION)),
                (tagmented_clean, zymo_waste))
    graph.discard(zymo_waste, "column flow-through and washes discarded")

    gap_polymerase = graph.input("gap_polymerase",
                                 "Nextera PCR Master Mix (NPM, 15 ul) polymerase",
                                 MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterologous Nextera fragments",
                           MaterialKind.MOLECULAR, state=gap_filled())
    graph.react("gap_fill", ("Fill both nine-base Tn5 gaps: 72 degrees C 3 min opening "
                             "the enrichment PCR program"),
                StepKind.EXTENSION,
                ((tagmented_clean, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))

    final = final_library()
    index_forward = graph.input("index_forward",
                                "Nextera Index Kit (FC-121-1011) N5xx index 2 primer, 5 ul",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse",
                                "Nextera Index Kit (FC-121-1011) N7xx index 1 primer, 5 ul",
                                MaterialKind.OLIGO)
    ppc = graph.input("ppc", "Nextera PCR Primer Cocktail (PPC), 5 ul (kit primers; bases not printed)",
                      MaterialKind.OLIGO)
    index_polymerase = graph.input(
        "index_polymerase",
        "NPM; 98 degrees C 30 s, then 5 cycles (98 degrees C 10 s, 63 degrees C 30 s, 72 degrees C 3 min)",
        MaterialKind.CATALYST)
    index_geometry = PCRGeometry.exact(
        final, index_f_primer, index_r_primer,
        statement=("The canonical N5xx and N7xx primers face inward from the s5 and s7 "
                   "entry points and install the indexed P5/P7 ends."))
    # Supplementary Table 5 costs a "TruSeq Dual Index sequencing primer kit for
    # Single Read runs" per Smart-seq2 sample: Read 1 plus both index reads.
    primers = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2"))
    pcr_library = graph.product("pcr_library", protocol + " enrichment-PCR product",
                                MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Amplify and install indexed P5/P7 ends",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (ppc, InputRole.STRUCTURAL),
                 (index_polymerase, InputRole.CATALYST)),
                (pcr_library,), geometry=index_geometry)
    final_beads = graph.input("final_beads", "AMPure XP beads 1:1", MaterialKind.SELECTION)
    main_library = graph.product("library", protocol + " sequencing-library molecules",
                                 MaterialKind.LIBRARY, state=final)
    final_waste = graph.product("final_supernatant", "bead supernatant after enrichment PCR",
                                MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify the library with 1:1 AMPure XP", StepKind.CLEANUP,
                ((pcr_library, InputRole.TEMPLATE), (final_beads, InputRole.SELECTION)),
                (main_library, final_waste))
    graph.discard(final_waste, "bead supernatant removed")
    # Online Methods "Read alignments": Nextera dual indexes (i7+i5), HiSeq 2000,
    # "43bp reads after demultiplexing"; index cycle counts are not printed.
    graph.sequence(main_library, SequencingPlan(
        library=final, primers=primers,
        acquired_roles=tuple(p.role for p in primers),
        limitation="The article reports 43-base reads after demultiplexing, and the deposit "
                   "contains 43-base single reads. Acquired Read 1 and index cycle counts "
                   "are not established by those processed/deposited lengths."))

    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "rt_conditions",
        "reverse_transcription",
        "tso", "template_switch", *primer_reviewed, "preamp_polymerase", "preamp",
        "preamp_beads", "preamp_cleanup", "transposome", "tn5", "tagment",
        "zymo", "tagment_cleanup", "gap_polymerase", "gap_fill",
        "index_forward", "index_reverse", "ppc", "index_polymerase", "index_pcr",
        "final_beads", "final_cleanup", "library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "template_switch", "Two published TSO linkers",
        (SourceClaim("Picelli et al. 2013, Supplementary Tables 2 and 4",
                     "The selected TSO is the 23-nt handle + AC + rGrG+G."),
         SourceClaim("Picelli et al. 2014 Nature Protocols (doi:10.1038/nprot.2014.006), "
                     "Reagent Setup",
                     "The Smart-seq2 TSO is 5'-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3'.")),
        "The drawn TSO uses ACAT.",
        ("The 2013 Online Methods print the oligo-dT and ISPCR primer but refer to a "
         "supplementary table for the TSO; the 2014 Nature Protocols Reagent Setup prints the "
         "longer form, which descendant protocols use. The 2013 article's pointer names "
         "Supplementary Table 1, while the TSO sequences are in Supplementary Table 2."),
        "The added AT is absent from the defining paper and is drawn inferred."))
    graph.annotate(Highlight(
        "index_pcr", "Kit identities and read plan",
        "The Online Methods use the Nextera DNA Sample Preparation kit (5 ng cDNA, 55 &deg;C "
        "5 min), a Zymo DNA Clean &amp; Concentrator-5 column, then enrichment PCR with "
        "NPM, N7xx, N5xx and PPC (72 &deg;C 3 min, 98 &deg;C 30 s, 5 cycles) and a 1:1 "
        "AMPure XP cleanup. Libraries were sequenced &ldquo;with Nextera dual indexes "
        "(i7+i5) on Illumina HiSeq 2000, giving 43bp reads after demultiplexing&rdquo;. "
        "Supplementary Table 5 costs a TruSeq Dual Index sequencing primer kit for Single "
        "Read runs, and all 76 GSE49321 runs (ENA PRJNA213629) are SINGLE-layout 43-base "
        "reads. This supports single-read acquisition with both indexes and no Read 2, "
        "but does not establish the acquired Read 1 cycle count. Index "
        "lengths of the Nextera 24-index kit are not printed in the sources read."))
    graph.annotate(Highlight(
        "tagment", "Selection after Tn5",
        "Heterologous s5/s7 fragments have both index-PCR sites. Same-end s5/s5 "
        "and s7/s7 candidate products are checked against the declared primer pair; "
        "their amplification and physical recovery fates remain unresolved."))
    return graph.require_publishable()

def preamp_template(insert_nt: int = 36, *, linker_inferred: bool = True) -> Construct:
    """Full-length cDNA before ISPCR.

    The ACAT linker is inferred when a descendant names Smart-seq2 without printing
    its TSO. Pass ``linker_inferred=False`` only when that source prints ACAT itself.
    """
    if insert_nt < 1:
        raise ValueError("SMART-seq2 cDNA insert must be non-empty")
    return Construct([
        Segment("5′ SMART handle", rt.SMART_HANDLE, "tso"),
        Segment("TSO linker", TSO_LINKER, inferred=linker_inferred),
        Segment("TSO G contribution", "GGG", "tso"),
        Segment("transcript cDNA", "X" * insert_nt, placeholder=True),
        Segment("poly(A) copy", "A" * DT_LEN),
        Segment("oligo-dT linker reverse complement", revcomp(DT_LINKER)),
        Segment("3′ SMART handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="SMART-seq2 full-length cDNA before ISPCR")


def ispcr_geometry(template: Construct | None = None) -> PCRGeometry:
    template = template or preamp_template()
    return PCRGeometry.symmetric(
        template, PCRPrimer("ISPCR", ISPCR),
        statement=("SMART-seq2 places the same 23-nt ISPCR handle at both cDNA ends, "
                   "so one primer species supplies both inward-facing sites."))


def nextera_index_geometry(final_library: Construct) -> PCRGeometry:
    """Invariant N5xx/N7xx landing sites; representative indices do not affect binding."""
    return PCRGeometry.exact(
        final_library,
        PCRPrimer("N5xx", nx.n5xx_primer("A" * 8), anneal_nt=len(nx.S5)),
        PCRPrimer("N7xx", nx.n7xx_primer("A" * 8), anneal_nt=len(nx.S7)),
        statement=("Nextera index primers terminate in s5 and s7; those 3′ regions face "
                   "inward on unlike-ended tagmented cDNA. Representative A8 index "
                   "windows show the invariant landing geometry."))


def oligo_section() -> OligoSection:
    """The page's oligos, assembled from this protocol's molecular constructors."""
    def seg(sequence: str, tag: str | None = None, **kw) -> Segment:
        return Segment("", sequence, tag, **kw)
    rows = [
        OligoEntry("oligo-dTVN", (seg(rt.SMART_HANDLE, "tso"), seg(DT_LINKER),
                                   seg("T" * DT_LEN, placeholder=True),
                                   seg(DT_ANCHOR, placeholder=True))),
        OligoEntry("Template Switching Oligo (TSO)",
                   (seg(rt.SMART_HANDLE, "tso"),
                    *[seg(s.top, s.tag, inferred=s.inferred) for s in defining_tso_linker()],
                    seg(rt.TSO_G_TAIL_LNA, placeholder=True))),
        OligoEntry("ISPCR primer", (seg(rt.SMART_HANDLE, "tso"),)),
        OligoEntry("Nextera mosaic end (ME)", (seg(nx.ME, "me"),)),
        OligoEntry("Nextera N/S5xx entry point (s5)", (seg(nx.S5, "s5"),)),
        OligoEntry("Nextera N7xx entry point (s7)", (seg(nx.S7, "s7"),)),
        OligoEntry("Nextera Index Kit N/S5xx index primer",
                   (seg(P5, "p5"), seg("[8-bp i5]", placeholder=True),
                    seg(nx.S5, "s5"))),
        OligoEntry("Nextera Index Kit N7xx index primer",
                   (seg(P7, "p7"), seg("[8-bp i7]", placeholder=True),
                    seg(nx.S7, "s7"))),
    ]
    return OligoSection(
        "Adapter and primer sequences", tuple(rows), multiline=True,
        note_html="""The oligo-dT primer, TSO and ISPCR primer share the same 23-nt
            handle. <code>rG</code> denotes RNA. <code>+G</code> denotes the SMART-seq2 LNA.""")


def final_section(graph: ProtocolGraph) -> str:
    con = graph.sequencing["library"].library
    duplex = Scene.duplex(con.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    rows = duplex.rows() + annotation_rows(con, prefix_width=4)
    return ("<h3>Index PCR, AMPure cleanup and the final library</h3>\n"
            + panel(rows, cls="long")
            + info(f"{len(con)}&nbsp;bp excluding the insert; no UMI or cell barcode."))


def source_section() -> str:
    read = [
        "Picelli et al. 2013 Online Methods: reverse transcription, PCR preamplification, "
        "tagmentation and enrichment PCR conditions, oligo-dT30VN and ISPCR bases, read length",
        "Picelli et al. 2014 <i>Nature Protocols</i> Reagent Setup: the TSO, oligo-dT30VN and "
        "ISPCR bases of that version",
        "Picelli et al. 2013 Supplementary Tables 2 and 4: every TSO tested, and which one "
        "<i>Smart-seq2</i> uses (<code>rGrG+G</code>)",
        "Picelli et al. 2013 Supplementary Tables 1 and 3: per-library reaction variables",
        "Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1: the Nature "
        "Protocols Smart-seq2 oligo-dT, TSO and forward PCR primer, printed in full",
        "Hahaut et al. 2022 Supplementary Note 2: the same handle oligo-dT and the "
        "biotinylated ISPCR primer, attributed to Picelli et al. 2013",
        "Illumina adapter sequences (via <code>lib/nextera.py</code>) for the Nextera back end",
    ]
    return ('<details class="sources"><summary>Sources read for these bases</summary><ul>'
            + ''.join(f"<li>{item}</li>" for item in read)
            + '</ul></details>'
            + '<details class="sources"><summary>Related protocols</summary><ul>'
            + '<li><a href="https://doi.org/10.1038/nbt.2282">SMART-seq</a></li>'
            + '</ul></details>')


def render_page() -> str:
    graph = protocol_graph()
    preamble = (f'<div class="wrap" data-protocol-complete="true">\n'
                f'<h1>SMART-seq2 &mdash; full-length single-cell mRNA</h1>\n\n'
                + info(f'Defining source: <a href="https://doi.org/{DOI}">doi:{DOI}</a>.')
                + '\n\n\n<details class="background">\n<summary>Background and interpretation</summary>\n'
                + "<p>Template switching installs the 5' handle when MMLV adds "
                  "untemplated Cs at the end of\nthe RNA template. Tn5 then fragments "
                  "the amplified cDNA and adds the sequencing-adapter\nentry points. "
                  "Cell identity is supplied by the well's i5/i7 index pair rather "
                  "than by a\ncell barcode in the molecule.</p>\n</details>\n")
    return assemble_page(graph, [head("SMART-seq2 library chemistry"), preamble,
                                 AnnotationSlot.INTRO, oligo_section(),
                                 WorkflowSection("Step-by-step library generation"),
                                 final_section(graph), AnnotationSlot.AFTER_WORKFLOW,
                                 PrimerBindingSection(
                                     "library", "Library sequencing",
                                     "Sequencing primers: SMART-seq2",
                                     """Read 1 and both index reads use the stock Nextera
            sequencing primers (43-base reads after demultiplexing with i7+i5 indexes, Online Methods; acquired cycle counts unspecified). Primer positions and the first bases read are shown below.""",
                                     info("The i5 + i7 index pair identifies the source well.")
                                     + source_section() + '</div>')])


def selftest() -> None:
    """Check exact source transcription and the inherited-linker boundary."""
    from checks import Check
    check = Check()
    check("Nature Protocols TSO, with ribo and LNA modifications", TSO,
          "AAGCAGTGGTATCAACGCAGAGTACATrGrG+G")
    check("Nature Protocols oligo-dT", OLIGO_DT,
          "AAGCAGTGGTATCAACGCAGAGTAC" + "T" * 30 + "VN")
    check("an inherited TSO linker is inferred by construction",
          preamp_template().get("TSO linker").inferred)
    check("a source-printed ACAT linker can be declared established",
          preamp_template(linker_inferred=False).get("TSO linker").inferred, False)
    check.report()
