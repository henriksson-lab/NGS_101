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
from protocol_graph import (AnnotationSlot, Highlight, InputRole, MaterialKind, MolecularRelation,
                            ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
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
    rt_label = "anchored oligo-dTVN with SMART handle"
    tso_label = "SMART-seq2 LNA template-switch oligo"
    preamp_template = drawn_preamp_template(BODY, POLYA)

    rt_primer = graph.input("rt_primer", rt_label, MaterialKind.OLIGO,
                            state=priming)
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "MMLV-derived template-switching reverse transcriptase",
        MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with untemplated CCC",
                          MaterialKind.MOLECULAR, state=first_scene)
    graph.react("reverse_transcription",
                "Reverse-transcribe from the anchored poly(dT) primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))

    tso = graph.input("tso", tso_label, MaterialKind.OLIGO, state=switch_scene)
    switched = graph.product("switched", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR, state=preamp_template)
    switch_label = TEMPLATE_SWITCH_CAPTION
    graph.react("template_switch", switch_label,
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))

    preamp_polymerase = graph.input("preamp_polymerase", "full-length cDNA PCR polymerase",
                                    MaterialKind.CATALYST)
    ispcr = graph.input("ispcr", "ISPCR primer used at both cDNA ends",
                        MaterialKind.OLIGO)
    preamp_geometry = PCRGeometry.symmetric(
        preamp_template, PCRPrimer("ISPCR", rt.SMART_HANDLE),
        statement=("The source-defined SMART handle occurs in opposite orientations "
                   "at the two cDNA ends, so one ISPCR reagent supplies both inward sites."))
    preamp_inputs = ((switched, InputRole.TEMPLATE),
                     (ispcr, InputRole.STRUCTURAL),
                     (preamp_polymerase, InputRole.CATALYST))
    primer_reviewed = ["ispcr"]
    amplified = graph.product("amplified", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=preamp_template)
    graph.react("preamp", "PCR-amplify full-length cDNA", StepKind.PCR,
                preamp_inputs, (amplified,), geometry=preamp_geometry)

    index_f_primer = PCRPrimer("N5xx", nx.n5xx_primer("N" * 8), anneal_nt=len(nx.S5))
    index_r_primer = PCRPrimer("N7xx", nx.n7xx_primer("N" * 8), anneal_nt=len(nx.S7))
    transposome = graph.input("transposome", "s5/s7-loaded Nextera transposome",
                              MaterialKind.ADAPTER, chemical_state=LoadedEndPool(
                                  (LoadedEnd("s5", nx.ADAPTOR_S5),
                                   LoadedEnd("s7", nx.ADAPTOR_S7)), nx.ME))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented_product = graph.product(
        "tagmented", "heterologous s5/s7-tagmented cDNA fragments",
        MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    tagment_waste = graph.product(
        "tagment_waste", "s5/s5, s7/s7 and other non-library Tn5 products",
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
    graph.discard(tagment_waste,
                  "lacks the heterologous s5/s7 end combination required for amplification")

    gap_polymerase = graph.input("gap_polymerase", "Nextera nine-base gap-fill polymerase",
                                 MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterologous Nextera fragments",
                           MaterialKind.MOLECULAR, state=gap_filled())
    graph.react("gap_fill", "Fill both nine-base Tn5 gaps at 72 degrees C",
                StepKind.EXTENSION,
                ((tagmented_product, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))

    final = final_library()
    index_forward = graph.input("index_forward", "Nextera N5xx index primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse", "Nextera N7xx index primer",
                                MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "Nextera index-PCR polymerase",
                                   MaterialKind.CATALYST)
    index_geometry = PCRGeometry.exact(
        final, index_f_primer, index_r_primer,
        statement=("The canonical N5xx and N7xx primers face inward from the s5 and s7 "
                   "entry points and install the indexed P5/P7 ends."))
    primers = SEQ_PRIMERS
    main_library = graph.product("library", protocol + " sequencing-library molecules",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Amplify and install indexed P5/P7 ends",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (index_polymerase, InputRole.CATALYST)),
                (main_library,), geometry=index_geometry)
    graph.sequence(main_library, SequencingPlan(
        library=final, primers=primers,
        acquired_roles=tuple(p.role for p in primers)))

    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "tso", "template_switch", *primer_reviewed, "preamp_polymerase", "preamp",
        "transposome", "tn5", "tagment", "gap_polymerase", "gap_fill",
        "index_forward", "index_reverse", "index_polymerase", "index_pcr",
        "library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "template_switch", "Two published TSO linkers",
        (SourceClaim("Picelli et al. 2013, Supplementary Tables 2 and 4",
                     "The selected TSO is the 23-nt handle + AC + rGrG+G."),
         SourceClaim("Nature Protocols form, transcribed in Hagemann-Jensen "
                     "et al. 2020, Supplementary Table 1, condition 1",
                     "The Smart-seq2 TSO has ACAT before rGrG+G.")),
        "The drawn TSO uses ACAT.",
        "The later protocol's oligo table explicitly prints this longer form.",
        "The added AT is absent from the defining paper and is drawn inferred."))
    graph.annotate(Highlight(
        "tagment", "Selection after Tn5",
        "Heterologous s5/s7 fragments have both index-PCR sites. Same-end s5/s5 "
        "and s7/s7 products are represented as discarded products and checked "
        "against the declared primer pair."))
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
        OligoEntry("Nextera (XT) N/S5xx index primer",
                   (seg(P5, "p5"), seg("[8-bp i5]", placeholder=True),
                    seg(nx.S5, "s5"))),
        OligoEntry("Nextera (XT) N7xx index primer",
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
    return ("<h3>(7&ndash;8) Index PCR, and the final library</h3>\n"
            + panel(rows, cls="long")
            + info(f"{len(con)}&nbsp;bp excluding the insert; no UMI or cell barcode."))


def source_section() -> str:
    read = [
        "Picelli et al. 2013 Supplementary Tables 2 and 4: every TSO tested, and which one "
        "<i>Smart-seq2</i> uses (<code>rGrG+G</code>)",
        "Picelli et al. 2013 Supplementary Tables 1 and 3: the oligo-dT named as the "
        "commercial <i>SMARTer dT30VN</i>, and the ISPCR primer named but not printed",
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
                                     """All four reads use the stock Nextera sequencing primers. Primer
            positions and the first bases read are shown below.""",
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
