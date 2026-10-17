#!/usr/bin/env python3
"""Shared drawing functions for the individual SMART-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import rt
import seqprimers as sp
import smartseq as ss
import smartseq2 as s2
from smartseq_states import (BODY, POLYA, SS2_HANDLE, SS2_TSO, _cdna,
                             _dt_primer, _rt_scenes, gap_filled, mrna as _mrna,
                             ss2_preamp_template, ss2_rt, ss2_tso,
                             ss2_tso_linker, tagmented)
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow, annotation_rows,
                      complement_segments, oligo, panel, protocol_graph_panels,
                      revcomp, workflow_panel)
from illumina import P5, P7
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, AnnotationSlot)
from protocol_views import assemble_page

PRE = 4
H = rt.SMART_HANDLE

# --------------------------------------------------------------- TSO linker provenance
# The SMART handle is followed by a short linker before the G tail, and that linker has
# two published forms:
#
#   AC    Picelli et al. 2013 (the defining SMART-seq2 paper) Supplementary Table 2 prints
#         the selected TSO -- Supplementary Table 4 names `rGrG+G` as the Smart-seq2 TSO --
#         as AAGCAGTGGTATCAACGCAGAGT + AC + rGrG+G. Ramskold et al. 2012 print the CDS
#         oligo-dT with the same AC linker, and Hahaut et al. 2022 Supplementary Note 2
#         prints the FLASH-seq oligo-dT and TSO with it too.
#   ACAT  the Nature Protocols version of Smart-seq2, transcribed verbatim in
#         Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1 ("Smartseq2 /
#         Smart-seq2 as of Nature Protocol"): ...ACATrGrG+G. smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER
#         holds this longer form, which several later SMART-seq2 derivatives inherit.
#
# The page models the longer form but marks the two bases it adds beyond the defining
# paper inferred, so no drawn base is asserted past its source.


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def chunks(*cs):
    return [tuple(list(c) + [None, False])[:3] for c in cs]


# =============================================================================== preamble
PAPERS = {
    "SMART-seq": ("10.1038/nbt.2282", "SMART-seq"),
    "SMART-seq2": (s2.DOI, "SMART-seq2"),
    "SMART-seq3": ("10.1038/s41587-020-0497-0", "SMART-seq3"),
    "SMART-seq3xpress": ("10.1038/s41587-022-01311-4", "SMART-seq3xpress"),
    "FLASH-seq": ("10.1038/s41587-022-01312-3", "FLASH-seq"),
}


CAVEATS = {
    "SMART-seq": ("<div class='caveat'><b>Inferred TSO and ISPCR region.</b> Ramsk&ouml;ld "
                  "et al. print the CDS oligo-dT primer in full but only <i>name</i> the "
                  "Clontech SMARTer II A template-switch oligo and the kit's IS PCR primer; "
                  "neither commercial oligo is published in any source read here. The TSO "
                  "and ISPCR bases below are the later SMART-seq2 / Nature Protocols design "
                  "(as transcribed in Hagemann-Jensen et al. 2020, Supplementary Table 1) "
                  "and are drawn dotted throughout. The paper's &ldquo;standard Nextera "
                  "PCR&rdquo; primer set is not printed either; canonical Nextera index "
                  "primers are shown as representative geometry.</div>"),
}


def preamble(protocol: str, *, complete: bool = False) -> str:
    doi, title = PAPERS[protocol]
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    caveat = CAVEATS.get(protocol, "")
    return f"""<div class="wrap"{complete_attr}>
<h1>{title} &mdash; full-length single-cell mRNA</h1>

{info(f'''Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.''')}
{caveat}

<details class="background">
<summary>Background and interpretation</summary>
<p>Template switching installs the 5' handle when MMLV adds untemplated Cs at the end of
the RNA template. Tn5 then fragments the amplified cDNA and adds the sequencing-adapter
entry points. Cell identity is supplied by the well's i5/i7 index pair rather than by a
cell barcode in the molecule.</p>
</details>
"""


# ====================================================== SMART-seq / SMART-seq2
def ss2_oligos(protocol: str) -> str:
    tail = rt.TSO_G_TAIL_LNA if protocol == "SMART-seq2" else rt.TSO_G_TAIL
    flash = protocol == "FLASH-seq"
    index_kind = ("Representative Nextera" if protocol == "SMART-seq"
                  else "Nextera (XT)")
    rows = [
        oligo("oligo-dTVN", [seg("", H, "tso"), seg("", s2.DT_LINKER),
                             seg("", "T" * s2.DT_LEN, None, placeholder=True),
                             seg("", s2.DT_ANCHOR, None, placeholder=True)],
              mods="/5Biosg/" if flash else ""),
        oligo("Template Switching Oligo (TSO)",
              [seg("", H, "tso", inferred=protocol == "SMART-seq"),
               *[seg("", s.top, s.tag, inferred=s.inferred)
                 for s in ss2_tso_linker(protocol)],
               seg("", tail, None, placeholder=True, inferred=protocol == "SMART-seq")],
              mods="/5Biosg/" if flash else ""),
        oligo("ISPCR primer", [seg("", H, "tso", inferred=protocol == "SMART-seq")],
              mods="/5Biosg/" if flash else ""),
        oligo("Nextera mosaic end (ME)", [seg("", nx.ME, "me")]),
        oligo("Nextera N/S5xx entry point (s5)", [seg("", nx.S5, "s5")]),
        oligo("Nextera N7xx entry point (s7)", [seg("", nx.S7, "s7")]),
        oligo(f"{index_kind} N/S5xx index primer",
              [seg("", P5, "p5"), seg("", "[8-bp i5]", None, placeholder=True), seg("", nx.S5, "s5")]),
        oligo(f"{index_kind} N7xx index primer",
              [seg("", P7, "p7"), seg("", "[8-bp i7]", None, placeholder=True), seg("", nx.S7, "s7")]),
    ]
    return ("<h2>Adapter and primer sequences</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info("""The oligo-dT primer, TSO and ISPCR primer share the same 23-nt
            handle. <code>rG</code> denotes RNA.""" +
                   (" <code>+G</code> denotes the SMART-seq2 LNA." if protocol == "SMART-seq2"
                    else "")))


# ------------------------------------------------ step drawings, placed by pairing
# Every duplex panel below is a chemdraw.Scene: strands are given 5'->3' as ordered or
# synthesised and placed by which segment pairs with which. No column is typed by hand,
# and a Scene that would draw an unpaired column raises instead of rendering.

def _ds_cdna(top: list[Segment], cdna: list[Segment], tso_handle: list[Segment],
             fwd: list[Segment], fwd_pair: tuple[str, str], rev_: list[Segment],
             rev_pair: tuple[str, str]) -> list[Row]:
    """Amplified cDNA: bottom = first strand + its copy of the TSO, top = its complement
    written out, both primers placed by the segment they anneal to."""
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", [*cdna, *complement_segments(tso_handle)], to="top",
              pair=("CCC", "GGG"), label="")
    sc.anneal("fwd", fwd, to="bottom", pair=fwd_pair, label="")
    sc.anneal("rev", rev_, to="top", pair=rev_pair, label="")
    sc.arrow("fwd", "")
    sc.arrow("rev", "")
    sc.stack("fwd", "top", "bottom", "rev")
    return sc.rows()


def ss2_pcr(protocol: str = "SMART-seq2") -> list[Row]:
    top = list(ss2_preamp_template(protocol))
    return _ds_cdna(top, _cdna(SS2_HANDLE, [seg("AC", s2.DT_LINKER)]),
                    ss2_tso(protocol)[:-1],
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H'"),
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H3"))


def ss2_steps(protocol: str = "SMART-seq2") -> str:
    s1, s2, s3 = ss2_rt(protocol)
    tagged, filled = tagmentation_states()
    initial = Scene()
    initial.strand("mRNA", _mrna())
    workflow = Workflow(MolecularState("Starting poly(A) RNA", tuple(initial.rows())))
    workflow.react("oligo-dTVN annealing", s1.rows(),
                   note="The anchored primer sits at the poly(A) junction.")
    workflow.react("MMLV reverse transcription + terminal transferase", s2.rows(),
                   note="Running off the RNA 5′ end adds three untemplated C.")
    workflow.react("template switching", s3.rows(),
                   note=f"The TSO's {'rGrG+G' if protocol == 'SMART-seq2' else 'rGrGrG'} "
                        "pairs with CCC and contributes the second handle.")
    workflow.react("ISPCR", ss2_pcr(protocol),
                   note="The identical terminal handles amplify full-length cDNA.")
    workflow.react("Nextera tagmentation", tagged,
                   note="Tn5 inserts s5/s7 entry points and leaves 9-nt gaps.")
    workflow.react("72 °C gap fill", filled,
                   note="The first Nextera PCR hold fills both gaps before denaturation.")
    return ("<h2>Step-by-step library generation</h2>\n"
            + workflow_panel(workflow, cls="long"))


def tagmentation_states() -> tuple[list[Row], list[Row]]:
    rows = []
    for a, b, amp, why in nx.TAGMENTATION_OUTCOMES:
        rows.append(Row(chunks=chunks((f"Product {a}/{b}: {why}",))))
        t = tagmented(a, b)
        t.same_line("top", "nt-bottom")
        t.same_line("bottom", "nt-top")
        rows.extend(t.rows())
        rows.append(Row())
    return rows, gap_filled().rows()


def ss2_final() -> str:
    con = s2.final_library()
    duplex = Scene.duplex(con.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    rows = duplex.rows() + annotation_rows(con, prefix_width=PRE)
    return ("<h3>(7&ndash;8) Index PCR, and the final library</h3>\n"
            + panel(rows, cls="long")
            + info(f"""{len(con)}&nbsp;bp excluding the insert; no UMI or cell barcode."""))


# ========================================================== SMART-seq3 family
def ss3_tso_segs(variant: str) -> list[Segment]:
    """The TSO drawn from the same constants as smartseq.ss3_tso()."""
    spacer = ss.SS3_TSO_SPACERS[variant]
    return [seg("", ss.SS3_TSO_ME3, "me"), seg("", ss.SS3_TSO_TAG, "r1"),
            seg("", f"[{ss.SS3_UMI_LEN}-bp UMI]", "umi", placeholder=True),
            *([seg("", spacer, None, placeholder=True)] if spacer else []),
            seg("", rt.TSO_G_TAIL, None, placeholder=True)]


def ss3_oligos(variant: str) -> str:
    tso_names = {"SMART-seq3": "Smartseq3_N8_TSO",
                 "SMART-seq3xpress": "Smartseq3xpress_TSO"}
    rows = [
        oligo("Smartseq3_OligodT30VN",
              [seg("", ss.SS3_OLIGO_DT_HANDLE, "r3"),
               seg("", "T" * s2.DT_LEN, None, placeholder=True),
               seg("", s2.DT_ANCHOR, None, placeholder=True)],
              mods="/5Biosg/"),
        oligo(tso_names[variant], ss3_tso_segs(variant), mods="/5Biosg/"),
        oligo("Fwd_PCR_primer",
              [seg("", nx.S5, "s5"), seg("", nx.ME, "me"),
               seg("", ss.SS3_TSO_TAG[:-3], "r1"),
               seg("", "*".join(ss.SS3_TSO_TAG[-3:]), "r1", placeholder=True)]),
        oligo("Rev_PCR_primer", [seg("", ss.SS3_OLIGO_DT_HANDLE[:-3], "r3"),
                                  seg("", "*".join(ss.SS3_OLIGO_DT_HANDLE[-3:]),
                                      "r3", placeholder=True)]),
    ]
    return ("<h2>Adapter and primer sequences</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info(f"""The forward and reverse PCR handles differ. The TSO contains the
            terminal 8 nt of the Nextera ME, the 11-bp 5' tag and the 8-bp UMI.
            <code>*</code> marks the two terminal phosphorothioate linkages on each PCR
            primer. The authors' protocol orders <i>custom</i> Nextera index primers and
            recommends a 10-bp barcode, noting that 8&nbsp;bp &mdash; Illumina's length,
            drawn below &mdash; works as well; no index member is asserted here."""))


SS3_DT_HANDLE = seg("dT handle", ss.SS3_OLIGO_DT_HANDLE, "r3")
SS3_TSO = [seg("ME3", ss.SS3_TSO_ME3, "me"), seg("tag", ss.SS3_TSO_TAG, "r1"),
           seg("UMI", "N" * ss.SS3_UMI_LEN, "umi", placeholder=True),
           seg("rGrG", "GGG", "tso")]                                  # rGrGrG, as GGG


def ss3_tso_scene(variant: str = "SMART-seq3") -> list[Segment]:
    spacer = ss.SS3_TSO_SPACERS[variant]
    return [*SS3_TSO[:3],
            *([seg("spacer", spacer, "r1", placeholder=True)] if spacer else []),
            SS3_TSO[-1]]


def ss3_rt(variant: str = "SMART-seq3") -> tuple[Scene, ...]:
    return _rt_scenes(SS3_DT_HANDLE, [], ss3_tso_scene(variant))


def ss3_pcr(variant: str = "SMART-seq3") -> list[Row]:
    top = list(ss3_preamp_template(variant))
    fwd = [seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"), seg("tag", ss.SS3_TSO_TAG, "r1")]
    return _ds_cdna(top, _cdna(SS3_DT_HANDLE, []), ss3_tso_scene(variant)[:-1],
                    fwd, ("tag", "tag'"),
                    [seg("Rev", ss.SS3_REV_PCR, "r3")], ("Rev", "handle3"))


def ss3_preamp_template(variant: str = "SMART-seq3") -> Construct:
    return Construct([
        *ss3_tso_scene(variant)[:-1], seg("GGG", "GGG", "tso"),
        seg("body", "X" * BODY, placeholder=True),
        seg("B", "B", placeholder=True), seg("polyA", "A" * POLYA),
        seg("handle3", revcomp(ss.SS3_OLIGO_DT_HANDLE), "r3")],
        name=f"{variant} full-length cDNA before PCR")


def ss3_steps(variant: str = "SMART-seq3") -> str:
    s1, s2, s3 = ss3_rt(variant)
    tagged, filled = tagmentation_states()
    initial = Scene()
    initial.strand("mRNA", _mrna())
    workflow = Workflow(MolecularState("Starting poly(A) RNA", tuple(initial.rows())))
    workflow.react("oligo-dT annealing", s1.rows(),
                   note="The anchored primer sits at the poly(A) junction.")
    workflow.react("MMLV reverse transcription + terminal transferase", s2.rows(),
                   note="At the RNA 5′ end, MMLV adds three untemplated C.")
    workflow.react("template switching", s3.rows(),
                   note="The TSO adds the 11-bp tag and 8-bp UMI before amplification.")
    workflow.react("two-primer cDNA PCR", ss3_pcr(variant),
                   note="The forward primer contributes its s5 + ME overhang.")
    workflow.react("Nextera tagmentation", tagged,
                   note="Tn5 inserts s5/s7 entry points and leaves 9-nt gaps.")
    workflow.react("72 °C gap fill", filled,
                   note="The first Nextera PCR hold fills both gaps before denaturation.")
    return ("<h2>Step-by-step library generation</h2>\n"
            + workflow_panel(workflow, cls="long"))


def ss3_final(variant: str = "SMART-seq3") -> str:
    five = ss.smartseq3_library(variant)
    internal = ss.smartseq3_library(five_prime=False)
    five_duplex = Scene.duplex(five.segments)
    five_duplex.strands["top"].label = five_duplex.strands["bottom"].label = ""
    internal_duplex = Scene.duplex(internal.segments)
    internal_duplex.strands["top"].label = internal_duplex.strands["bottom"].label = ""
    r5 = five_duplex.rows() + annotation_rows(five, prefix_width=PRE)
    ri = internal_duplex.rows()
    return ("<h3>(9) Final library structures</h3>\n"
            + "<h3>5' fragments &mdash; these retained the TSO</h3>\n" + panel(r5, cls="long")
            + info("""Contains the <r1>11-bp 5' tag</r1> followed by the
            <umi>8-bp UMI</umi>.""")
            + "<h3>Internal fragments &mdash; these did not</h3>\n" + panel(ri, cls="long"))


# =============================================================== sequencing
def seq_primer(primer: list[Segment], on: str, pair: tuple[str, str]) -> list[Row]:
    """The SMART-seq2 library duplex with one sequencing primer annealed to strand `on`."""
    con = s2.final_library()
    sc = Scene()
    sc.strand("top", list(con), label="")
    sc.anneal("bottom", complement_segments(list(con)), to="top", pair=("cDNA'", "cDNA"),
              label="")
    sc.anneal("primer", primer, to=on, pair=pair, label="")
    sc.arrow("primer", "")
    sc.stack(*(("primer", "top", "bottom") if on == "bottom" else ("top", "bottom", "primer")))
    return sc.rows()


def seq_primer_drawings() -> dict[str, tuple[list[Segment], str, tuple[str, str]]]:
    """Segment drawings of the four stock Nextera primers: (segments, strand, pair).
    Built from the same nextera.py pieces; the self-test asserts each spells exactly the
    sp.NEXTERA primer of the same key."""
    s5, me, s7 = seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"), seg("s7", nx.S7, "s7")
    me_rc = seg("ME_RC", nx.ME_RC, "me")
    return {"R1": ([s5, me], "bottom", ("s5", "s5'")),
            "I1": ([me_rc, seg("s7rc", nx.S7_RC, "s7")], "bottom", ("ME_RC", "ME '")),
            "I2": ([me_rc, seg("s5rc", nx.S5_RC, "s5")], "top", ("ME_RC", "ME")),
            "R2": ([s7, me], "top", ("s7", "s7"))}


def sequencing(protocol: str, *, acquisition_unspecified: bool = False) -> str:
    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        tables = sp.section(s2.final_library(), s2.SEQ_PRIMERS,
                            heading=f"Sequencing primers: {protocol}",
                            show_read_content=True,
                            acquisition_unspecified=acquisition_unspecified)
    else:
        tables = (sp.section(ss.smartseq3_library(protocol), ss.SS3_SEQ_PRIMERS,
                             heading=f"Sequencing primers: {protocol} 5' fragments",
                             show_read_content=True)
                  + sp.section(ss.smartseq3_library(protocol, five_prime=False),
                               ss.SS3_SEQ_PRIMERS,
                               heading=f"Sequencing primers: {protocol} internal fragments",
                               show_read_content=True))
    if protocol == "SMART-seq":
        intro = info("""Primer positions on the canonical Nextera library are shown
            below. Ramsk&ouml;ld et al. Supplementary Table 1 gives only a read length
            (50 on the GAIIx, 40 on the MiSeq for Tn5 libraries) and does not state
            single- or paired-end acquisition or index reads.""")
    else:
        intro = info("""All four reads use the stock Nextera sequencing primers. Primer
            positions and the first bases read are shown below.""")
    return ('<h2 id="seq-primers">Library sequencing</h2>\n'
            + intro
            + tables
            + info("""The i5 + i7 index pair identifies the source well.""")
            + sources(protocol)
            + "</div>")


# Documents actually read for the oligo bases on each page. Every entry is a primary
# publication, author protocol or vendor-grade specification; secondary redrawings of
# these chemistries are not sources and are deliberately absent.
READ_FOR_OLIGOS = {
    "SMART-seq": [
        "Ramsk&ouml;ld et al. 2012 full text (PMC3467340): the CDS oligo-dT primer, verbatim",
        "Picelli et al. 2013 Supplementary Table 2: the SMART-handle TSO series, including "
        "the all-DNA <code>dGCGGG</code> mimic of the SMARTer II A oligo",
        "Illumina adapter sequences (via <code>lib/nextera.py</code>) for ME, s5, s7 and "
        "the N5xx/N7xx index primers",
    ],
    "SMART-seq3": [
        "Hagemann-Jensen et al. 2020 Supplementary Table 1: oligo-dT, TSO, forward and "
        "reverse PCR primers, printed in full",
        "Hagemann-Jensen et al. 2020 author protocol (protocols.io "
        "<code>smart-seq3-protocol-36wgq5rjxgk5</code>): the ordered oligo table with "
        "5' biotin and terminal phosphorothioates",
        "The preprint Methods (bioRxiv 817924v1) restate the same four oligos",
        "Illumina adapter sequences (via <code>lib/nextera.py</code>) for the Nextera back end",
    ],
    "SMART-seq3xpress": [
        "Hagemann-Jensen et al. 2022 Methods (PMC9546772): the improved TSO "
        "<code>5'-Biotin-AGAGACAGATTGCGCAATGNNNNNNNNWWrGrGrG-3'</code> and both PCR primers",
        "Hagemann-Jensen et al. 2022 author protocol (protocols.io "
        "<code>smart-seq3xpress-bwh4pb8w</code>): the ordered oligo table, including the "
        "oligo-dT and the terminal phosphorothioates",
        "Illumina adapter sequences (via <code>lib/nextera.py</code>) for the Nextera back end",
    ],
    "FLASH-seq": [
        "Hahaut et al. 2022 Supplementary Note 2, &ldquo;List of all the oligonucleotides "
        "used in this study&rdquo;: FS dT30VN, FS TSO, ISPCR and the FS-UMI TSO series",
        "Hahaut et al. 2022 Methods (PMC9546769): the same oligos inside the reaction recipes",
        "Illumina adapter sequences (via <code>lib/nextera.py</code>) for the Nextera back end",
    ],
}


def sources(protocol: str) -> str:
    related = {
        "SMART-seq": [],
        "SMART-seq2": [("SMART-seq", PAPERS["SMART-seq"][0])],
        "SMART-seq3": [("SMART-seq2", PAPERS["SMART-seq2"][0])],
        "SMART-seq3xpress": [("SMART-seq3", PAPERS["SMART-seq3"][0])],
        "FLASH-seq": [("SMART-seq2", PAPERS["SMART-seq2"][0])],
    }[protocol]
    read = ''.join(f"<li>{item}</li>" for item in READ_FOR_OLIGOS[protocol])
    links = ''.join(f'<li><a href="https://doi.org/{doi}">{name}</a></li>'
                    for name, doi in related)
    out = ('<details class="sources"><summary>Sources read for these bases</summary><ul>'
           + read + '</ul></details>')
    if links:
        out += ('<details class="sources"><summary>Related protocols</summary><ul>'
                + links + '</ul></details>')
    return out


def protocol_graph() -> ProtocolGraph:
    """The original SMART-seq graph, Tn5 ("Tn5 protocol") branch only.

    Ramskold et al. 2012 (PMC3467340) Methods, "Generation and amplification of
    Smart-Seq cDNA" and "Construction and sequencing of Smart-Seq sequencing
    libraries": 12/15/18 preamplification cycles; ~5 ng cDNA tagmented at 55 C for
    5 min (0.25 ul transposase, 4 ul 5x HMW Nextera buffer, 20 ul) in a modified
    Epicentre Nextera DNA sample preparation protocol; 35 ul PB strips the
    transposase; SPRI XP 1:1.6; 9 cycles of standard Nextera PCR. Commercial
    TSO/ISPCR bases and the Nextera PCR primer set are not printed there.
    """
    source = "Defining source: doi:10.1038/nbt.2282"
    graph = ProtocolGraph("SMART-seq", authority=source)
    mrna = input_mrna(graph)
    first, enzyme, switch_scene, template = reverse_transcribe(
        graph, "SMART-seq", mrna, rt_label="anchored oligo-dTVN with SMART handle")
    switched = template_switch(
        graph, first, enzyme, switch_scene, template,
        tso_label="SMART-seq template-switch oligo",
        step_label=("INFERRED — Template-switch at the transcript 5-prime end; the "
                    "commercial SMARTer II A oligo's bases are not published"))
    amplified = preamplify(
        graph, switched, template, family="ss2", inferred=True,
        label=("PCR-amplify full-length cDNA: 12 cycles (1 ng total RNA), 15 "
               "(100 pg) or 18 (10 pg or single cell)"),
        statement=("INFERRED — The original SMARTer II A oligo and the kit's IS PCR "
                   "primer are commercial reagents whose complete sequences no source "
                   "read here prints; the later SMART-seq2 handle, printed in Picelli "
                   "et al. 2013 Supplementary Table 2 and in Hahaut et al. 2022 "
                   "Supplementary Note 2, gives this representative symmetric geometry."))
    tagged = tagment(graph, amplified, label=(
        "Tagment ~5 ng cDNA at 55 degrees C for 5 min (0.25 ul transposase, 5x HMW "
        "Nextera buffer; modified Epicentre Nextera protocol)"))
    beads = graph.input("tn5_strip_beads", "PB buffer, then SPRI XP beads 1:1.6",
                        MaterialKind.SELECTION)
    cleaned = graph.product("tagmented_purified", "purified s5/s7-tagmented cDNA fragments",
                            MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    cleanup_waste = graph.product("tn5_strip_waste",
                                  "stripped transposase and bead supernatant",
                                  MaterialKind.WASTE)
    graph.react("tn5_strip_cleanup",
                "Strip transposase with 35 ul PB; purify with 88 ul SPRI XP beads (1:1.6)",
                StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (cleaned, cleanup_waste))
    graph.discard(cleanup_waste, "removed with the bead supernatant and washes")
    filled = fill_tn5_gaps(graph, cleaned, label=(
        "Fill both nine-base Tn5 gaps (first extension of the standard Nextera PCR; "
        "hold conditions not printed in the defining paper)"))
    index_pcr(graph, filled, s2.final_library(), tuple(s2.SEQ_PRIMERS), name="SMART-seq",
              label="9 cycles of standard Nextera PCR install P5/P7 ends",
              primer_labels=("Nextera PCR primer, P5 side (2012 primer set not printed; "
                             "canonical N5xx geometry drawn)",
                             "Nextera PCR primer, P7 side (2012 primer set not printed; "
                             "canonical N7xx geometry drawn)"),
              acquisition_unspecified=True)
    graph.annotate(Highlight(
        "tagment", "Library branch drawn, and what the source leaves open",
        "Ramsk&ouml;ld et al. made libraries two ways: the Illumina &ldquo;PE&rdquo; "
        "protocol (Covaris shearing, end repair, A-tailing, PE adaptor ligation, 12&ndash;18 "
        "PCR cycles) and the Tn5 protocol drawn here. Their Methods state that all figure "
        "data used the PE protocol unless a legend says otherwise. Supplementary Table 1 "
        "lists 41 Tn5 libraries read at 50 on the GAIIx and 2 at 40 on the MiSeq; it does "
        "not state the index primers used, so the 8-base index windows below are the "
        "canonical Nextera geometry. Which reads were acquired is marked as unstated in "
        "the read table."))
    return review_smart_graph(graph, source, family="ss2", extra=("tn5_strip_beads", "tn5_strip_cleanup",
                                                                 "tagmented_purified"))


def input_mrna(graph: ProtocolGraph):
    initial = Scene(); initial.strand("mRNA", _mrna())
    return graph.input("mrna", "single-cell poly(A) RNA", MaterialKind.MOLECULAR,
                       state=initial, primary=True)


def reverse_transcribe(graph: ProtocolGraph, protocol: str, mrna, *, rt_label: str):
    if protocol in ("SMART-seq", "FLASH-seq"):
        priming, first_scene, switch_scene = ss2_rt(protocol)
        template = ss2_preamp_template(protocol)
    else:
        priming, first_scene, switch_scene = ss3_rt(protocol)
        template = ss3_preamp_template(protocol)
    primer = graph.input("rt_primer", rt_label, MaterialKind.OLIGO, state=priming)
    enzyme = graph.input("reverse_transcriptase",
                         "MMLV-derived template-switching reverse transcriptase",
                         MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with untemplated CCC",
                          MaterialKind.MOLECULAR, state=first_scene)
    graph.react("reverse_transcription", "Reverse-transcribe from the anchored poly(dT) primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (first,))
    return first, enzyme, switch_scene, template


def template_switch(graph: ProtocolGraph, first, enzyme, switch_scene, template,
                    *, tso_label: str, step_label: str):
    tso = graph.input("tso", tso_label, MaterialKind.OLIGO, state=switch_scene)
    switched = graph.product("switched", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR, state=template)
    graph.react("template_switch", step_label, StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (switched,))
    return switched


def preamplify(graph: ProtocolGraph, switched, template, *, family: str,
               statement: str, inferred: bool = False, label: str = "",
               polymerase_label: str = "full-length cDNA PCR polymerase"):
    polymerase = graph.input("preamp_polymerase", polymerase_label,
                             MaterialKind.CATALYST)
    if family == "ss2":
        ispcr = graph.input("ispcr", "ISPCR primer used at both cDNA ends", MaterialKind.OLIGO)
        geometry = PCRGeometry.symmetric(template, PCRPrimer("ISPCR", H), statement=statement)
        inputs = ((switched, InputRole.TEMPLATE), (ispcr, InputRole.STRUCTURAL),
                  (polymerase, InputRole.CATALYST))
    elif family == "ss3":
        forward = graph.input("preamp_forward", "s5-ME-tag forward PCR primer", MaterialKind.OLIGO)
        reverse = graph.input("preamp_reverse", "oligo-dT-handle reverse PCR primer", MaterialKind.OLIGO)
        geometry = PCRGeometry.exact(template, PCRPrimer("Fwd PCR primer", ss.SS3_FWD_PCR),
                                     PCRPrimer("Rev PCR primer", ss.SS3_REV_PCR), statement=statement)
        inputs = ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                  (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST))
    else:
        raise ValueError("unknown SMART-seq preamplification family")
    amplified = graph.product("amplified", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=template)
    label = label or "PCR-amplify full-length cDNA"
    graph.react("preamp", ("INFERRED — " + label if inferred else label), StepKind.PCR, inputs,
                (amplified,), geometry=geometry)
    return amplified


def tagment(graph: ProtocolGraph, amplified, *,
            label: str = "Tagment amplified cDNA with Nextera s5/s7 adapters"):
    transposome = graph.input("transposome", "s5/s7-loaded Nextera transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented", "heterologous s5/s7-tagmented cDNA fragments",
                           MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    waste = graph.product("tagment_waste", "s5/s5, s7/s7 and other non-library Tn5 products",
                          MaterialKind.WASTE)
    graph.react("tagment", label,
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "lacks the heterologous s5/s7 end combination required for amplification")
    return tagged


def fill_tn5_gaps(graph: ProtocolGraph, tagged, *,
                  label: str = "Fill both nine-base Tn5 gaps at 72 degrees C"):
    polymerase = graph.input("gap_polymerase", "Nextera nine-base gap-fill polymerase",
                             MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterologous Nextera fragments",
                           MaterialKind.MOLECULAR, state=gap_filled())
    graph.react("gap_fill", label, StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,))
    return filled


def index_pcr(graph: ProtocolGraph, filled, final: Construct, primers: tuple,
              *, name: str, internal: Construct | None = None,
              label: str = "Amplify and install indexed P5/P7 ends",
              primer_labels: tuple[str, str] = ("Nextera N5xx index primer",
                                                "Nextera N7xx index primer"),
              acquisition_unspecified: bool = False):
    forward = graph.input("index_forward", primer_labels[0], MaterialKind.OLIGO)
    reverse = graph.input("index_reverse", primer_labels[1], MaterialKind.OLIGO)
    polymerase = graph.input("index_polymerase", "Nextera index-PCR polymerase", MaterialKind.CATALYST)
    geometry = PCRGeometry.exact(
        final, PCRPrimer("N5xx", nx.n5xx_primer("N" * 8)),
        PCRPrimer("N7xx", nx.n7xx_primer("N" * 8)),
        statement="The canonical N5xx and N7xx primers face inward from the s5 and s7 entry points and install the indexed P5/P7 ends.")
    library = graph.product("library", name + " sequencing-library molecules",
                            MaterialKind.LIBRARY, state=final)
    outputs = [library]
    if internal is not None:
        outputs.append(graph.product("internal_library", name + " internal-fragment library molecules",
                                     MaterialKind.LIBRARY, state=internal))
    graph.react("index_pcr", label, StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                outputs, geometry=geometry)
    for material, construct in zip(outputs, (final, internal)):
        graph.sequence(material, SequencingPlan(
            library=construct, primers=primers,
            acquired_roles=(() if acquisition_unspecified
                            else tuple(primer.role for primer in primers)),
            acquisition="unspecified" if acquisition_unspecified else "declared"))


def review_smart_graph(graph: ProtocolGraph, source: str, *, family: str,
                       internal: bool = False, extra: tuple[str, ...] = ()):
    reviewed = ["mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
                "tso", "template_switch", "preamp_polymerase", "preamp", "transposome",
                "tn5", "tagment", "gap_polymerase", "gap_fill", "index_forward",
                "index_reverse", "index_polymerase", "index_pcr", "library"]
    reviewed.extend(("ispcr",) if family == "ss2" else ("preamp_forward", "preamp_reverse"))
    if internal:
        reviewed.append("internal_library")
    reviewed.extend(extra)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render_page(protocol: str, *, graph_factory=None) -> str:
    if protocol == "SMART-seq2":
        raise ValueError("SMART-seq2 page is rendered by its own protocol module")
    if graph_factory is None and protocol != "SMART-seq":
        raise ValueError("derived SMART pages must supply their protocol-owned graph")
    graph = graph_factory() if graph_factory else protocol_graph()
    graph_workflow = ("<h2>Step-by-step library generation</h2>\n"
                      + protocol_graph_panels(graph, cls="long"))
    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        body = [ss2_oligos(protocol), graph_workflow, ss2_final()]
    else:
        body = [ss3_oligos(protocol), graph_workflow, ss3_final(protocol)]
    return assemble_page(graph, [
        head(f"{protocol} library chemistry"), preamble(protocol, complete=True),
        AnnotationSlot.INTRO, *body, AnnotationSlot.AFTER_WORKFLOW,
        sequencing(protocol, acquisition_unspecified=any(
            plan.acquisition_unspecified for plan in graph.sequencing.values()))])
