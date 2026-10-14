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
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow, annotation_rows,
                      complement_segments, oligo, panel, protocol_graph_panels,
                      revcomp, workflow_panel)
from illumina import P5, P7
from page import head, info, table
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

PRE = 4
H = rt.SMART_HANDLE


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def chunks(*cs):
    return [tuple(list(c) + [None, False])[:3] for c in cs]


# =============================================================================== preamble
PAPERS = {
    "SMART-seq": ("10.1038/nbt.2282", "SMART-seq"),
    "SMART-seq2": ("10.1038/nmeth.2639", "SMART-seq2"),
    "SMART-seq3": ("10.1038/s41587-020-0497-0", "SMART-seq3"),
    "SMART-seq3xpress": ("10.1038/s41587-022-01311-4", "SMART-seq3xpress"),
    "FLASH-seq": ("10.1038/s41587-022-01312-3", "FLASH-seq"),
}


def preamble(protocol: str, *, complete: bool = False) -> str:
    doi, title = PAPERS[protocol]
    complete_attr = ' data-protocol-complete="true"' if complete else ""
    caveat = ("<div class='caveat'><b>Inferred TSO region.</b> The original SMART-seq "
              "paper names the SMARTer II A template-switch oligo but does not print its "
              "bases; the TSO linker below follows the later SMART-seq2 design.</div>"
              if protocol == "SMART-seq" else "")
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
    rows = [
        oligo("oligo-dTVN", [seg("", H, "tso"), seg("", ss.SS2_DT_LINKER),
                             seg("", "T" * ss.DT_LEN, None, placeholder=True),
                             seg("", ss.DT_ANCHOR, None, placeholder=True)],
              mods="/5Biosg/" if flash else ""),
        oligo("Template Switching Oligo (TSO)",
              [seg("", H, "tso", inferred=protocol == "SMART-seq"),
               seg("", ss.SS2_DT_LINKER if flash else ss.SS2_TSO_LINKER,
                   inferred=protocol == "SMART-seq"),
               seg("", tail, None, placeholder=True, inferred=protocol == "SMART-seq")],
              mods="/5Biosg/" if flash else ""),
        oligo("ISPCR primer", [seg("", H, "tso", inferred=protocol == "SMART-seq")],
              mods="/5Biosg/" if flash else ""),
        oligo("Nextera mosaic end (ME)", [seg("", nx.ME, "me")]),
        oligo("Nextera N/S5xx entry point (s5)", [seg("", nx.S5, "s5")]),
        oligo("Nextera N7xx entry point (s7)", [seg("", nx.S7, "s7")]),
        oligo("Nextera (XT) N/S5xx index primer",
              [seg("", P5, "p5"), seg("", "[8-bp i5]", None, placeholder=True), seg("", nx.S5, "s5")]),
        oligo("Nextera (XT) N7xx index primer",
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

BODY = 20           # drawn transcript body
POLYA = 30          # drawn poly(A): exactly the oligo-dT's (T)30, so the handle overhangs


def _mrna() -> list[Segment]:
    return [seg("body", "X" * BODY, placeholder=True),
            seg("B", "B", placeholder=True),
            seg("polyA", "A" * POLYA)]


def _dt_primer(handle: Segment, link: list[Segment]) -> list[Segment]:
    """Anchored oligo-dT, 5'->3': handle, (linker), T30, V, N."""
    return [handle, *link, seg("dT", "T" * POLYA),
            seg("V", "V", placeholder=True), seg("N", "N", placeholder=True)]


def _cdna(handle: Segment, link: list[Segment]) -> list[Segment]:
    """First strand 5'->3': the oligo-dT, the copied body (the N already paired the last
    body base), then MMLV's untemplated CCC."""
    return [*_dt_primer(handle, link),
            *complement_segments([seg("body", "X" * (BODY - 1), placeholder=True)]),
            seg("CCC", rt.UNTEMPLATED_TAIL, "tso")]


def _rt_scenes(handle: Segment, link: list[Segment], tso: list[Segment]) -> tuple[Scene, ...]:
    """(priming, first strand + CCC, template switch) for one chemistry."""
    s1 = Scene()
    s1.strand("mRNA", _mrna())
    s1.anneal("oligo-dT", _dt_primer(handle, link), to="mRNA", pair=("dT", "polyA"))
    s1.mark("oligo-dT", "N", "VN anchor: sits at the poly(A) junction", through="V")
    s1.arrow("oligo-dT", "MMLV")

    s2 = Scene()
    s2.strand("mRNA", _mrna())
    s2.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s2.mark("cDNA", "CCC", "untemplated CCC, overhanging the mRNA 5' end")

    s3 = Scene()
    s3.strand("mRNA", _mrna())
    s3.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s3.anneal("TSO", tso, to="cDNA", pair=("rGrG", "CCC"), above=True)
    s3.mark("cDNA", "CCC", "CCC pairs the TSO's G tail")
    s3.arrow("cDNA", "RT switches template")
    return s1, s2, s3


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


SS2_HANDLE = seg("H", H, "tso")
SS2_TSO = [SS2_HANDLE, seg("ACAT", ss.SS2_TSO_LINKER), seg("rGrG", "GGG", "tso")]   # rGrG+G, as GGG


def ss2_rt(protocol: str = "SMART-seq2") -> tuple[Scene, ...]:
    tso = ([SS2_HANDLE, seg("AC", ss.SS2_DT_LINKER), seg("rGrG", "GGG", "tso")]
           if protocol == "FLASH-seq" else SS2_TSO)
    if protocol == "SMART-seq":
        tso = [seg("H", H, "tso", inferred=True),
               seg("ACAT", ss.SS2_TSO_LINKER, inferred=True),
               seg("rGrG", "GGG", "tso", inferred=True)]
    return _rt_scenes(SS2_HANDLE, [seg("AC", ss.SS2_DT_LINKER)], tso)


def ss2_pcr() -> list[Row]:
    top = list(ss2_preamp_template())
    return _ds_cdna(top, _cdna(SS2_HANDLE, [seg("AC", ss.SS2_DT_LINKER)]), SS2_TSO[:2],
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H'"),
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H3"))


def ss2_preamp_template(protocol: str = "SMART-seq2") -> Construct:
    tso_linker = ss.SS2_DT_LINKER if protocol == "FLASH-seq" else ss.SS2_TSO_LINKER
    return Construct([
        seg("H5", H, "tso", inferred=protocol == "SMART-seq"),
        seg("TSO linker", tso_linker, inferred=protocol == "SMART-seq"),
        seg("GGG", "GGG", "tso", inferred=protocol == "SMART-seq"),
        seg("body", "X" * BODY, placeholder=True),
        seg("B", "B", placeholder=True), seg("polyA", "A" * POLYA),
        seg("GT", revcomp(ss.SS2_DT_LINKER)),
        seg("H3", revcomp(H), "tso")],
        name="SMART-seq full-length cDNA before ISPCR")


def ss2_steps(protocol: str = "SMART-seq2") -> str:
    s1, s2, s3 = ss2_rt()
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
    workflow.react("ISPCR", ss2_pcr(),
                   note="The identical terminal handles amplify full-length cDNA.")
    workflow.react("Nextera tagmentation", tagged,
                   note="Tn5 inserts s5/s7 entry points and leaves 9-nt gaps.")
    workflow.react("72 °C gap fill", filled,
                   note="The first Nextera PCR hold fills both gaps before denaturation.")
    return ("<h2>Step-by-step library generation</h2>\n"
            + workflow_panel(workflow, cls="long"))


def _entry(which: str) -> Segment:
    return seg(which, nx.S5 if which == "s5" else nx.S7, which)


def tagmented(a: str, b: str) -> Scene:
    """One Tn5 product, before gap fill.

    Each fragment strand carries a transferred adaptor (entry + ME) on its 5' end. The
    non-transferred strand (ME_RC) stays annealed to each ME, and Tn5's 9-bp stagger leaves
    each fragment 3' end 9 nt short of it: a 9-nt gap opposite the other strand's first 9
    insert bases.
    """
    g = nx.TAGMENTATION_GAP
    me = seg("ME", nx.ME, "me")
    top = [_entry(a), me, seg("ss", "X" * g, placeholder=True),
           seg("core", "XXXXXX...XXXXXX", placeholder=True)]
    bot = [_entry(b), me, seg("ss", "x" * g, placeholder=True),
           *complement_segments([top[-1]])]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", bot, to="top", pair=("core'", "core"), label="")
    sc.anneal("nt-top", complement_segments([me], suffix=""), to="top", pair=("ME", "ME"),
              label="")
    sc.anneal("nt-bottom", complement_segments([me], suffix=""), to="bottom",
              pair=("ME", "ME"), label="", above=True)
    return sc


def gap_filled() -> Scene:
    g = nx.TAGMENTATION_GAP
    top = [_entry("s5"), seg("ME", nx.ME, "me"), seg("ss", "X" * g, placeholder=True),
           seg("core", "XXXXXX...XXXXXX", placeholder=True),
           seg("fill", "X" * g, placeholder=True),
           seg("ME_RC", nx.ME_RC, "me"), seg("s7rc", nx.S7_RC, "s7")]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", complement_segments(top), to="top", pair=("core'", "core"), label="")
    sc.mark("top", "fill", "top 3' end: gap filled, ME' displaced, copied to the end",
            through="s7rc")
    sc.mark("bottom", "ss'", "bottom 3' end: likewise", through="s5'")
    return sc


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
    con = ss.smartseq2_library()
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
               seg("", "T" * ss.DT_LEN, None, placeholder=True),
               seg("", ss.DT_ANCHOR, None, placeholder=True)],
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
            <code>*</code> marks the two terminal phosphorothioate linkages on each PCR primer."""))


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
    con = ss.smartseq2_library()
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


def sequencing(protocol: str) -> str:
    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        tables = sp.section(ss.smartseq2_library(), ss.SS2_SEQ_PRIMERS,
                            heading=f"Sequencing primers: {protocol}",
                            show_read_content=True)
    else:
        tables = (sp.section(ss.smartseq3_library(protocol), ss.SS3_SEQ_PRIMERS,
                             heading=f"Sequencing primers: {protocol} 5' fragments",
                             show_read_content=True)
                  + sp.section(ss.smartseq3_library(protocol, five_prime=False),
                               ss.SS3_SEQ_PRIMERS,
                               heading=f"Sequencing primers: {protocol} internal fragments",
                               show_read_content=True))
    return ('<h2 id="seq-primers">Library sequencing</h2>\n'
            + info("""All four reads use the stock Nextera sequencing primers. Primer
            positions and the first bases read are shown below.""")
            + tables
            + info("""The i5 + i7 index pair identifies the source well.""")
            + sources(protocol)
            + "</div>")


def sources(protocol: str) -> str:
    related = {
        "SMART-seq": [],
        "SMART-seq2": [("SMART-seq", PAPERS["SMART-seq"][0])],
        "SMART-seq3": [("SMART-seq2", PAPERS["SMART-seq2"][0])],
        "SMART-seq3xpress": [("SMART-seq3", PAPERS["SMART-seq3"][0])],
        "FLASH-seq": [("SMART-seq2", PAPERS["SMART-seq2"][0])],
    }[protocol]
    links = ''.join(f'<li><a href="https://doi.org/{doi}">{name}</a></li>'
                    for name, doi in related)
    links += ('<li><a href="https://teichlab.github.io/scg_lib_structs/methods_html/'
              'SMART-seq_family.html">scg_lib_structs schematic</a></li>')
    return ('<details class="sources"><summary>Related sources</summary><ul>'
            + links + '</ul></details>')


def protocol_graph(protocol: str) -> ProtocolGraph:
    doi, _ = PAPERS[protocol]
    source = f"Defining source: doi:{doi}"
    graph = ProtocolGraph(protocol, authority=source)
    initial = Scene(); initial.strand("mRNA", _mrna())
    mrna = graph.input("mrna", "single-cell poly(A) RNA", MaterialKind.MOLECULAR,
                       state=initial, primary=True)

    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        priming, first_scene, switch_scene = ss2_rt(protocol)
        rt_label = ("biotinylated anchored oligo-dTVN with SMART handle"
                    if protocol == "FLASH-seq" else
                    "anchored oligo-dTVN with SMART handle")
        tso_label = ("SMART-seq2 LNA template-switch oligo" if protocol == "SMART-seq2"
                     else protocol + " template-switch oligo")
        preamp_template = ss2_preamp_template(protocol)
    else:
        priming, first_scene, switch_scene = ss3_rt(protocol)
        rt_label = "biotinylated anchored oligo(dT) primer"
        tso_label = protocol + " UMI template-switch oligo"
        preamp_template = ss3_preamp_template(protocol)

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
    graph.react("template_switch", ("INFERRED — Template-switch at the transcript 5-prime end"
                                    if protocol == "SMART-seq" else
                                    "Template-switch at the transcript 5-prime end"),
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))

    preamp_polymerase = graph.input("preamp_polymerase", "full-length cDNA PCR polymerase",
                                    MaterialKind.CATALYST)
    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        ispcr = graph.input("ispcr", "ISPCR primer used at both cDNA ends",
                            MaterialKind.OLIGO)
        preamp_geometry = PCRGeometry.symmetric(
            preamp_template, PCRPrimer("ISPCR", H),
            statement=(("INFERRED — The original SMARTer II A oligo and cDNA PCR primer "
                        "are kit reagents without printed complete sequences; the later "
                        "SMART-seq2 handle gives this representative symmetric geometry.")
                       if protocol == "SMART-seq" else
                       ("The source-defined SMART handle occurs in opposite orientations "
                        "at the two cDNA ends, so one ISPCR reagent supplies both inward sites.")))
        preamp_inputs = ((switched, InputRole.TEMPLATE),
                         (ispcr, InputRole.STRUCTURAL),
                         (preamp_polymerase, InputRole.CATALYST))
        primer_reviewed = ["ispcr"]
    else:
        forward = graph.input("preamp_forward", "s5-ME-tag forward PCR primer",
                              MaterialKind.OLIGO)
        reverse = graph.input("preamp_reverse", "oligo-dT-handle reverse PCR primer",
                              MaterialKind.OLIGO)
        preamp_geometry = PCRGeometry.exact(
            preamp_template,
            PCRPrimer("Fwd PCR primer", ss.SS3_FWD_PCR),
            PCRPrimer("Rev PCR primer", ss.SS3_REV_PCR),
            statement=("The published asymmetric forward and reverse primers face inward "
                       "around the full-length template-switched cDNA."))
        preamp_inputs = ((switched, InputRole.TEMPLATE),
                         (forward, InputRole.STRUCTURAL),
                         (reverse, InputRole.STRUCTURAL),
                         (preamp_polymerase, InputRole.CATALYST))
        primer_reviewed = ["preamp_forward", "preamp_reverse"]
    amplified = graph.product("amplified", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=preamp_template)
    graph.react("preamp", ("INFERRED — PCR-amplify full-length cDNA" if protocol == "SMART-seq"
                            else "PCR-amplify full-length cDNA"), StepKind.PCR,
                preamp_inputs, (amplified,), geometry=preamp_geometry)

    transposome = graph.input("transposome", "s5/s7-loaded Nextera transposome",
                              MaterialKind.ADAPTER)
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
                (tagmented_product, tagment_waste))
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

    final = (ss.smartseq2_library() if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq")
             else ss.smartseq3_library(protocol))
    index_forward = graph.input("index_forward", "Nextera N5xx index primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse", "Nextera N7xx index primer",
                                MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "Nextera index-PCR polymerase",
                                   MaterialKind.CATALYST)
    index_geometry = PCRGeometry.exact(
        final, PCRPrimer("N5xx", nx.n5xx_primer("N" * 8)),
        PCRPrimer("N7xx", nx.n7xx_primer("N" * 8)),
        statement=("The canonical N5xx and N7xx primers face inward from the s5 and s7 "
                   "entry points and install the indexed P5/P7 ends."))
    primers = tuple(ss.SS2_SEQ_PRIMERS if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq")
                    else ss.SS3_SEQ_PRIMERS)
    finals = []
    main_library = graph.product("library", protocol + " sequencing-library molecules",
                                 MaterialKind.LIBRARY, state=final)
    finals.append((main_library, final))
    outputs = [main_library]
    if protocol not in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        internal = ss.smartseq3_library(protocol, five_prime=False)
        internal_library = graph.product(
            "internal_library", protocol + " internal-fragment library molecules",
            MaterialKind.LIBRARY, state=internal)
        finals.append((internal_library, internal)); outputs.append(internal_library)
    graph.react("index_pcr", "Amplify and install indexed P5/P7 ends",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (index_polymerase, InputRole.CATALYST)),
                outputs, geometry=index_geometry)
    for material, construct in finals:
        graph.sequence(material, SequencingPlan(
            library=construct, primers=primers,
            acquired_roles=tuple(p.role for p in primers)))

    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "tso", "template_switch", *primer_reviewed, "preamp_polymerase", "preamp",
        "transposome", "tn5", "tagment", "gap_polymerase", "gap_fill",
        "index_forward", "index_reverse", "index_polymerase", "index_pcr",
        "library",
    ]
    if protocol not in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        reviewed.append("internal_library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def render_page(protocol: str) -> str:
    graph = protocol_graph(protocol)
    graph_workflow = ("<h2>Step-by-step library generation</h2>\n"
                      + protocol_graph_panels(graph, cls="long"))
    if protocol in ("SMART-seq", "SMART-seq2", "FLASH-seq"):
        body = [ss2_oligos(protocol), graph_workflow, ss2_final()]
    else:
        body = [ss3_oligos(protocol), graph_workflow, ss3_final(protocol)]
    return "\n".join([head(f"{protocol} library chemistry"), preamble(protocol, complete=True),
                      *body, sequencing(protocol)])
