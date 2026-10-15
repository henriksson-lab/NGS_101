#!/usr/bin/env python3
"""
Self-test for the Markdown renderer and the docs build.

lib/mdrender.py is a hand-written renderer rather than a library dependency, so the
subset it supports is pinned here by example. The last section renders every real note
in the repo and asserts structural invariants against the Markdown source, which is
what catches a regression on content none of the unit cases cover.

Run:  python3 tools/selftest.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT))

from checks import Check, run_common  # noqa: E402
from adapter_ligation import cohesive_duplex  # noqa: E402
from chemdraw import Row, Segment, oligo, panel, protocol_graph_panels  # noqa: E402
from mdrender import render, rewrite_link, slug  # noqa: E402
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            MolecularRelation,
                            ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem,
                            Step, StepKind, Terminal)  # noqa: E402
from chemistry_validation import (Status, check_graph as check_chemistry,
                                  warnings as chemistry_warnings)  # noqa: E402

def _raise_fact():
    """Expanding a fact that names a missing attribute must raise (regression guard)."""
    from mdfacts import expand as _x
    try:
        _x("{{= chemdraw.NOT_A_REAL_NAME }}", "probe.md")
    except ValueError as e:
        return e
    return None


check = Check()
run_common(check)

check.section("restriction adapter compatibility")
from restriction import ECORI, MSPI  # noqa: E402
check("cohesive ends are derived from cut coordinates",
      (ECORI.overhang, MSPI.overhang), ("AATT", "CG"))
check.raises("a mismatched adapter overhang is rejected at construction",
             lambda: ECORI.adapter_end("CG"))
sticky = cohesive_duplex("round", "TTCGGCTTGGTTGCGCTTGGTTGC",
                         "GTCAGCAACCAAGCGCAACCAAGC")
check("source-written sticky duplex derives both 5-prime ends",
      (sticky.top_overhang, sticky.paired_core, sticky.bottom_overhang),
      ("TTCG", "GCTTGGTTGCGCTTGGTTGC", "GTCA"))

check.section("typed protocol graph")


def _complete_graph():
    graph = ProtocolGraph("example", authority="example authority")
    graph.add_input(Material("sample", "sample DNA", MaterialKind.BIOLOGICAL,
                             primary=True, source="example authority"))
    graph.add_input(Material("adapter", "loaded adapter", MaterialKind.ADAPTER,
                             source="example authority"))
    graph.add_input(Material("tn5", "Tn5", MaterialKind.CATALYST,
                             source="example authority"))
    library = Material("library", "library", MaterialKind.LIBRARY,
                       source="example authority")
    waste = Material("waste", "non-library products", MaterialKind.WASTE,
                     source="example authority")
    graph.add_step(Step(
        "tagment", "Tagment", StepKind.TAGMENTATION,
        (Input("sample", InputRole.TEMPLATE), Input("adapter", InputRole.STRUCTURAL),
         Input("tn5", InputRole.CATALYST)), ("library", "waste"),
        source="example authority"), (library, waste))
    graph.mark_terminal(Terminal("library", Disposition.SEQUENCING),
                        SequencingPlan(limitation="example source withholds primer sites"))
    graph.mark_terminal(Terminal("waste", Disposition.DISCARD, "not sequenceable"))
    for item, target in (("sample_source", "sample"), ("adapter_source", "adapter"),
                         ("tn5_source", "tn5"), ("tagmentation", "tagment"),
                         ("sequencing_output", "library")):
        graph.add_source_item(SourceItem(item, item.replace("_", " "),
                                         "example authority"))
        graph.cover_source_item(SourceCoverage(item, "modelled", target))
    return graph


check("a connected, source-covered graph with a sequencing terminal is publishable",
      _complete_graph().publication_problems(), [])
check.raises("tagmentation cannot hide its non-library output", lambda: ProtocolGraph(
    "bad", authority="source")._validate_step_shape(Step(
        "tagment", "Tagment", StepKind.TAGMENTATION,
        (Input("sample", InputRole.TEMPLATE), Input("adapter", InputRole.STRUCTURAL),
         Input("tn5", InputRole.CATALYST)), ("library",), source="source")))


def _uncovered_source():
    graph = _complete_graph()
    graph.add_source_item(SourceItem("cleanup", "cleanup", "example authority"))
    return graph.publication_problems()


check("unaccounted authoritative operations prevent publication",
      _uncovered_source(), ["authoritative source item 'cleanup' is not accounted for"])

check.raises("magic steps require an explanation", lambda: Step(
    "unknown", "Unknown chemistry", StepKind.MAGIC,
    (Input("sample", InputRole.TEMPLATE),), ("product",)))
check.raises("magic steps cannot assert a checked molecular relation", lambda: Step(
    "unknown", "Unknown chemistry", StepKind.MAGIC,
    (Input("sample", InputRole.TEMPLATE),), ("product",), reason="Gap in methods",
    relations=(MolecularRelation.STRAND_ORDER,)))

magic = ProtocolGraph("explained gaps", authority="example source")
raw = magic.input("raw", "raw material", MaterialKind.BIOLOGICAL, primary=True)
intermediate = magic.product("intermediate", "intermediate", MaterialKind.MOLECULAR)
library = magic.product("library", "library", MaterialKind.LIBRARY)
magic.react("gap_one", "First gap", StepKind.MAGIC,
            ((raw, InputRole.TEMPLATE),), (intermediate,),
            reason="Paper omits the conversion <mechanism>.")
magic.react("gap_two", "Second gap", StepKind.MAGIC,
            ((intermediate, InputRole.TEMPLATE),), (library,),
            reason="Product structure is not specified.")
magic.sequence(library, SequencingPlan(limitation="Primer sites are unavailable"))
magic.declare_reviewed_scope(("raw", "gap_one", "gap_two", "library"),
                             "example source")
check("explained gaps remain publishable with reviewed provenance",
      magic.publication_problems(), [])
check("validation reports every magic warning",
      [(w.target, w.message) for w in chemistry_warnings(magic)],
      [("gap_one", "Paper omits the conversion <mechanism>."),
       ("gap_two", "Product structure is not specified.")])
check("magic transitions remain unresolved",
      [c.status for c in check_chemistry(magic) if c.target.startswith("gap_")],
      [Status.UNRESOLVED, Status.UNRESOLVED])
magic_html = protocol_graph_panels(magic)
check("page shows and escapes every warning",
      (magic_html.count("<li>") == 2
       and "Paper omits the conversion &lt;mechanism&gt;." in magic_html
       and "Product structure is not specified." in magic_html), True)

check.section("independent molecular claims")
from chemdraw import Construct, Scene  # noqa: E402
from circular import CircularizationGeometry, circularize_ssdna, rolling_circle  # noqa: E402
from dataclasses import replace  # noqa: E402

linear = Construct([Segment("body", "ACGTTGCA")], "test strand")
circle = circularize_ssdna(linear, five_prime_phosphate=True)
chem = ProtocolGraph("synthetic chemistry")
substrate = chem.input("substrate", "linear strand", MaterialKind.MOLECULAR,
                       state=linear, primary=True)
ligase = chem.input("ligase", "ligase", MaterialKind.CATALYST)
closed = chem.product("closed", "circle", MaterialKind.MOLECULAR,
                      chemical_state=circle)
chem.react("close", "close", StepKind.CIRCULARIZATION,
           ((substrate, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
           (closed,), relations=(MolecularRelation.CIRCLE_CLOSURE,),
           geometry=CircularizationGeometry.known(
               "single-strand intramolecular ligation", statement="synthetic test"))
check("circle order is checked against the input strand",
      check_chemistry(chem)[0].status, Status.ESTABLISHED)

wrong_circle = circularize_ssdna(Construct([Segment("body", "ACGTTGCC")],
                                           "wrong strand"), five_prime_phosphate=True)
chem.materials["closed"] = replace(closed, chemical_state=wrong_circle)
check("a circular product with a changed base is contradictory",
      check_chemistry(chem)[0].status, Status.CONTRADICTORY)
chem.materials["closed"] = closed

polymerase = chem.input("polymerase", "polymerase", MaterialKind.CATALYST)
rca = rolling_circle(circle, "TGC", copies=2)
product = chem.product("product", "concatemer", MaterialKind.MOLECULAR,
                       chemical_state=rca)
chem.react("copy", "rolling copy", StepKind.EXTENSION,
           ((closed, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (product,),
           relations=(MolecularRelation.ROLLING_CIRCLE_COPY,))
check("RCA repeat derives from input circle",
      [c.status for c in check_chemistry(chem) if c.target == "copy"],
      [Status.ESTABLISHED])
chem.materials["product"] = replace(product, chemical_state=replace(rca, unit="A" * len(rca.unit)))
check("wrong RCA repeat is contradictory",
      [c.status for c in check_chemistry(chem) if c.target == "copy"],
      [Status.CONTRADICTORY])

scene = Scene(); scene.strand("a", [Segment("a", "AAAA")])
scene.anneal("b", [Segment("b", "TTTT")], to="a", pair=("b", "a"))
pair_graph = ProtocolGraph("pairing")
pair_graph.input("pair", "duplex", MaterialKind.MOLECULAR, state=scene)
check("exact scene pairing is established",
      check_chemistry(pair_graph)[0].status, Status.ESTABLISHED)
scene.strands["b"].segs[0].top = "CCCC"
check("checker detects a changed strand after scene construction",
      check_chemistry(pair_graph)[0].status, Status.CONTRADICTORY)

from restriction import Digest, MBOI  # noqa: E402
from topology import HairpinTopology  # noqa: E402

digest = Digest(MBOI, "AAAA", "CCCC")
fill = digest.fill_in(biotin_base="A")
cut_graph = ProtocolGraph("restriction fill")
cut = cut_graph.input("cut", "MboI ends", MaterialKind.MOLECULAR,
                      chemical_state=digest)
enzyme = cut_graph.input("enzyme", "polymerase", MaterialKind.CATALYST)
filled = cut_graph.product("filled", "blunt ends", MaterialKind.MOLECULAR,
                           chemical_state=fill)
cut_graph.react("fill", "fill", StepKind.EXTENSION,
                ((cut, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (filled,),
                relations=(MolecularRelation.RESTRICTION_FILL,))
check("restriction cut conserves both strand paths",
      [c.status for c in check_chemistry(cut_graph) if c.target == "cut"],
      [Status.ESTABLISHED])
class BrokenDigest(Digest):
    @property
    def left_top(self):
        return "T" + super().left_top

cut_graph.materials["cut"] = replace(cut, chemical_state=BrokenDigest(MBOI, "AAAA", "CCCC"))
check("an inconsistent restriction cut is contradictory",
      [c.status for c in check_chemistry(cut_graph) if c.target == "cut"],
      [Status.CONTRADICTORY])
cut_graph.materials["cut"] = cut
check("restriction fill derives both recessed tracts",
      [c.status for c in check_chemistry(cut_graph) if c.target == "fill"],
      [Status.ESTABLISHED])
cut_graph.materials["filled"] = replace(
    filled, chemical_state=replace(fill, left_new_top="AAAA"))
check("wrong fill tract is contradictory",
      [c.status for c in check_chemistry(cut_graph) if c.target == "fill"],
      [Status.CONTRADICTORY])

fold_graph = ProtocolGraph("symbolic fold")
unfolded = fold_graph.input("unfolded", "unknown strand", MaterialKind.MOLECULAR,
                            state=Construct([Segment("whole", "N" * 10,
                                                     placeholder=True)]))
temperature = fold_graph.input("temperature", "folding temperature",
                               MaterialKind.CONDITION)
topology = HairpinTopology("hairpin", "N" * 10, 0, 7, 3)
folded = fold_graph.product("folded", "hairpin", MaterialKind.MOLECULAR,
                            chemical_state=topology)
fold_graph.react("fold", "fold", StepKind.FOLDING,
                 ((unfolded, InputRole.TEMPLATE),
                  (temperature, InputRole.CONDITION)), (folded,),
                 relations=(MolecularRelation.HAIRPIN_STEM,))
check("unknown stem bases make folding conditional",
      check_chemistry(fold_graph)[0].status, Status.CONDITIONAL)

order_graph = ProtocolGraph("relation independent of step label")
order_input = order_graph.input("input", "DNA", MaterialKind.MOLECULAR, state=linear)
heat = order_graph.input("heat", "heat", MaterialKind.CONDITION)
order_output = order_graph.product("output", "heated DNA", MaterialKind.MOLECULAR,
                                  state=linear)
order_graph.react("heat_dna", "Heat DNA", StepKind.THERMAL_TREATMENT,
                  ((order_input, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                  (order_output,), relations=(MolecularRelation.STRAND_ORDER,))
check("same relation checks a different displayed step kind",
      [c.status for c in check_chemistry(order_graph) if c.target == "heat_dna"],
      [Status.ESTABLISHED])
order_graph.materials["output"] = replace(order_output, state=None)
check("declared relation without sufficient states stays unresolved and warns",
      ([c.status for c in check_chemistry(order_graph) if c.target == "heat_dna"],
       [(w.target, w.message.startswith("Unresolved"))
        for w in chemistry_warnings(order_graph)]),
      ([Status.UNRESOLVED], [("heat_dna", True)]))

check.section("molecular diagram renderers")
svg_panel = panel([Row(chunks=[("ACGT", "r1", False)])])
check("SVG is the default diagram renderer", '<svg class="chem-svg' in svg_panel)
check("SVG keeps sequence as selectable text", "<tspan class=\"chem-r1\">ACGT</tspan>" in svg_panel)
check("the legacy renderer remains explicitly available",
      panel([Row(chunks=[("ACGT", None, False)])], renderer="legacy").startswith("<pre>"))
text_oligo = oligo("PCR primer", [Segment("binding", "ACGTACGT")])
check("ordering oligos stay raw selectable text", "<svg" not in text_oligo)
check("ordering oligos retain explicit strand ends", "5'- ACGTACGT -3'" in text_oligo)


def r(md: str) -> str:
    return render(md)


check.section("headings and anchors")
check("an h2 gets an id", r("## Hello there"), '<h2 id="hello-there">Hello there</h2>')
check("levels 1-6 are all recognised",
      [r(f"{'#' * i} x").split(">")[0][1:3] for i in range(1, 7)],
      ["h1", "h2", "h3", "h4", "h5", "h6"])
check("a trailing closing hash run is not part of the text", r("## x ##"),
      '<h2 id="x">x</h2>')
check("seven hashes is not a heading", r("####### x").startswith("<p>"))
check("emoji and punctuation are dropped from the slug",
      slug("## 9. ✅ MIP worked on gDNA (why)"), "9-mip-worked-on-gdna-why")
check("a slug keeps digits and hyphens", slug("lentiCRISPR v1 - #49535"),
      "lenticrispr-v1---49535")
check("a slug never comes back empty", slug("✅ ⚠ 🔴"), "section")
check("inline markup in a heading still renders",
      r("## a `b` **c**"),
      '<h2 id="a-b-c">a <code>b</code> <strong>c</strong></h2>')

check.section("inline markup")
check("bold", r("**x**"), "<p><strong>x</strong></p>")
check("italic with asterisks", r("*x*"), "<p><em>x</em></p>")
check("italic with underscores", r("_x_"), "<p><em>x</em></p>")
check("bold wins over italic", r("**x**"), "<p><strong>x</strong></p>")
check("bold inside a word boundary is left alone", r("a_b_c"), "<p>a_b_c</p>")
check("an unmatched ** stays literal -- the notes use it as a footnote marker",
      r("Wash Buffer (WB)**"), "<p>Wash Buffer (WB)**</p>")
check("inline code is escaped", r("`a<b>c`"), "<p><code>a&lt;b&gt;c</code></p>")
check("markup inside code is NOT interpreted", r("`**x**`"),
      "<p><code>**x**</code></p>")
check("a code span may contain backticks via a longer fence", r("``a ` b``"),
      "<p><code>a ` b</code></p>")
check("asterisks inside code never leak out", "<em>" not in r("`a * b * c`"))

check.section("links")
check("a plain link", r("[t](http://x/)"), '<p><a href="http://x/">t</a></p>')
check("a .md link is rewritten to .html", r("[t](a/b.md)"),
      '<p><a href="a/b.html">t</a></p>')
check("an anchor on a .md link survives", rewrite_link("a.md#frag"), "a.html#frag")
check("a .md inside a longer name is not rewritten", rewrite_link("a.mdx"), "a.mdx")
check("http links are never rewritten", rewrite_link("http://x/a.md"), "http://x/a.md")
check("a bare anchor is left alone", rewrite_link("#x"), "#x")
check("a non-md repo path is left alone", rewrite_link("lib/umimodel.py"),
      "lib/umimodel.py")
check("bold inside a link label renders", r("[**t**](u.md)"),
      '<p><a href="u.html"><strong>t</strong></a></p>')
check("an autolink", r("<https://x/>"), '<p><a href="https://x/">https://x/</a></p>')
check("a bare URL is linkified", r("see https://x/y"),
      '<p>see <a href="https://x/y">https://x/y</a></p>')
check("trailing sentence punctuation stays outside the href",
      r("see https://x/y."), '<p>see <a href="https://x/y">https://x/y</a>.</p>')
check("a URL already inside a link is not linkified twice",
      r("[a](https://x/)").count("<a "), 1)

check.section("raw HTML: allowlist in, metavariables escaped")
for tag in ("br", "sub", "sup", "i", "b", "em", "strong"):
    check(f"<{tag}> passes through", f"<{tag}>" in r(f"a<{tag}>b"))
check("a closing allowlisted tag passes through", "</sub>" in r("T<sub>m</sub>"))
check("<cbc> is escaped, not treated as markup", r("a<cbc>b"), "<p>a&lt;cbc&gt;b</p>")
for mv in ("seqid", "align", "inf", "name", "id", "p5", "t7", "pre", "script"):
    check(f"<{mv}> is escaped", f"&lt;{mv}&gt;" in r(f"x<{mv}>y"))
check("an ampersand is escaped once, not twice", r("a & b"), "<p>a &amp; b</p>")
check("an already-escaped entity is not double-escaped into gibberish",
      "&amp;amp;" in r("a &amp; b"))

check.section("tables")
t = r("| a | b |\n|---|---|\n| 1 | 2 |")
check("a table renders as a table", t.count("<table>"), 1)
check("it is wrapped so it can scroll instead of widening the page",
      '<div class="tw">' in t)
check("the header row becomes th", t.count("<th>"), 2)
check("the body row becomes td", t.count("<td>"), 2)
check("left alignment adds no style", 'style=' not in r("|a|\n|:--|\n|1|"))
check("right alignment is honoured",
      'style="text-align:right"' in r("|a|\n|--:|\n|1|"))
check("centre alignment is honoured",
      'style="text-align:center"' in r("|a|\n|:-:|\n|1|"))
check("alignment applies to body cells too",
      r("|a|\n|--:|\n|1|").count('text-align:right'), 2)
check("inline markup works inside a cell",
      "<code>x</code>" in r("| a |\n|---|\n| `x` |"))
check("a <br> works inside a cell -- the notes use it for stacked values",
      "<br>" in r("| a |\n|---|\n| x<br>y |"))
# The notes use an empty header row for key/value tables with no column titles.
hl = r("| | |\n|---|---|\n| k | v |")
check("a headerless table (empty header row) still renders as one table",
      hl.count("<table>"), 1)
check("...with the empty cells as th, not as an alignment row",
      hl.count("<th>"), 2)
check("...and the data row intact", "<td>k</td><td>v</td>" in hl)
check("a pipe row with no alignment row is just a paragraph",
      r("| a | b |").startswith("<p>"))
check("a --- line is a rule, not a one-column alignment row",
      r("text\n\n---\n\nmore").count("<hr>"), 1)
check("a table ends at the first blank line",
      r("|a|\n|---|\n|1|\n\nafter").count("<table>"), 1)

check.section("code blocks")
check("a fence becomes pre/code", r("```\nx\n```"), "<pre><code>x</code></pre>")
check("a language is recorded as a class",
      r("```python\nx\n```"), '<pre><code class="language-python">x</code></pre>')
check("content inside a fence is escaped",
      r("```\n<cbc>&\n```"), "<pre><code>&lt;cbc&gt;&amp;</code></pre>")
check("blank lines and indentation inside a fence are preserved exactly",
      r("```\na\n\n  b\n```"), "<pre><code>a\n\n  b</code></pre>")
check("markdown inside a fence is not interpreted",
      r("```\n## x\n| a |\n|---|\n```"), "<pre><code>## x\n| a |\n|---|</code></pre>")
check("an indented fence inside a list item is still a fence -- "
      "scg_lib_structs_style.md:110 relies on this",
      r("- item\n\n   ```html\n   <b>\n   ```").count("<pre><code"), 1)
check("a tilde fence works too", r("~~~\nx\n~~~"), "<pre><code>x</code></pre>")

check.section("lists")
check("a tight unordered list", r("- a\n- b"), "<ul>\n<li>a</li>\n<li>b</li>\n</ul>")
check("an ordered list", r("1. a\n2. b"), "<ol>\n<li>a</li>\n<li>b</li>\n</ol>")
check("an ordered list starting at 3 keeps its numbering",
      'start="3"' in r("3. a\n4. b"))
check("a nested list nests", r("- a\n  - b").count("<ul>"), 2)
check("a list ends at a following paragraph",
      r("- a\n\nafter").count("<li>"), 1)
check("a lazy continuation line joins its item, not a new one",
      r("1. a\n   still a").count("<li>"), 1)
check("...and its text is joined into one paragraph",
      r("1. a\n   still a"), "<ol>\n<li>a still a</li>\n</ol>")
check("a 5-space continuation joins too -- ira1.md:302 relies on this",
      r("1. a\n     still a").count("<li>"), 1)
check("switching marker type starts a separate list",
      r("- a\n\n1. b").count("<ol>"), 1)
check("inline markup works in a list item", "<strong>x</strong>" in r("- **x**"))
check("a table inside a list item renders",
      r("- item\n\n  | a |\n  |---|\n  | 1 |").count("<table>"), 1)

check.section("blockquotes and rules")
check("a blockquote", r("> x"), "<blockquote>\n<p>x</p>\n</blockquote>")
check("a multi-line blockquote is one quote",
      r("> a\n> b").count("<blockquote>"), 1)
check("markup inside a blockquote renders", "<strong>x</strong>" in r("> **x**"))
check("a table inside a blockquote renders",
      r("> | a |\n> |---|\n> | 1 |").count("<table>"), 1)
check("--- is a horizontal rule", r("---"), "<hr>")
check("*** is a horizontal rule too", r("***"), "<hr>")
check("a rule does not swallow the text after it",
      r("---\nafter"), "<hr>\n<p>after</p>")

check.section("paragraphs")
check("wrapped source lines join into one paragraph",
      r("a\nb\nc"), "<p>a b c</p>")
check("a blank line separates paragraphs", r("a\n\nb").count("<p>"), 2)
check("a heading interrupts a paragraph without needing a blank line",
      r("a\n## b").count("<h2"), 1)
check("empty input gives empty output", r(""), "")
check("whitespace-only input gives empty output", r("\n\n  \n"), "")

check.section("no placeholder ever escapes")
for sample in ("`a` **b** [c](d.md) <https://x/> T<sub>m</sub> `**x**`",
               "| `a` | **b** |\n|---|---|\n| [c](d.md) | <br> |",
               "> `a` and [b](c.md)",
               "- `a`\n- [b](c.md)"):
    out = r(sample)
    check(f"no \\x00/\\x01 left in {sample[:28]!r}...",
          "\x00" not in out and "\x01" not in out)

check.section("every real note in the repo renders consistently")
SKIP = {".git", "_data", "pdf", "__pycache__"}
sys.path.insert(0, str(ROOT))
import build_docs  # noqa: E402

notes = sorted(p for p in ROOT.rglob("*.md") if not SKIP & set(p.relative_to(ROOT).parts)
               and str(p.relative_to(ROOT)) not in build_docs.NOT_NOTES)   # agent instructions
check("the repo's notes are found", len(notes) >= 30)

PAIRED = ("table", "tr", "td", "th", "ul", "ol", "li", "blockquote", "p", "pre",
          "code", "strong", "em", "a", "h1", "h2", "h3", "h4")
bad_balance, bad_tables, bad_fences, leaked, unconsumed = [], [], [], [], []

for p in notes:
    src = p.read_text(encoding="utf-8")
    out = render(src)
    rel = p.relative_to(ROOT)
    for tag in PAIRED:
        if len(re.findall(rf"<{tag}[ >]", out)) != out.count(f"</{tag}>"):
            bad_balance.append(f"{rel} <{tag}>")
    # count source tables and fences, skipping anything inside a fence
    n_tbl = n_fence = 0
    infence = False
    for line in src.split("\n"):
        if re.match(r"^\s*(```|~~~)", line):
            infence = not infence
            n_fence += 1
            continue
        if infence:
            continue
        # A real alignment row has a dash in EVERY cell. `| | |` is the notes'
        # headerless-table idiom -- an empty header row, not an alignment row.
        if re.match(r"^\s*\|?(\s*:?-+:?\s*\|)+\s*:?-*:?\s*$", line):
            n_tbl += 1
    if n_tbl != out.count("<table>"):
        bad_tables.append(f"{rel}: {n_tbl} vs {out.count('<table>')}")
    if n_fence // 2 != out.count("<pre><code"):
        bad_fences.append(f"{rel}: {n_fence // 2} vs {out.count('<pre><code')}")
    if "\x00" in out or "\x01" in out:
        leaked.append(str(rel))
    # strip code, then nothing markdown-ish may remain
    body = re.sub(r"<pre><code[^>]*>.*?</code></pre>|<code>.*?</code>", "", out, flags=re.S)
    if re.search(r"^\s*#{1,6} ", body, re.M) or re.search(r"\]\(\S+\)", body) \
            or re.search(r"\*\*[^*\n]+\*\*", body):
        unconsumed.append(str(rel))

check("every note's tags balance", bad_balance, [])
check("every source table becomes exactly one HTML table", bad_tables, [])
check("every source fence becomes exactly one code block", bad_fences, [])
check("no internal placeholder leaks into any note", leaked, [])
check("no unconsumed markdown is left in any note", unconsumed, [])

check.section("the build script and its output")
import build_docs  # noqa: E402

check("build_docs finds the same notes", len(build_docs.notes()), len(notes))
check("it excludes _data and pdf",
      not any({"_data", "pdf"} & set(p.parts) for p in build_docs.notes()))
check("plain() strips a link to its label", build_docs.plain("a [b](c) d"), "a b d")
check("plain() strips code and emphasis ticks",
      build_docs.plain("`a` **b** *c*"), "a b c")
check("a TOC appears only once there are enough headings",
      build_docs.toc_html([(2, "a", "a")]), "")
check("...and lists h2 and h3 but not h4",
      build_docs.toc_html([(2, "a", "a"), (3, "b", "b"), (2, "c", "c"),
                           (2, "d", "d"), (4, "e", "e")]).count("<li"), 4)
check("a TOC entry is plain text, never a nested <a>",
      "](" not in build_docs.toc_html([(2, f"h{i} [x](y)", f"h{i}") for i in range(5)]))

one = build_docs.build_page(ROOT / "README.md")
check("a built page declares a title", one.startswith("<title>"))
check("it sets an explicit body background, so the host ground cannot show through",
      "background:var(--bg)" in one)
check("it defines the light palette on bare :root", "\n:root {" in one)
check("it redefines the palette for system dark",
      '@media (prefers-color-scheme: dark)' in one and ':root:not([data-theme="light"])' in one)
check("it redefines the palette for an explicit dark choice",
      ':root[data-theme="dark"]' in one)
check("tables are wrapped so the page never scrolls sideways", 'class="tw"' in one)
check("it links back to the front page", 'href="index.html"' in one)
check("a nested note's breadcrumb walks up to the root",
      'href="../index.html"' in build_docs.build_page(
          ROOT / "crispr-mip__10.1101+2024.03.28.587082" / "01_crispr-mip.md"))
check("a twice-nested note walks up twice",
      'href="../../index.html"' in build_docs.build_page(
          ROOT / "crispr-mip__10.1101+2024.03.28.587082" / "ref" / "MANIFEST.md"))
check("the rendered HTML is on disk for every note",
      [str(p.relative_to(ROOT)) for p in notes if not p.with_suffix(".html").exists()], [])
check("rendering is deterministic -- building twice gives the same bytes",
      build_docs.build_page(ROOT / "README.md"), one)

import build_index  # noqa: E402

_front, _ = build_index.build(run_checks=False)
_public = [p for p in build_index.collect(run_checks=False) if p["section"] == "published"]
check("the public index retains the complete published catalogue",
      _front.count('<tr class="pr"'), len(_public))
check("a protocol title opens its schematic, not its first reference note",
      ('href="smart-seq2__10.1038+nmeth.2639/smart-seq2.html">SMART-seq2</a>' in _front,
       'href="smart-seq__10.1038+nbt.2282/01_smart-seq.html"' in _front),
      (True, False))
check("work in progress and the all-notes index are absent from the public front page",
      ("Work in progress" in _front, 'href="notes.html"' in _front), (False, False))
check("reference-note bodies are not embedded in public search data",
      "Evidence marking" in _front, False)

# ------------------------------------------------- computed facts inside the notes
check.section("computed facts in notes expand (lib/mdfacts.py)")
from mdfacts import FACT, expand  # noqa: E402

_notes = build_docs.notes()
_withfacts = [p for p in _notes if FACT.search(p.read_text())]
check("at least one note carries computed facts", len(_withfacts) > 0)
for _p in _withfacts:
    _rel = _p.relative_to(ROOT)
    try:
        expand(_p.read_text(), str(_rel))
        _err = None
    except ValueError as e:
        _err = str(e)
    check(f"{_rel}: every {{{{fact}}}} resolves", _err, None)
check("a fact naming something that no longer exists is an error, not silent",
      any("failed" in str(_e) for _e in [_raise_fact()]))

# ---------------------------------------------------- drawn duplexes must base-pair
check.section("every generated page: drawn duplex columns base-pair (tools/lint_pairing.py)")
sys.path.insert(0, str(HERE))
import lint_pairing  # noqa: E402

for page in sorted(p for p in ROOT.rglob("*.html")
                   if not (set(p.relative_to(ROOT).parts) & lint_pairing.SKIP)):
    bad = lint_pairing.lint_text(page.read_text())
    check(f"{page.relative_to(ROOT)}: no unpaired duplex columns",
          [f"{c} | {a.strip()} / {b.strip()}" for c, a, b, _ in bad], [])

# --------------------------------------- third-party material has provenance recorded
check.section("every ref/ directory documents where its files came from")
for _ref in sorted(ROOT.glob("*/ref")) + sorted(ROOT.glob("*/ref/*")):
    if not _ref.is_dir():
        continue
    _files = [f for f in _ref.iterdir() if f.is_file() and f.suffix not in (".html",)]
    if not _files or all(f.name == "MANIFEST.md" for f in _files):
        continue
    check(f"{_ref.relative_to(ROOT)}: has a MANIFEST.md",
          (_ref / "MANIFEST.md").exists()
          or (_ref.parent / "MANIFEST.md").exists())
check("no MANIFEST hard-codes an absolute path from this machine",
      [str(m.relative_to(ROOT)) for m in ROOT.rglob("ref/**/MANIFEST.md")
       if "/Users/" in m.read_text()], [])

# ------------------------------------------------- sequence triage (tools/scrape_primers.py)
check.section("oligo scraping from source documents (tools/scrape_primers.py)")
import scrape_primers as spx  # noqa: E402

_doc = ("Adapters. RA3 (rAppTGGAATTCTCGGGTGCCAAGG-ddC) is pre-adenylated.\n"
        "RA5 (NH2-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrA) is RNA.\n"
        "RP1 (AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTAC\nAGTCCGA) is the forward primer.\n"
        "TSO: 5'-/5Biosg/AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3'\n"
        "Vector ctctagaGATCGGAAGAGCACACGT and A*C*G*TACGTACGT/iSp18/ACGTACGTAC/3Phos/.\n"
        "TAGGED CATS GATHER DATA, and acgtacgtacgtacgt in lower case.\n"
        "Again RP1 AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA here.\n")
_hits = {h.label or h.seq: h for h in spx.scan(Path("doc"), _doc, spx.known_sequences())}
check("RA3: bases without the modifications", _hits["RA3"].seq, "TGGAATTCTCGGGTGCCAAGG")
check("RA3: written with them", _hits["RA3"].written, "rAppTGGAATTCTCGGGTGCCAAGG-ddC")
check("RA3: 5' adenylation and 3' dideoxy C are named", _hits["RA3"].mods,
      ["5' adenylated (rApp)", "3' dideoxy (ddC)"])
check("RA5: per-base RNA marks rejoin into one sequence", _hits["RA5"].seq,
      "GUUCAGAGUUCUACAGUCCGA")
check("RA5: amino linker and RNA bases are named", _hits["RA5"].mods,
      ["5' amine (NH2)", "RNA bases x21"])
check("RP1: a sequence broken by PDF extraction is rejoined", len(_hits["RP1"].seq), 50)
check("RP1: and recognised as containing P5", "illumina.P5" in _hits["RP1"].known)
check("RP1: a repeat is listed as a further line, not a second hit", _hits["RP1"].also, [8])
check("TSO: IDT 5' biotin and the rGrG+G 3' end are kept", _hits["TSO"].mods,
      ["5' biotin (/5Biosg/)", "3' rGrG+G"])
check("mixed case: vector context in lower case stays attached",
      any(h.seq == "ctctagaGATCGGAAGAGCACACGT" for h in _hits.values()))
check("phosphorothioates and internal spacers are named",
      [h.mods for h in _hits.values() if h.seq.startswith("ACGTACGT")],
      [["3' phosphate (/3Phos/)", "phosphorothioate x3", "internal /iSp18/"]])
check("upper-case English words and all-lower-case runs are not sequences",
      [h.seq for h in _hits.values() if "CATS" in h.seq or h.seq.islower()], [])
check("--find matches the reverse strand across a line break",
      spx.find(_doc, "TCGGACTGTAGAACTCTGAACG"), [3, 8])

# Regression corpus for the source forms described in tools/scrape_primers.py.
_d2 = ("STRT-V3-T30 (5'-biotin-AAGCAGTGGTATCAACGCAGAGTCGACT30VN-3')\n"
       "STRT-V2-n (5'-AAGCAGTGGTATCAACGCAGAGTGCAGTGCTXXXXXXrGrGrG-3')\n"
       "DI-P1A-idx (5'Bio-AATGATACGGCGACCACCGAGATCTACAC-XXXXX-CTACACGACGCTCTTCCGATC)\n"
       "RT (5'-/5Phos/CAGAGCNNNNNNNN[10bp barcode]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3')\n"
       "CEL GAGTTCTACAGTCCGACGATC[8 base\nbarcode]TTTTTTTTTTTTTTTTTTTTTTTTV\n"
       "BEAD CTACACGACGCTCTTCCGATCT-N16-N12-TTTCTTATATrGrGrG and (dT)30\n"
       "Idx CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG\n"
       "C1-P1-RNA-TSO\tBio-AAUGAUACGGCGACCACCGAUNNNNNGGG\n"
       "P7-1\tCAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG\n"
       "P5-gMac\taatgatacggcgaccaccgagatctacattgtatagaattcgcggccgctcgcgaT*A*c\n"
       "Randomer (5' TCA GAC GTG TGC TCT TCC GAT CTNNNNNNNNN 3')\n"
       "TSO-P\tAAGCAGTGGTATCAACGCAGAGT\n"
       "TSO AAGCAGTGGTATCAACGCAG\nAG TGAATrGrGrG\n"
       "S-P7\tCAAGCAGAAGACGGCATACGAGAT[NNNNNN]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T\n"
       "oligo1 ([Btn]CTACACGACGCTCTTCCGATCTNNNNNNNNN) and E5V6NEXT 5'-iCiGiCACACTCTTTCCCTACACGACGCrGrGrG-3'\n"
       "TTTCTTATATGGGcDNA (101 nt) and T30 oligo in prose\n")
_h2 = {h.label: h for h in spx.scan(Path("doc"), _d2, spx.known_sequences())}
check("T30VN shorthand expands inside the oligo; written stays as in the source",
      (_h2["STRT-V3-T30"].seq.count("T") >= 30, _h2["STRT-V3-T30"].seq.endswith("VN"),
       _h2["STRT-V3-T30"].written), (True, True, "5'-biotin-AAGCAGTGGTATCAACGCAGAGTCGACT30VN"))
check("XXXXXX placeholder is one N each, and the rGrGrG 3' end is kept",
      (_h2["STRT-V2-n"].seq[-6:], _h2["STRT-V2-n"].mods), ("NNNNNN", ["3' rGrGrG"]))
check("a hyphenated XXXXX index does not split the oligo",
      _h2["DI-P1A-idx"].length, 55)
check("[10bp barcode] expands to 10 N; the oligo stays whole with its 5' phosphate",
      (_h2["RT"].length, _h2["RT"].mods), (54, ["5' phosphate (/5Phos/)"]))
check("[8 base\\nbarcode] wrapped over a line is still one placeholder",
      "NNNNNNNNTTTT" in _h2["CEL"].seq)
check("N16-N12 shorthand keeps a bead oligo whole", _h2["BEAD"].length, 60)
check("an uncounted placeholder ([i7]) is a gap '…' and the length reads N+ nt",
      (spx.GAP in _h2["Idx"].seq, _h2["Idx"].length_text), (True, "39+ nt"))
check("U without r prefix is flagged as RNA; 'Bio-' is named biotin",
      _h2["C1-P1-RNA-TSO"].mods, ["5' biotin (Bio)", "RNA (written with U)"])
check("lower-case index inside a table row is part of the oligo",
      _h2["P7-1"].seq, "CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG")
check("an all-lower-case table cell with marks is an oligo", _h2["P5-gMac"].length, 58)
check("codon-spaced triplets are one sequence", _h2["Randomer"].length, 32)
check("'TSO-P' is a name, not a 5' phosphate", _h2["TSO-P"].mods, [])
check("a stray space and a 2-base piece after a wrap are bridged", _h2["TSO"].length, 27)
check("-s- phosphorothioate is looked through and [NNNNNN] read as bases",
      (_h2["S-P7"].length, _h2["S-P7"].mods), (64, ["phosphorothioate x1"]))
check("[Btn] is a 5' biotin, not a name or bases", _h2["oligo1"].mods, ["5' biotin ([Btn])"])
check("iCiGiC iso-bases are a 5' modification, not a C",
      _h2["E5V6NEXT"].mods, ["5' iso-dC/iso-dG (iCiGiC)", "3' rGrGrG"])
check("'cDNA' glued to a sequence is not bases; T30 in prose is not an oligo",
      sorted(h.seq for h in _h2.values() if h.seq.startswith(("TTTCTTATAT", "TTTTTTTTTT"))),
      ["TTTCTTATATGGG"])
_lst = "".join(f"AAAAAAAAAAAANNNNNN{b}CGATCGTGTCACCGA\n" for b in ("CCTAAA", "TGGAAT", "TTGAGC"))
check("one-sequence-per-line lists of equal length are not joined",
      [h.length for h in spx.scan(Path("l"), _lst, {})], [39, 39, 39])
check("a drawn duplex (second strand indented under the first) is not joined",
      [h.seq for h in spx.scan(Path("d"), "5'- TCGTCGGCAGCGTCAGATGTGTAT\n"
                                         "        AGCAGTCTACACATA -5'\n", {})],
      ["TCGTCGGCAGCGTCAGATGTGTAT", "AGCAGTCTACACATA"])
_fam = "".join(f"P7-{i}\tCAAGCAGAAGACGGCATACGAGAT{b}GTCTCGTGGGCTCGG\n"
               for i, b in enumerate(["ccgaatccga", "ataagccgga", "ccggcggcga", "ggcttgccaa"], 1))
_c = spx.collapse(spx.scan(Path("f"), _fam, {}))
check("a barcode table collapses into one family hit with its variable window",
      [(len(h.family), h.windows) for h in _c], [(4, [(25, 33)])])
check("--find: query N and text N, U = T, lower case, IUPAC S, '*' and line wraps",
      [spx.find("x\tAAUGAUACGGCGACCACCGAUNNNNNGGG\n", "AATGATACGGCGACCACCGATNNNNNGGG"),
       spx.find("P7-1\tCAAGCAGAAGACGGCATACGAGATccgaatccgaGTC\n", "caagcagaagacggcatacgagatCCGAATCCGAgtc"),
       spx.find("SB84 GCCAGASACGTTAGGCAGGACCTAACGT\n", "GCCAGACACGTTAGGCAGGACCTAACGT"),
       spx.find("AATGATACGGCGACC\nACCGAGATCTACACGCCTGTCCGCG-GAAGCAG TGGTATCAACGCAGAGT∗A∗C\n",
                "AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")],
      [[2 - 1], [1], [1], [1]])
check("--find: a pasted query with ends, mods, marks and a placeholder",
      spx.find(_d2, "5'-/5Phos/CAGAGCNNNNNNNN[10bp barcode]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3'"), [4])
check("--find: rGrGrG / +rG marks and -s- in the text are looked through",
      [spx.find("ACACTCTTTCCCTACACGACGC+rGrGrG\n", "ACACTCTTTCCCTACACGACGCGGG"),
       spx.find(_d2, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")], [[1], [16]])
check("--find: a following rG modification does not extend a completed primer hit",
      spx.find("ACACTCTTTCCCTACACGACGCrGrGrG\n", "ACACTCTTTCCCTACACGACGC"), [1])
check("--find: internal biotin-dT counts as its T base",
      [(strand, exact) for _, strand, exact, _ in spx.find(
          "A1\tGATCGAGCTCGAGAA/iBiodT/T\n", "GATCGAGCTCGAGAATT", detail=True)],
      [("+", "exact")])
check("--find: adjacent (N)10(T)30 shorthands stay one oligo",
      spx.find("E3V6NEXT\tACACTCTTTCCCTACACGACGCTCTTCCGATCT[BC6](N)10(T)30VN\n",
               "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" + "N" * 16 + "T" * 30 + "VN"),
      [1])
check("--find: tags exact vs only through IUPAC letters; mostly-N text is not a hit",
      [t[2] for t in spx.find("AAGCAGTGGTATCAACGCAGAGTGANNNGG\nNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n",
                              "AAGCAGTGGTATCAACGCAGAGTGAATGG", detail=True)], ["iupac"])
check("--find: a strand written 3'->5' is found as 'rev'",
      [t[1] for t in spx.find("3'-TCTCTTACTCCTTGGGCCCCGTC-5'\n", "CTGCCCCGGGTTCCTCATTCTCT",
                              detail=True)], ["rev"])
from source_oligo_check import _exact_hit  # noqa: E402
check("source guard accepts a printed 3'-to-5' strand as an ordered 5'-to-3' input",
      _exact_hit("P2 lower,3'-TCT CTT ACT CCT TGG GCC CCG TC-5'\n",
                 "CTGCCCCGGGTTCCTCATTCTCT"), True)

check.section("documents to text (tools/doctext.py)")
import tempfile  # noqa: E402
import zipfile  # noqa: E402

import doctext  # noqa: E402

with tempfile.TemporaryDirectory() as _tmp:
    _x = Path(_tmp) / "t.xlsx"
    with zipfile.ZipFile(_x, "w") as z:
        z.writestr("xl/workbook.xml", '<workbook><sheets><sheet name="Primers" sheetId="1"/>'
                                      '</sheets></workbook>')
        z.writestr("xl/sharedStrings.xml", "<sst><si><t>RP1</t></si>"
                   "<si><t>AATGATACGGCGACCACCGAGATCTACAC</t></si></sst>")
        z.writestr("xl/worksheets/sheet1.xml", '<worksheet><sheetData><row r="1">'
                   '<c r="A1" t="s"><v>0</v></c><c r="B1" t="s"><v>1</v></c>'
                   '<c r="C1"><v>50</v></c></row></sheetData></worksheet>')
    _d = Path(_tmp) / "t.docx"
    with zipfile.ZipFile(_d, "w") as z:
        z.writestr("word/document.xml", "<w:document><w:body><w:p><w:r><w:t>Oligos</w:t>"
                   "</w:r></w:p><w:tbl><w:tr><w:tc><w:p><w:r><w:t>RA3</w:t></w:r></w:p></w:tc>"
                   "<w:tc><w:p><w:r><w:t>rAppTGGAATTCTCGGGTGCCAAGG/3ddC/</w:t></w:r></w:p>"
                   "</w:tc></w:tr></w:tbl></w:body></w:document>")
    check("xlsx: each sheet under a heading, cells tab-separated",
          doctext.read(_x), "## sheet Primers\nRP1\tAATGATACGGCGACCACCGAGATCTACAC\t50\n")
    check("docx: paragraphs on their own lines, table cells tab-separated",
          [ln for ln in doctext.read(_d).splitlines() if ln.strip()],
          ["Oligos", "RA3 \trAppTGGAATTCTCGGGTGCCAAGG/3ddC/ \t"])
    _t = doctext.convert(_x)
    check("convert writes FILE.txt beside FILE", _t, Path(_tmp) / "t.xlsx.txt")
    check("...and the scraper reads the twin, not the original",
          [f.name for f in spx.files([_tmp])], ["t.docx", "t.xlsx.txt"])

check.report()
