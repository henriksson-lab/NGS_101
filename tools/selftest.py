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
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT))

from checks import Check, run_common  # noqa: E402
from adapter_ligation import cohesive_duplex  # noqa: E402
from chemdraw import Construct, Row, Segment, oligo, panel, protocol_graph_panels, revcomp  # noqa: E402
from mdrender import render, rewrite_link, slug  # noqa: E402
from protocol_graph import (AnnotationSlot, DetailView, Disposition, Highlight,
                            Input, InputRole, Material, MaterialKind, MolecularRelation,
                            ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem,
                            SourceClaim, SourceDiscrepancy, Step, StepKind,
                            Terminal)  # noqa: E402
from protocol_views import (OffsetPrimerSection, OligoEntry, OligoSection,
                            PrimerBindingSection, ReadLayoutSection, SourceRead,
                            WorkflowSection,
                            assemble_page, render_annotations)  # noqa: E402
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

from unittest.mock import patch
import build_index
with patch.object(build_index.subprocess, "run", return_value=type("FailedSuite", (), {
        "stdout": "all 3 checks passed\n", "returncode": 1})()):
    check("a successful check banner cannot hide a failing suite exit status",
          build_index.checks("char-seq__10.7554+eLife.27024"), (0, False))

check.section("IUPAC nucleotide-class complements")
from chemdraw import iupac_complement  # noqa: E402
check("all IUPAC classes complement in place and preserve case",
      iupac_complement("ACGTURYSWKMBDHVNacgturyswkmbdhvn"),
      "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
_iupac = "ACGTRYSWKMBDHVN"
check("DNA class complement is an involution",
      iupac_complement(iupac_complement(_iupac)), _iupac)
check.raises("a schematic placeholder is not a nucleotide class",
             lambda: iupac_complement("[insert]"))
check("an explicit IUPAC bottom is preserved by segment complement",
      Segment("degenerate bases", "HV", placeholder=True, length_bp=2,
              bottom=iupac_complement("HV")).bottom_text(),
      "DB")

check.section("sequencing preview across an unknown-length insert")
import seqprimers as sp  # noqa: E402
left_site = "ACGTTGCAACGTTGCAAGGA"
right_site = "GATCGTACCTGATCGTACCT"
preview_lib = Construct([
    Segment("left site", left_site), Segment("left flank", "AG"),
    Segment("variable insert", "[variable insert]", placeholder=True),
    Segment("right flank", "CT"), Segment("right site", right_site),
], "unknown insert preview")
left_read = sp.locate(preview_lib, sp.custom("Read 1", "left", left_site, "test"))
right_read = sp.locate(preview_lib, sp.custom("Read 2", "right", revcomp(right_site), "test"))
check("Read 1 preview stops before unknown insert", left_read.reads, "AG…")
check("Read 2 preview stops before unknown insert without reversing its label",
      right_read.reads, "AG…")
check("primer landing still reports the flanking segment",
      (left_read.reads_from, right_read.reads_from), ("left flank", "right flank"))
preview_html = sp.section(
    preview_lib,
    (sp.custom("Read 1", "left", left_site, "test"),
     sp.custom("Read 2", "right", revcomp(right_site), "test")),
    required_roles=("Read 1", "Read 2"), show_read_content=False)
check("first-bases table renders one unknown-region marker",
      (preview_html.count("<code>AG…</code> (from"),
       "<code>AG…</code>&hellip;" in preview_html), (2, False))

check.section("read contents across an index window of unstated length")
import illumina as il  # noqa: E402
from chemdraw import feature  # noqa: E402
unstated_i7_lib = Construct([
    Segment("P5", il.P5[:-4]), Segment("Read 1", il.TRUSEQ_READ1),
    Segment("insert", "ACGTACGTAC"), Segment("Read 2 arm", il.INDEX1_PRIMER),
    Segment("i7 index", "[i7 index; length not stated]", placeholder=True,
            feature=feature("i7", "sample_index", "unknown")),
    Segment("P7", il.P7_RC)], "unstated i7 window")
unstated_rows = {r.read: r for r in sp.read_contents(
    unstated_i7_lib, (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"]))}
check("a featured role-token index window stops Index 1 without a cycle count",
      (unstated_rows["Index 1 (i7)"].cycles, unstated_rows["Index 1 (i7)"].content),
      ("not fixed", "sample index (from cycle 1; length not stated)"))
unknown_insert_lib = Construct([
    Segment("Read 1", il.TRUSEQ_READ1),
    Segment("variable insert", "X" * 20, placeholder=True),
    Segment("downstream identifier", "NNNN", placeholder=True,
            feature=feature("downstream", "umi", "random"))], "unknown insert length")
check("identifier cycles cannot be placed beyond a variable-length insert",
      sp.feature_spans(unknown_insert_lib, (sp.TRUSEQ["R1"],), {"Read 1": 100}), [])
known_insert_lib = Construct([
    unknown_insert_lib.segments[0],
    Segment("fixed insert", "X" * 20, placeholder=True, length_bp=40),
    unknown_insert_lib.segments[2]], "known insert length")
check("a source-known insert length allows downstream identifier placement",
      [(s.cycle_start, s.cycle_end) for s in sp.feature_spans(
          known_insert_lib, (sp.TRUSEQ["R1"],), {"Read 1": 100})], [(41, 44)])
check("read content uses molecular length rather than compressed drawing width",
      "UMI (cycles 41–44)" in sp.read_contents(
          known_insert_lib, (sp.TRUSEQ["R1"],), {"Read 1": 100})[0].content, True)

check.section("restriction adapter compatibility")
branch_graph = ProtocolGraph("two contributing preparations", authority="test source")
branch_sample = branch_graph.input("sample", "input DNA", MaterialKind.BIOLOGICAL, primary=True)
branch_condition = branch_graph.input("condition", "physical fragmentation condition", MaterialKind.CONDITION)
branch_a = branch_graph.product("a", "first aliquot", MaterialKind.MOLECULAR,
                                state=Construct([Segment("first sequence", "ACGTACGTACGT")]))
branch_b = branch_graph.product("b", "second aliquot", MaterialKind.MOLECULAR,
                                state=Construct([Segment("second sequence", "TTGCTTGCTTGC")]))
branch_graph.react("prep_a", "Prepare first aliquot", StepKind.FRAGMENTATION,
                   ((branch_sample, InputRole.TEMPLATE), (branch_condition, InputRole.CONDITION)), (branch_a,))
branch_graph.react("prep_b", "Prepare second aliquot", StepKind.FRAGMENTATION,
                   ((branch_sample, InputRole.TEMPLATE), (branch_condition, InputRole.CONDITION)), (branch_b,))
branch_library = branch_graph.product("library", "pooled library", MaterialKind.LIBRARY,
                                      state=branch_a.state)
branch_graph.react("pool", "Pool preparations", StepKind.POOL,
                   ((branch_a, InputRole.TEMPLATE), (branch_b, InputRole.TEMPLATE)), (branch_library,))
branch_graph.sequence(branch_library, SequencingPlan(limitation="Test endpoint; primer sites unresolved"))
branch_graph.declare_reviewed_scope(("sample", "condition", "prep_a", "prep_b", "pool", "library"), "test source")
check("matching the displayed output to one input must not hide another input's preparation",
      "Prepare second aliquot" in protocol_graph_panels(branch_graph), True)
from types import SimpleNamespace
from batch_page import _render_terminal_sequencing
mixed_graph = ProtocolGraph("mixed read plans", authority="test source")
mixed_graph.input("known", "represented branch", MaterialKind.LIBRARY, state=known_insert_lib)
mixed_graph.input("unknown", "unresolved branch", MaterialKind.LIBRARY)
mixed_parts = []
_render_terminal_sequencing(
    SimpleNamespace(TITLE="mixed read plans"), mixed_graph,
    [(["known"], SequencingPlan(library=known_insert_lib, primers=(sp.TRUSEQ["R1"],),
                              acquired_roles=("Read 1",), read_lengths={"Read 1":100})),
     (["unknown"], SequencingPlan(limitation="Branch-specific primer identity unresolved"))],
    mixed_parts, lambda slot: None)
mixed_html = "\n".join(mixed_parts)
check("mixed exact/unresolved plans retain both terminal sections and the local limitation",
      ("Final sequencing library: represented branch" in mixed_html,
       "Final sequencing library: unresolved branch" in mixed_html,
       "Branch-specific primer identity unresolved" in mixed_html), (True, True, True))

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
annotated = _complete_graph()
annotated.annotate(Highlight("tagment", "Branch", "Retained < discarded.",
                             slot=AnnotationSlot.AFTER_WORKFLOW))
annotated.annotate(SourceDiscrepancy(
    "adapter", "Conflicting sources",
    (SourceClaim("source A, table 1", "AC"),
     SourceClaim("source B, table 2", "ACAT")),
    "ACAT", "source B prints the full oligo", "AT is absent from source A"))
check("graph annotations render at their declared slots and escape source text",
      ("Retained &lt; discarded." in render_annotations(
          annotated, AnnotationSlot.AFTER_WORKFLOW),
       "source A, table 1" in render_annotations(annotated, AnnotationSlot.INTRO),
       "Conflicting sources" in render_annotations(
           annotated, AnnotationSlot.AFTER_WORKFLOW)),
      (True, True, False))
check.raises("annotation target must exist in the graph",
             lambda: annotated.annotate(Highlight("missing", "Bad", "Unknown node")))
check.raises("detail view needs a molecular state",
             lambda: annotated.annotate(DetailView("waste", "Waste", "No drawing")))
check.raises("source discrepancy requires two located claims",
             lambda: annotated.annotate(SourceDiscrepancy(
                 "adapter", "Incomplete", (SourceClaim("source A", "AC"),),
                 "AC", "printed", "uncertain")))
check.raises("page assembly rejects an omitted graph annotation slot",
             lambda: assemble_page(annotated, ["<h1>Example</h1>",
                                               AnnotationSlot.INTRO]))
check("page assembly includes each requested graph annotation slot once",
      assemble_page(annotated, [AnnotationSlot.INTRO,
                                AnnotationSlot.AFTER_WORKFLOW]).count(
                                    'data-graph-target="'), 2)
view_page = assemble_page(annotated, [
    AnnotationSlot.INTRO,
    OligoSection("Source oligos", (OligoEntry("adapter", (Segment("end", "AC"),)),)),
    WorkflowSection("Reaction workflow"), AnnotationSlot.AFTER_WORKFLOW])
check("declared oligo and workflow sections render from one graph",
      ("Source oligos" in view_page and "AC" in view_page,
       "Tagment" in view_page and "Reaction workflow" in view_page),
      (True, True))
branch_graph = _complete_graph()
branch_text = protocol_graph_panels(branch_graph)
check("an unmodelled discarded output is still listed with its limitation",
      ('data-graph-target="waste" data-graph-disposition="discard"' in branch_text,
       'Molecular structure is not represented in this graph.' in branch_text),
      (True, True))
branch_graph.materials["waste"] = replace(
    branch_graph.materials["waste"],
    state=Construct([Segment("off-product", "ACGT")]))
branch_drawing = protocol_graph_panels(branch_graph)
branch_section = branch_drawing.split(
    'data-graph-target="waste" data-graph-disposition="discard"', 1)[1].split(
        '</section>', 1)[0]
check("a modelled discarded output is drawn directly from its graph state",
      ('<svg' in branch_section, 'ACGT' in branch_section), (True, True))
check.raises("an empty declared oligo section is rejected",
             lambda: OligoSection("Empty", ()).render())
check.raises("exact-primer view rejects a graph without exact sequencing geometry",
             lambda: PrimerBindingSection("library", "Reads", "Primers").render(
                 _complete_graph()))
check.raises("source run rows reject a graph without exact sequencing geometry",
             lambda: ReadLayoutSection("library", (
                 SourceRead("Read 1", "50", "insert", "Read 1", "include"),),
                 "source run").render(_complete_graph()))
check.raises("offset primers require a declared sequencing-output limitation",
             lambda: OffsetPrimerSection("library", (("N", "ACGT"),),
                 "offsets", "product", "limitation").render(_complete_graph()))
check("tagmentation permits one source-declared product class without inventing waste", ProtocolGraph(
    "single-product class", authority="source")._validate_step_shape(Step(
        "tagment", "Tagment", StepKind.TAGMENTATION,
        (Input("sample", InputRole.TEMPLATE), Input("adapter", InputRole.STRUCTURAL),
         Input("tn5", InputRole.CATALYST)), ("library",), source="source")), None)
check("binding permits a declared structural complex without inventing a separation", ProtocolGraph(
    "bound complex", authority="source")._validate_step_shape(Step(
        "bind", "Bind SSB to DNA", StepKind.BINDING,
        (Input("dna", InputRole.TEMPLATE), Input("ssb", InputRole.STRUCTURAL)),
        ("complex",), source="source")), None)
check.raises("binding still requires a declared partner", lambda: ProtocolGraph(
    "missing partner", authority="source")._validate_step_shape(Step(
        "bind", "Bind", StepKind.BINDING,
        (Input("dna", InputRole.TEMPLATE),), ("complex",), source="source")))


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

from end_combinations import (EndCombinationGeometry, EndProductGroup,
                              LoadedEnd, LoadedEndPool)  # noqa: E402
import nextera as nx  # noqa: E402
from pcr import PCRGeometry, PCRPrimer  # noqa: E402

end_graph = ProtocolGraph("two-ended outcomes")
end_sample = end_graph.input("sample", "DNA", MaterialKind.MOLECULAR, primary=True)
end_pool = end_graph.input("end_pool", "loaded adapters", MaterialKind.ADAPTER,
                           chemical_state=LoadedEndPool(
                               (LoadedEnd("s5", nx.ADAPTOR_S5),
                                LoadedEnd("s7", nx.ADAPTOR_S7)), nx.ME))
end_enzyme = end_graph.input("enzyme", "Tn5", MaterialKind.CATALYST)
end_kept = end_graph.product("kept", "dual-site fragments", MaterialKind.MOLECULAR)
end_waste = end_graph.product("waste", "same-end fragments", MaterialKind.WASTE)
end_primers = (PCRPrimer("N5xx", nx.n5xx_primer("N" * 8), len(nx.S5)),
               PCRPrimer("N7xx", nx.n7xx_primer("N" * 8), len(nx.S7)))
end_geometry = EndCombinationGeometry(
    end_primers,
    (EndProductGroup("kept", (("s5", "s7"),)),
     EndProductGroup("waste", (("s5", "s5"), ("s7", "s7")))),
    "index_pcr")
end_graph.react("attach", "Attach ends", StepKind.TAGMENTATION,
                ((end_sample, InputRole.TEMPLATE), (end_pool, InputRole.STRUCTURAL),
                 (end_enzyme, InputRole.CATALYST)), (end_kept, end_waste),
                geometry=end_geometry,
                relations=(MolecularRelation.END_COMBINATIONS,))
end_pcr_template = Construct([
    Segment("s5", nx.S5), Segment("me", nx.ME),
    Segment("body", "X" * 30, placeholder=True),
    Segment("me2", nx.ME_RC), Segment("s7", nx.S7_RC)])
end_f = end_graph.input("forward", "N5xx", MaterialKind.OLIGO)
end_r = end_graph.input("reverse", "N7xx", MaterialKind.OLIGO)
end_pol = end_graph.input("polymerase", "PCR polymerase", MaterialKind.CATALYST)
end_library = end_graph.product("library", "library", MaterialKind.LIBRARY)
end_graph.react("index_pcr", "Index PCR", StepKind.PCR,
                ((end_kept, InputRole.TEMPLATE), (end_f, InputRole.STRUCTURAL),
                 (end_r, InputRole.STRUCTURAL), (end_pol, InputRole.CATALYST)),
                (end_library,), geometry=PCRGeometry.exact(
                    end_pcr_template, *end_primers, statement="synthetic test"))
check("all three Tn5 end classes follow from two loaded ends and primer sites",
      [c.status for c in check_chemistry(end_graph) if c.target == "attach"],
      [Status.ESTABLISHED])
end_graph.steps["attach"] = replace(
    end_graph.steps["attach"], geometry=EndCombinationGeometry(
        end_primers,
        (EndProductGroup("kept", (("s5", "s5"),)),
         EndProductGroup("waste", (("s5", "s7"), ("s7", "s7")))),
        "index_pcr"))
check("a same-end class claimed as dual-primer eligible is contradictory",
      [c.status for c in check_chemistry(end_graph) if c.target == "attach"],
      [Status.CONTRADICTORY])
end_graph.steps["attach"] = replace(
    end_graph.steps["attach"], geometry=EndCombinationGeometry(
        end_primers,
        (EndProductGroup("kept", (("s5", "s7"),)),
         EndProductGroup("waste", (("s5", "s5"),))),
        "index_pcr"))
check("omitting an off-product class is contradictory",
      [c.status for c in check_chemistry(end_graph) if c.target == "attach"],
      [Status.CONTRADICTORY])
end_graph.steps["attach"] = replace(end_graph.steps["attach"], geometry=end_geometry)
end_graph.materials["end_pool"] = replace(
    end_pool, chemical_state=LoadedEndPool(
        (LoadedEnd("s5", nx.ADAPTOR_S5),
         LoadedEnd("s7", nx.ADAPTOR_S5)), nx.ME))
check("a loaded end missing its asserted primer footprint is contradictory",
      [c.status for c in check_chemistry(end_graph) if c.target == "attach"],
      [Status.CONTRADICTORY])
end_graph.materials["end_pool"] = end_pool
altered_primer = replace(end_primers[0],
                         sequence="C" + end_primers[0].sequence[1:])
end_graph.steps["index_pcr"] = replace(
    end_graph.steps["index_pcr"], geometry=PCRGeometry.exact(
        end_pcr_template, altered_primer, end_primers[1],
        statement="synthetic changed downstream primer"))
check("end classes must use the primers of the actual downstream PCR step",
      [c.status for c in check_chemistry(end_graph) if c.target == "attach"],
      [Status.CONTRADICTORY])

check.section("molecular diagram renderers")
svg_panel = panel([Row(chunks=[("ACGT", "r1", False)])])
check("SVG is the default diagram renderer", '<svg class="chem-svg' in svg_panel)
check("SVG keeps sequence as selectable text", "<tspan class=\"chem-r1\">ACGT</tspan>" in svg_panel)
check("the legacy renderer remains explicitly available",
      panel([Row(chunks=[("ACGT", None, False)])], renderer="legacy").startswith("<pre>"))
text_oligo = oligo("PCR primer", [Segment("binding", "ACGTACGT")])
check("ordering oligos stay raw selectable text", "<svg" not in text_oligo)
check("ordering oligos retain explicit strand ends", "5'- ACGTACGT -3'" in text_oligo)

from types import SimpleNamespace
from multimodal_page import render as render_multimodal
source_readout=SimpleNamespace(TITLE="Source-declared endpoint",NOTES="notes.html",
    SOURCE="Original source",SUMMARY="One library",sections=lambda:[],
    FINAL_LIBRARIES=(("Library",Construct([Segment("insert","ACGTACGT")]),(),"Library","Readout"),),
    SEQUENCING_UNAVAILABLE={"Library":"Loaded primer sequence remains unverified."},
    READ_SPECS={"Library":(sp.read_spec("Read 1",8,"source-declared payload","include"),)})
source_readout_html=render_multimodal(source_readout)
check("multimodal source reads render without invented primer sequences",
      "source-declared payload" in source_readout_html and "Loaded primer sequence remains unverified." in source_readout_html)
check.raises("missing primer geometry requires an explicit source read layout",
    lambda:render_multimodal(SimpleNamespace(**{**vars(source_readout),"READ_SPECS":{}})))
check.raises("source read layout cannot hide an unexplained missing primer",
    lambda:render_multimodal(SimpleNamespace(**{**vars(source_readout),"SEQUENCING_UNAVAILABLE":{}})))


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
check("tables are wrapped so the page never scrolls sideways",
      'class="tw"' in build_docs.build_page(ROOT / "catalogue" / "README.md"))
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
# Many protocols name their module `chem`; a note must reach its own, never another's.
_chem_dirs = sorted(p.parent.parent for p in ROOT.glob("*/tools/chem.py"))[:2]
check("a note's bare module name resolves to its own protocol's tools/ module",
      [expand("{{= chem.__file__ }}", f"{d.name}/01_note.md") for d in _chem_dirs],
      [str(d / "tools" / "chem.py") for d in _chem_dirs])
try:
    expand("{{= chem.ANY_NAME }}", "README.md")
    _amb = None
except ValueError as e:
    _amb = str(e)
check("a protocol module name shared by several directories is an error outside them",
      "ambiguous module name 'chem'" in (_amb or ""), True)

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

# A filename prefix does not establish provenance: source bundles sometimes contain
# an original publisher supplement under upstream_* as well as an upstream drawing.
from contextlib import redirect_stdout
from io import StringIO
from tempfile import TemporaryDirectory

with TemporaryDirectory() as temp:
    _bundle = Path(temp)
    _bundle.joinpath("MANIFEST.tsv").write_text(
        "file\turl\tbytes\tmd5\twhat\n"
        "upstream_schematic.html\thttps://example.org/drawn\t\t\tscg_lib_structs method page\n"
        "upstream_author_supplement.xlsx\t\t\t\t(on disk, not fetched by get_sources: added by hand?)\n"
        "upstream_vendor_original.pdf\thttps://example.org/vendor\t\t\tvendor original oligo table\n"
        "vendor_original.pdf\thttps://example.org/guide\t\t\tvendor original oligo table\n")
    for _name, _seq in (
        ("upstream_schematic.html.txt", "AATGATACGGCGACCACCGAGATCTACAC"),
        ("upstream_author_supplement.xlsx.txt", "CAAGCAGAAGACGGCATACGAGATGTGACT"),
        ("upstream_vendor_original.pdf.txt", "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"),
        ("vendor_original.pdf.txt", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"),
    ):
        _bundle.joinpath(_name).write_text("Primer\t" + _seq + "\n")
    _output = StringIO()
    with redirect_stdout(_output):
        spx.main([str(_bundle), "--tsv"])
    _lines = [line.split("\t") for line in _output.getvalue().splitlines()]
    check("scraper labels each hit from manifest and ranks the drawing last",
          {Path(row[0]).name: row[-2] for row in _lines[1:]},
          {"vendor_original.pdf.txt": "source document",
           "upstream_vendor_original.pdf.txt": "source document",
           "upstream_author_supplement.xlsx.txt": "unclassified — verify original",
           "upstream_schematic.html.txt": "upstream schematic"})
    check("original supplement under upstream_ retains its manifest description",
          next(row[-1] for row in _lines[1:]
               if Path(row[0]).name == "upstream_author_supplement.xlsx.txt"),
          "(on disk, not fetched by get_sources: added by hand?)")

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

# NEB E7600 v6.0_6/24 p.3 prints unbracketed CdUA inside a codon-spaced oligo.
# Keep the complete molecule and its marks, rather than two hits around dU.
_neb_du = ("NEBNext Adaptor\n"
           "5´-/5Phos/GAT CGG AAG AGC ACA CGT CTG AAC TCC AGT CdUA CAC TCT TTC "
           "CCT ACA CGA CGC TCT TCC GAT C-s-T-3´\n")
_neb_du_hits = spx.scan(Path("original_neb_oligo"), _neb_du, {})
check("unbracketed internal dU keeps the complete source-printed hairpin",
      [h.seq for h in _neb_du_hits],
      ["GATCGGAAGAGCACACGTCTGAACTCCAGTCUACACTCTTTCCCTACACGACGCTCTTCCGATCT"])
check("internal dU preserves written phosphate, deoxyuridine and phosphorothioate",
      [h.written for h in _neb_du_hits],
      ["/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCdUACACTCTTTCCCTACACGACGCTCTTCCGATC-s-T"])
check("explicit dU is recorded as deoxyuridine rather than RNA",
      [h.mods for h in _neb_du_hits],
      [["5' phosphate (/5Phos/)", "phosphorothioate x1", "deoxyuridine x1"]])
check("--find crosses the internal dU and retains the original source line",
      spx.find(_neb_du, _neb_du_hits[0].seq), [2])
check("terminal and consecutive explicit dU keep their nucleotide positions",
      [(h.seq, h.mods) for h in spx.scan(Path("du_ends"),
          "start\tdUdUACGTACGTACGT\nend\tACGTACGTACGTdU\n", {})],
      [("UUACGTACGTACGT", ["deoxyuridine x2"]),
       ("ACGTACGTACGTU", ["deoxyuridine x1"])])
check("an isolated prose dU mention is not rewritten as a sequence",
      spx.scan_copy("the dU lesion").text, "the dU lesion")

_bracketed = ("ME-ddC\t[Phos]CTGTCTCTTATACACATCT[ddc]\n"
              "Blocker\tGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[invdT]\n")
_bracketed_hits = {h.label: h for h in spx.scan(Path("sum_oligos"), _bracketed,
                                                spx.known_sequences())}
check("bracketed terminal blocks survive source scraping",
      [(_bracketed_hits[n].written, _bracketed_hits[n].mods) for n in ("ME-ddC", "Blocker")],
      [("[Phos]CTGTCTCTTATACACATCT[ddc]",
        ["5' phosphate ([Phos])", "3' dideoxy ([ddc])"]),
       ("GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[invdT]",
        ["3' inverted dT ([invdT])"])])

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
       "TTTCTTATATGGGcDNA (101 nt) and T30 oligo in prose\n"
       "RevBC 5'-CAAGCAGAAGACGGCATACGAGAT-barcode-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-3'\n"
       "UMIhy ACACTCTTTCCCTACACGACGCTCTTCCGATCT-UMI-TTTTTTTTTTTTTTTTTTTTVN\n"
       "Neg GTGACTGGAGTTCAGACGTGTGC-barcoded-CAAGCAGAAGACGGCATACGAGAT\n")
_h2 = {h.label: h for h in spx.scan(Path("doc"), _d2, spx.known_sequences())}
check("a hyphen-delimited bare placeholder word (-barcode-, -UMI-) is an uncounted gap",
      [(_h2[k].seq.count(spx.GAP), _h2[k].length_text) for k in ("RevBC", "UMIhy")],
      [(1, "57+ nt"), (1, "55+ nt")])
check("an ordinary hyphenated word ('-barcoded-') is not a placeholder",
      (_h2["Neg"].seq, spx.GAP in _h2["Neg"].seq), ("GTGACTGGAGTTCAGACGTGTGC", False))
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
check("--find: bracketed biotin-dT preserves its base and original notation",
      [(strand, exact, written) for _, strand, exact, written in spx.find(
          "A1\tGATCGAGCTCGAGAA[BtndT]T\n", "GATCGAGCTCGAGAATT", detail=True)],
      [("+", "exact", "GATCGAGCTCGAGAA[BtndT]T")])
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
check("terminal biotin-dT and blocked ddC retain all three nucleotide positions",
      spx._letters("/5App/CTGACNNNNNNNNNNNNNNNTGGAATTCTCGGGTGCCAAGGC/iBiodT//iBiodT//3ddC/"),
      "CTGAC" + "N"*15 + "TGGAATTCTCGGGTGCCAAGGCTTC")
check("source guard accepts a printed 3'-to-5' strand as an ordered 5'-to-3' input",
      _exact_hit("P2 lower,3'-TCT CTT ACT CCT TGG GCC CCG TC-5'\n",
                 "CTGCCCCGGGTTCCTCATTCTCT"), True)

# Source forms found by curation rounds 7-20 (ChAIR, easySHARE, FIND-seq, NEBNext E6420,
# OK-seq, UMI-4C, scRibo-seq, RAMPAGE, STRT-seq, sci-RNA-seq).
check("internal modified nucleotides keep their base; spacers and dyes add none",
      [spx.internal_base(x) for x in ("iBIOdT", "iAmMC6T", "iFluorT", "i5MedC", "iMe-dC",
                                      "iBiodG", "ideoxyU", "ideoxyI", "iSp18", "iSpC3",
                                      "iSpPC", "iCy3", "iBHQ1")],
      ["T", "T", "T", "C", "C", "G", "U", "N", "", "", "", "", ""])
check("an internal token after a 5' /mod/ is not hidden inside a slash-delimited base run",
      spx.scan_copy("5'-/5Phos/CGTGATATT/iBIOdT/CACGACTCT-3'").text,
      "5'-/5Phos/CGTGATATTTCACGACTCT-3'")
check("source guard reads internal biotin-dT and amino-dT as their T",
      [_exact_hit("x\t/5Phos/GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNNT/iBiodT/TTTTTTTTTTTTTTTTTTTTTTTTVN\n",
                  "GGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNNN" + "T" * 26 + "VN"),
       _exact_hit("bridge CGTGATATT/iAmMC6T/CACGACTCT here", "CGTGATATTTCACGACTCT")],
      [True, True])
check("a quencher glued to the 3' base is a modification; the probe bases are kept",
      [(h.seq, h.mods) for h in spx.scan(Path("find"),
          "probe: 5′-FAM-TTTTGGCGTACTCACCAGTMGBNHQ-3′; 6,000\n", {})],
      [("TTTTGGCGTACTCACCAGT", ["5' FAM", "3' MGBNHQ"])])
check("codon triplets wrapped over a line stay one oligo up to the 3' end",
      [h.seq for h in spx.scan(Path("neb"), "5´-AAG CAG TGG TAT CAA CGC AGA GTA CTT TTT TTT "
                               "TTT TTT TTT TTT TTT\nTTT TTT TV-3´\n", {})],
      ["AAGCAGTGGTATCAACGCAGAGTAC" + "T" * 30 + "V"])
check("a run of asterisks between bases is an index placeholder, not phosphorothioates",
      [(h.seq, h.mods) for h in spx.scan(Path("okseq"), "5′-CAAGCAGAAGACGGCATACGAGAT-3′******"
                                         "5′-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′\n", {})],
      [("CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT", [])])
check("(NNNNNN), IDT r(N1:...) mixes, base(N)15 and XXXXXXrGrGrG are variable windows",
      [_exact_hit("(PO4)-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC(NNNNNN)ATCTCGTATG\n",
                  "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATG"),
       _exact_hit("X\trCrUrArCrArCrGrArCrGrCrUrCrUrUrCrCrGrArUrCrUr(N1:25252525)r(N1)r(N1)\n",
                  "CUACACGACGCUCUUCCGAUCUNNN"),
       _exact_hit("RT: 5 -TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCT(N)15\n",
                  "TAGTCGAACGAAGGTCTCCGAACCGCTCTTCCGATCT" + "N" * 15),
       _exact_hit("TSO 5'-AAGCAGTGGTATCAACGCAGAGTGCAGTGCTXXXXXXrGrGrG-3'\n",
                  "AAGCAGTGGTATCAACGCAGAGTGCAGTGCTNNNNNNGGG")],
      [True, True, True, True])
check("a printed [i7] of unstated length stands only for a model N window",
      [_exact_hit("P7 CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG\n",
                  "CAAGCAGAAGACGGCATACGAGAT" + "N" * 10 + "GTCTCGTGGGCTCGG"),
       _exact_hit("P7 CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG\n",
                  "CAAGCAGAAGACGGCATACGAGATNNNNACGTNNGTCTCGTGGGCTCGG")],
      [True, False])

check("a 3' run of slash-delimited modified bases reads as bases; A/G stays an alternative",
      [spx.internal_base("iXNA_G"), spx.internal_base("Index_A"),
       _exact_hit("TSO LNA\t AAGCAGTGGTATCAACGCAGAGTGAAT/rG/rG/iXNA_G \t\n",
                  "AAGCAGTGGTATCAACGCAGAGTGAATGGG"),
       _exact_hit("x CTGCTGACGTACTGAGAGGC/rG//rG//iXNA_G/ y\n", "CTGCTGACGTACTGAGAGGCGGG"),
       spx.scan_copy("x GCTAGCTAGCTAGCTA/G y").text],
      ["G", "", True, True, "x GCTAGCTAGCTAGCTA/G y"])
import zipfile  # noqa: E402
with TemporaryDirectory() as _xtmp:
    _xl = Path(_xtmp) / "rows.xlsx"
    with zipfile.ZipFile(_xl, "w") as _z:
        _z.writestr("xl/worksheets/sheet1.xml",
                    '<worksheet><sheetData><row r="1"/><row r="2"><c t="inlineStr"><is><t>A</t>'
                    '</is></c></row><row r="3"><c t="inlineStr"><is><t>B</t></is></c></row>'
                    '</sheetData></worksheet>')
    check("xlsx: a self-closing empty <row/> does not swallow the next row",
          spx.xlsx_text(_xl).splitlines()[1:], ["A", "B"])

from chemdraw import Construct as _GC, Scene as _GS, Segment as _GSeg, feature as _gfeat  # noqa: E402
from source_oligo_check import _sequences, manifest_locators, source_texts  # noqa: E402
import source_oligo_check as source_guard  # noqa: E402
_bc = _gfeat("probe_bc", "sample_index", "unknown")
_scenes = []
for _top in ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"):
    _s = _GS()
    _s.strand("top", [_GSeg("adapter", _top)])
    _scenes.append(_s)
check("a tuple of loaded-adapter Scenes yields every strand of every Scene",
      _sequences(tuple(_scenes)),
      ("TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"))
check("identifier and IUPAC-class placeholders are variable; other stand-ins are not",
      [_sequences(_GC([_GSeg("P7", "CAAGCAGAAGACGGCATACGAGAT"),
                       _GSeg("i7", "XXXXXXXX", placeholder=True, feature=_bc),
                       _GSeg("arm", "GTCTCGTGGGCTCGG")])),
       _sequences(_GC([_GSeg("dT", "T" * 20), _GSeg("end", "VN", placeholder=True)])),
       _sequences(_GC([_GSeg("handle", "AAGCAGTGGTATCAACGCAGAGT"),
                       _GSeg("insert", "NNNNNNNNNN", placeholder=True)])),
       _sequences(_GC([_GSeg("handle", "AAGCAGTGGTATCAACGCAGAGT"),
                       _GSeg("guess", "ACGT", inferred=True)]))],
      [("CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGG",), ("T" * 20 + "NN",), (), ()])
check("a placeholder in source modification notation contributes its nucleotides only",
      [_sequences(_GC([_GSeg("arm", "CGTCAGACAGGGAGCG"), _GSeg("C3", "[C3]", placeholder=True),
                       _GSeg("arm 2", "GGAAGTCGCC"), _GSeg("end", "[3SpC3]", placeholder=True)])),
       _sequences(_GC([_GSeg("handle", "CTGCTGACGTACTGAGAGGC"),
                       _GSeg("ends", "/rG//rG//iXNA_G/", placeholder=True, length_bp=3)])),
       _sequences(_GC([_GSeg("handle", "CTGCTGACGTACTGAGAGGC"),
                       _GSeg("ends", "/rG//rG//iXNA_G/", placeholder=True, length_bp=2)])),
       _sequences(_GC([_GSeg("handle", "CTGCTGACGTACTGAGAGGC"),
                       _GSeg("bead", "[bead]", placeholder=True)]))],
      [("CGTCAGACAGGGAGCGGGAAGTCGCC",), ("CTGCTGACGTACTGAGAGGCGGG",), (), ()])
with TemporaryDirectory() as temp:
    _bundle = Path(temp) / "proto"
    _bundle.mkdir()
    _bundle.joinpath("MANIFEST.tsv").write_text(
        "file\turl\tbytes\tmd5\twhat\n"
        "preprint_media-2.xlsx\thttps://biorxiv.example/media-2.xlsx\t\t\tbiorxiv supplementary file\n"
        "(duplicate) s99999-000-0001-x_MOESM2_ESM.xlsx\thttps://springer.example/MOESM2.xlsx"
        "\t\t\tSpringer supplementary file; identical to preprint_media-2.xlsx\n"
        "vendor_primers.fasta\thttps://vendor.example/primers.fasta\t\t\tvendor primer FASTA\n")
    _bundle.joinpath("preprint_media-2.xlsx.txt").write_text("P1\tACGTACGTACGTACGT\n")
    _bundle.joinpath("vendor_primers.fasta").write_text(">fwd\nACGTTGCAACGTTGCA\n")
    _bundle.joinpath("scg_derived.fa").write_text(">x\nACGTTGCAACGTTGCA\n")
    _bundle.joinpath("MANIFEST.tsv").write_text(
        _bundle.joinpath("MANIFEST.tsv").read_text()
        + "repository-full.txt\thttps://repo.example/s99999-000-0001-x/full.txt\t\t\tdeposited text\n")
    _bundle.joinpath("repository-full.txt").write_text("RT\tACGTTGCAACGTTGCAAC\n")
    check("a (duplicate) manifest row lends its published name and URL to the on-disk twin",
          "s99999-000-0001-x_moesm2_esm.xlsx" in " ".join(
              manifest_locators(_bundle / "MANIFEST.tsv")["preprint_media-2.xlsx"]).lower(), True)
    check("the guard reads a published supplement's twin and a plain-text original by its own manifest URL",
          sorted(p.name for p, _ in source_texts(_bundle, "10.1038/s99999-000-0001-x")),
          ["preprint_media-2.xlsx.txt", "repository-full.txt"])
    check("a declared FASTA original is read; an scg_ FASTA is not",
          sorted(p.name for p, _ in source_texts(_bundle, "")),
          ["preprint_media-2.xlsx.txt", "repository-full.txt", "vendor_primers.fasta"])
    from unittest.mock import patch
    _catalogue = Path(temp) / "catalogue"
    _catalogue.mkdir()
    _catalogue.joinpath("scg_lib_structs.tsv").write_text(
        "our_dir\trole\tdoi\n"
        "proto\tprimary\t10.1101/123456\n"
        "proto\tprimary\t10.1038/s99999-000-0001-x\n")
    check("Genome Research DOIs are published articles, not bioRxiv identifiers",
          source_guard._biorxiv_doi("10.1101/gr.161679.113"), False)
    _bundle.joinpath("123456_v1.full.pdf.txt").write_text("Preprint oligo\tACGTTGCAACGTTGCA\n")
    with patch.object(source_guard, "ROOT", Path(temp)):
        check("linked preprints do not validate a published protocol's oligos",
              sorted(p.name for p, _ in source_texts(_bundle, "10.1038/s99999-000-0001-x")),
              ["preprint_media-2.xlsx.txt", "repository-full.txt"])
        check("a preprint-named directory uses its linked publication as authority",
              sorted(p.name for p, _ in source_texts(_bundle, "10.1101/123456")),
              ["preprint_media-2.xlsx.txt", "repository-full.txt"])
    _big = Path(temp) / "big"
    _big.mkdir()
    _big.joinpath("table.txt").write_text("P1\tACGTACGTACGTACGTAC\n")
    _big.joinpath("whitelist.csv").write_text("x" * 64)
    _saved = spx.LARGE_BYTES
    spx.LARGE_BYTES = 32
    _err, _out = StringIO(), StringIO()
    from contextlib import redirect_stderr
    with redirect_stdout(_out), redirect_stderr(_err):
        spx.main([str(_big), "--tsv"])
    spx.LARGE_BYTES = _saved
    check("scan mode skips an over-size file and names it rather than ignoring it",
          ("whitelist.csv" in _err.getvalue().splitlines()[-1], "table.txt" in _out.getvalue()),
          (True, True))

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
