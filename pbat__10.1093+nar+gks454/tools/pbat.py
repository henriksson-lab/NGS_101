"""Miura et al. amplification-free PBAT model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Segment
from chromatin_epigenetics import pbat_priming_scene
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "PBAT — post-bisulfite adaptor tagging"
NOTES = "01_pbat.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gks454">Miura et al., <i>Nucleic Acids Research</i> (2012)</a>.'
SUMMARY = "Bisulfite-convert genomic DNA first, then install the two sequencing handles by successive random-primer extensions, avoiding conversion-induced loss of completed library molecules."
BIOPEA2 = il.TRUSEQ_READ1
PE_REVERSE = il.P7
PRIMER3 = il.P5 + il.TRUSEQ_READ1[4:]
FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"), Segment("Read 1 arm remainder", BIOPEA2[4:], "r1"),
    Segment("first random N4", "N" * 4, placeholder=True),
    Segment("directional PBAT insert", "X" * 36, placeholder=True),
    Segment("second random N4 reverse complement", "N" * 4, placeholder=True),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
], name="amplification-free PBAT sequencing template")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("invalid PBAT read layout: " + "; ".join(problems))
FINAL_CAPTION = "Amplification-free PBAT template after Primer-3 extension. The two N4 tracts record the random priming sites; the original protocol has no sample index."
SEQUENCING_INTRO = "The original GAIIx/HiSeq protocol produces a single Read 1 from the strand complementary to bisulfite-converted DNA; reads are consequently generally G-poor."

def sections():
    return [
        ("Bisulfite-convert untagged genomic DNA", [Row(chunks=[("genomic duplex → denature + bisulfite → fragmented C→U single strands", "w1", False)])],
         "Conversion precedes adapter tagging, the defining reversal relative to conventional WGBS."),
        ("First random-primer extension and capture", pbat_priming_scene("first-strand", "BioPEA2", BIOPEA2, mod5="biotin").rows(),
         "A 5′-biotinylated adapter primer ending in random N4 copies the converted strand; streptavidin immobilizes the product."),
        ("Second random-primer extension", pbat_priming_scene("second-strand", "PE-reverse", PE_REVERSE).rows(),
         "After alkaline denaturation, PE-reverse-N4 copies the bead-bound first strand; Primer-3 then adds P5 and completes the sequencing template."),
    ]


def protocol_graph() -> ProtocolGraph:
    """Amplification-free PBAT graph with explicit capture and strand fates."""
    source = "Miura et al. 2012, doi:10.1093/nar/gks454"
    graph = ProtocolGraph(TITLE, authority=source)

    def add_input(key, label, kind, *, state=None, primary=False):
        value = Material(key, label, kind, state=state, primary=primary, source=source)
        graph.add_input(value)
        return value

    dna = add_input("genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
                    state=(Row(chunks=[("genomic DNA duplex", None, False)]),), primary=True)
    bisulfite = add_input("bisulfite_reagent", "bisulfite conversion reagent",
                          MaterialKind.CONDITION)
    converted = Material(
        "converted_ssdna", "fragmented bisulfite-converted single-stranded DNA",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[("fragmented C→U single-stranded DNA", "w1", False)]),),
        source=source)
    graph.add_step(Step(
        "bisulfite_conversion", "Bisulfite-convert untagged genomic DNA",
        StepKind.CHEMICAL_CONVERSION,
        (Input(dna.id, InputRole.TEMPLATE), Input(bisulfite.id, InputRole.CONDITION)),
        (converted.id,), source=source), (converted,))

    biopea2 = add_input("biopea2_n4", "5-prime-biotinylated BioPEA2-N4 primer",
                        MaterialKind.OLIGO)
    polymerase = add_input("strand_synthesis_polymerase", "strand-synthesis polymerase",
                           MaterialKind.CATALYST)
    first_product = Material(
        "first_extension_product", "biotinylated first random-primer extension product",
        MaterialKind.MOLECULAR,
        state=pbat_priming_scene("first-strand", "BioPEA2", BIOPEA2, mod5="biotin"),
        source=source)
    graph.add_step(Step(
        "first_extension", "First random-primer extension", StepKind.EXTENSION,
        (Input(converted.id, InputRole.TEMPLATE),
         Input(biopea2.id, InputRole.STRUCTURAL),
         Input(polymerase.id, InputRole.CATALYST)),
        (first_product.id,), source=source), (first_product,))

    streptavidin = add_input("streptavidin_beads", "streptavidin magnetic beads",
                             MaterialKind.SELECTION)
    captured = Material("captured_first_product", "bead-captured biotinylated product",
                        MaterialKind.MOLECULAR, state=first_product.state, source=source)
    uncaptured = Material("uncaptured_molecules", "uncaptured first-round molecules",
                          MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "streptavidin_capture", "Capture the biotinylated first extension product",
        StepKind.AFFINITY_CAPTURE,
        (Input(first_product.id, InputRole.TEMPLATE),
         Input(streptavidin.id, InputRole.SELECTION)),
        (captured.id, uncaptured.id), source=source), (captured, uncaptured))
    graph.mark_terminal(Terminal(
        uncaptured.id, Disposition.DISCARD,
        "not retained by the streptavidin capture"))

    alkali = add_input("alkaline_denaturation", "alkaline denaturation condition",
                       MaterialKind.CONDITION)
    bead_strand = Material(
        "bead_bound_first_strand", "bead-bound first extension strand",
        MaterialKind.MOLECULAR,
        state=(Row(chunks=[("5′-biotin—BioPEA2—N4—copied converted strand", None, False)]),),
        source=source)
    released_parent = Material(
        "released_parental_strand", "released bisulfite-converted parental strand",
        MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "strand_separation", "Alkaline-denature the captured duplex",
        StepKind.STRAND_SEPARATION,
        (Input(captured.id, InputRole.TEMPLATE), Input(alkali.id, InputRole.CONDITION)),
        (bead_strand.id, released_parent.id), source=source),
        (bead_strand, released_parent))
    graph.mark_terminal(Terminal(
        released_parent.id, Disposition.DISCARD,
        "released by alkaline denaturation while the biotinylated copy remains bead-bound"))

    pe_reverse = add_input("pe_reverse_n4", "PE-reverse-N4 second random primer",
                           MaterialKind.OLIGO)
    second_product = Material(
        "second_extension_product", "second random-primer extension product",
        MaterialKind.MOLECULAR,
        state=pbat_priming_scene("second-strand", "PE-reverse", PE_REVERSE),
        source=source)
    graph.add_step(Step(
        "second_extension", "Second random-primer extension", StepKind.EXTENSION,
        (Input(bead_strand.id, InputRole.TEMPLATE),
         Input(pe_reverse.id, InputRole.STRUCTURAL),
         Input(polymerase.id, InputRole.CATALYST)),
        (second_product.id,), source=source), (second_product,))

    primer3 = add_input("primer_3", "Primer-3 P5 completion primer", MaterialKind.OLIGO)
    library = Material("pbat_library", "amplification-free PBAT sequencing template",
                       MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "primer3_extension", "Primer-3 extension completes the sequencing template",
        StepKind.EXTENSION,
        (Input(second_product.id, InputRole.TEMPLATE),
         Input(primer3.id, InputRole.STRUCTURAL),
         Input(polymerase.id, InputRole.CATALYST)),
        (library.id,), source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS)))

    graph.declare_reviewed_scope((
        "genomic_dna", "bisulfite_reagent", "bisulfite_conversion", "biopea2_n4",
        "strand_synthesis_polymerase", "first_extension", "streptavidin_beads",
        "streptavidin_capture", "alkaline_denaturation", "strand_separation",
        "pe_reverse_n4", "second_extension", "primer_3", "primer3_extension",
        "pbat_library",
    ), source)
    return graph.require_publishable()
