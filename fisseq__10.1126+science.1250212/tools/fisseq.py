"""Molecular model, graph, and page for transcriptome-wide RNA FISSEQ."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from chemdraw import Construct, Scene, Segment, annotation_rows
from circular import CircularizationGeometry, circularize_ssdna, rolling_circle
from page import head, info
from protocol_graph import (AnnotationSlot, DetailStyle, DetailView, InputRole, MaterialKind,
                            MolecularRelation, ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_views import (OffsetPrimerSection, PlanLimitationSection, WorkflowSection,
                            assemble_page)
from topology import CircularPrimerTopology


TITLE = "FISSEQ — transcriptome sequencing inside fixed cells"
PAGE_FILE = "fisseq.html"
SOURCE = ("Lee et al., Science (2014), doi:10.1126/science.1250212; "
          "Lee et al., Nature Protocols (2015), doi:10.1038/nprot.2014.191")


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


ADAPTER = "TCTCGGGAACGCTGAAGA"
RANDOM_HEXAMER = "N" * 6
CDNA = "N" * 36
RT_PRIMER_WRITTEN = "/5phos/TCTCGGGAACGCTGAAGANNNNNN"
RCA_PRIMER = "TCTTCAGCGTTCCCGAGA"
RCA_PRIMER_WRITTEN = "TCTTCAGCGTTCCCGA*G*A"
SEQUENCING_PRIMERS = (
    ("N", ADAPTER),
    ("N-1", ADAPTER[1:]),
    ("N-2", ADAPTER[2:]),
    ("N-3", ADAPTER[3:]),
    ("N-4", ADAPTER[4:]),
)
CIRCLIGASE_C = 60
CIRCLIGASE_MIN = 60


def rt_scene() -> Scene:
    target = [seg("RNA upstream", "N" * 18, placeholder=True),
              seg("random landing site", "N" * 6, placeholder=True),
              seg("RNA downstream", "N" * 18, placeholder=True)]
    primer = [seg("5-prime adapter", ADAPTER, "r1"),
              seg("random hexamer", RANDOM_HEXAMER, "r3", placeholder=True)]
    sc = Scene(); sc.strand("RNA", target, label="cellular RNA")
    sc.anneal("RT primer", primer, to="RNA",
              pair=("random hexamer", "random landing site"),
              unpaired=("5-prime adapter",), label="5'-p RT primer")
    sc.arrow("RT primer", "M-MuLV copies RNA; aminoallyl-dUTP marks cDNA")
    return sc


def linear_cdna() -> Construct:
    return Construct([
        seg("adapter", ADAPTER, "r1"),
        seg("random hexamer", RANDOM_HEXAMER, "r3", placeholder=True),
        seg("cDNA", CDNA, placeholder=True),
    ], name="cross-linked single-stranded cDNA")


CIRCLE = circularize_ssdna(linear_cdna(), five_prime_phosphate=True)
RCA = rolling_circle(CIRCLE, RCA_PRIMER, copies=3)
RCA_TOPOLOGY = CircularPrimerTopology(
    "FISSEQ covalently closed cDNA",
    CIRCLE.sequence,
    "phosphorothioate-protected RCA primer",
    RCA_PRIMER,
)


def rca_priming_scene() -> Scene:
    template = list(linear_cdna())
    primer = [seg("RCA primer", RCA.primer, "r2")]
    sc = Scene(); sc.strand("circle opened for drawing", template, label="cDNA circle")
    sc.anneal("RCA primer", primer, to="circle opened for drawing",
              pair=("RCA primer", "adapter"), label="RCA primer")
    sc.arrow("RCA primer", "Phi29 copies repeatedly around the circle")
    return sc


def rca_primed_rows():
    """Show circle topology and an exact opened-out view of the primer site together."""
    return [*CIRCLE.rows("CircLigase II closure"),
            *rca_priming_scene().rows()]


def rolony() -> Construct:
    p, c = len(RCA.primer), len(CDNA)
    unit = RCA.unit
    one = (unit[:p], unit[p:p + c], unit[p + c:])
    parts = []
    for i in range(1, RCA.copies + 1):
        parts.extend([
            seg(f"adapter complement {i}", one[0], "r1"),
            seg(f"transcript copy {i}", one[1], placeholder=True),
            seg(f"random hexamer copy {i}", one[2], "r3", placeholder=True),
        ])
    product = Construct(parts, name="RCA concatemer (three repeats drawn)")
    if product.top() != RCA.sequence:
        raise AssertionError("drawn rolony must come from the RCA product")
    return product


SEQ_N = sp.custom("SOLiD sequencing", "FISSEQ primer N", ADAPTER,
                  "Lee et al., Nature Protocols 2015")


def sequencing_scene() -> Scene:
    product = rolony()
    hits = sp.locate_all(product, SEQ_N)
    if len(hits) != RCA.copies:
        raise ValueError("primer N must bind once per RCA repeat")
    target = hits[len(hits) // 2].covers[0]
    sc = Scene(); sc.strand("rolony", list(product), label="RCA concatemer")
    for name, sequence in SEQUENCING_PRIMERS:
        strand = f"primer {name}"
        segment = f"sequencing primer {name}"
        sc.anneal(strand, [seg(segment, sequence, "r2")], to="rolony",
                  pair=(segment, target), label=strand, mod5="p")
    # Lee et al. 2015 Fig. 1b/Box 1: probe cleavage leaves a 5'-phosphorylated ligation
    # complex, and the offsets are 5' truncations with a shared 3' end, so probes ligate
    # at each primer's 5'-phosphate end (toward the transcript copy), not by 3' extension.
    sc.note("primer N", "probes ligate at the 5'-phosphate ends; the 3' ends coincide")
    return sc


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("cellular_rna", "cellular RNA in fixed cells or tissue",
                      MaterialKind.BIOLOGICAL, state=rt_scene(), primary=True)
    rt_primer = graph.input(
        "rt_primer", "5-prime-phosphorylated adapter–random-hexamer RT primer",
        MaterialKind.OLIGO,
        state=Construct(list(linear_cdna())[:2], name="FISSEQ RT primer"))
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "M-MuLV reverse transcriptase", MaterialKind.CATALYST)
    aminoallyl_dutp = graph.input(
        "aminoallyl_dutp", "aminoallyl-dUTP", MaterialKind.MOLECULAR)
    hybrid = graph.product(
        "rna_cdna_hybrid", "aminoallyl-bearing RNA–cDNA hybrid",
        MaterialKind.MOLECULAR, state=linear_cdna())
    graph.react(
        "reverse_transcribe",
        "Random-prime RT with aminoallyl-dUTP (10 min 4 °C, then 37 °C 1–2 h to overnight)",
        StepKind.REVERSE_TRANSCRIPTION,
        ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
         (reverse_transcriptase, InputRole.CATALYST),
         (aminoallyl_dutp, InputRole.STRUCTURAL)), (hybrid,))

    crosslinker = graph.input("crosslinker", "BS(PEG)9 cross-linker",
                              MaterialKind.CONDITION)
    immobilized_hybrid = graph.product(
        "immobilized_hybrid", "immobilized aminoallyl-bearing RNA–cDNA hybrid",
        MaterialKind.MOLECULAR, state=linear_cdna())
    graph.react("immobilize_cdna",
                "Cross-link the aminoallyl-bearing cDNA in place (BS(PEG)9 1 h, Tris quench)",
                StepKind.CROSS_LINKING,
                ((hybrid, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (immobilized_hybrid,), relations=(MolecularRelation.STRAND_ORDER,))

    ribonucleases = graph.input("ribonucleases", "RNase plus RNase H",
                                MaterialKind.CATALYST)
    linear = graph.product("linear_cdna", "immobilized linear single-stranded cDNA",
                           MaterialKind.MOLECULAR, state=linear_cdna())
    rna_waste = graph.product("rna_waste", "digested RNA", MaterialKind.WASTE)
    graph.react("remove_rna", "Remove the RNA template (RNase + RNase H, 1 h 37 °C)",
                StepKind.DIGESTION,
                ((immobilized_hybrid, InputRole.TEMPLATE),
                 (ribonucleases, InputRole.CATALYST)), (linear, rna_waste))
    graph.discard(rna_waste, "RNase and RNase H digestion products are washed away")

    circligase = graph.input("circligase", "CircLigase II", MaterialKind.CATALYST)
    circle = graph.product(
        "cdna_circle", "covalently closed single-stranded cDNA circle",
        MaterialKind.MOLECULAR,
        chemical_state=CIRCLE,
        state=tuple([*CIRCLE.rows("CircLigase II"),
                     *annotation_rows(linear_cdna(), prefix_width=6)]))
    graph.react(
        "circularize",
        "Join the cDNA 3-prime hydroxyl to its 5-prime phosphate (MnCl2, betaine, 60 °C 1 h)",
        StepKind.CIRCULARIZATION,
        ((linear, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,),
        geometry=CircularizationGeometry.known(
            "single-strand intramolecular ligation",
            statement=("CircLigase II directly closes the phosphorylated single-stranded "
                       "cDNA by joining its 3-prime and 5-prime ends.")),
        relations=(MolecularRelation.CIRCLE_CLOSURE,))

    rca_primer = graph.input(
        "rca_primer", "phosphorothioate-protected RCA primer", MaterialKind.OLIGO,
        state=Construct([seg("RCA primer", RCA_PRIMER, "r2")], name="RCA primer"))
    primed_circle = graph.product(
        "primed_cdna_circle", "RCA-primer-bound cDNA circle", MaterialKind.MOLECULAR,
        chemical_state=CIRCLE,
        state=tuple(rca_primed_rows()))
    graph.react(
        "bind_rca_primer",
        "Anneal 500 nM phosphorothioate-protected RCA primer (60 °C 1 h, 60 °C wash)",
        StepKind.HYBRIDIZATION,
        ((circle, InputRole.TEMPLATE), (rca_primer, InputRole.STRUCTURAL)),
        (primed_circle,))
    phi29 = graph.input("phi29", "phi29 strand-displacing DNA polymerase",
                        MaterialKind.CATALYST)
    rolony_product = graph.product(
        "rolony_product", "aminoallyl-bearing rolling-circle concatemer",
        MaterialKind.MOLECULAR, state=rolony(), chemical_state=RCA)
    graph.react(
        "rolling_circle",
        "Amplify the cDNA circle into a tandem-repeat rolony (Phi29, 30 °C overnight)",
        StepKind.EXTENSION,
        ((primed_circle, InputRole.TEMPLATE), (phi29, InputRole.CATALYST),
         (aminoallyl_dutp, InputRole.STRUCTURAL)),
        (rolony_product,), relations=(MolecularRelation.ROLLING_CIRCLE_COPY,))

    library = graph.product("sequencing_library", "cross-linked FISSEQ rolony",
                            MaterialKind.LIBRARY, state=rolony())
    graph.react("immobilize_rolony",
                "Cross-link the amplified rolony in place (BS(PEG)9 1 h, Tris quench)",
                StepKind.CROSS_LINKING,
                ((rolony_product, InputRole.TEMPLATE),
                 (crosslinker, InputRole.CONDITION)), (library,))
    graph.sequence(library, SequencingPlan(
        acquired_roles=("Sequencing output",),
        limitation=("FISSEQ uses five offset SOLiD sequencing-by-ligation primers on "
                    "every repeat of an in-situ rolony, rather than one unique linear "
                    "R1/I1/I2/R2 primer landing site. Each primer gets seven T4 DNA ligase "
                    "cycles with SOLiD ToP eight-base probes and fluorophore cleavage. The "
                    "Science paper reports 27-base fibroblast reads; no acquired read "
                    "length is modelled.")))
    graph.declare_reviewed_scope((
        "cellular_rna", "rt_primer", "reverse_transcriptase", "aminoallyl_dutp",
        "reverse_transcribe", "crosslinker", "immobilize_cdna", "ribonucleases",
        "remove_rna", "circligase", "circularize", "rca_primer", "bind_rca_primer", "phi29",
        "rolling_circle", "immobilize_rolony", "sequencing_library"), SOURCE)
    graph.annotate(DetailView(
        "cdna_circle", "Circular template and RCA primer",
        "The complete cDNA is covalently closed. The protected RCA primer binds "
        "the adapter antiparallel; Phi29 extends from its 3′ end, and the resulting "
        "concatemer is shown as the next molecular state in the workflow.",
        state=RCA_TOPOLOGY, slot=AnnotationSlot.BEFORE_WORKFLOW,
        style=DetailStyle.FOLDED))
    graph.annotate(SourceDiscrepancy(
        "sequencing_library", "Ligation start positions listed in Box 1",
        (SourceClaim("Lee et al. 2015 Nature Protocols, Fig. 1b legend and procedure overview",
                     "each primer's seven ligation cycles interrogate every fifth "
                     "dinucleotide pair"),
         SourceClaim("Lee et al. 2015 Nature Protocols, Box 1",
                     "primer N positions 1-6-11-16-20-26-31 and primer N-4 positions "
                     "2-7-12-16-22-27-32; primer N-1 is listed as 0-5-10-15-20-25-30")),
        "The page shows the five primers at their binding site and does not draw a "
        "per-position read map.",
        "A five-base stride predicts 21 (not 20) for primer N and 17 (not 16) for "
        "primer N-4; as printed, position 20 is assigned to two primers and 21 and 17 "
        "to none.",
        "Box 1 does not explain the two values, so whether they are typographical "
        "errors is not established by the source.", AnnotationSlot.SEQUENCING))
    graph.annotate(DetailView(
        "sequencing_library", "Sequencing primers on the finished rolony",
        "All five sliding primers are shown on the middle repeat. The three primer-N "
        "sites in this three-repeat drawing are found computationally; every physical "
        "repeat carries the same nested sites.",
        state=sequencing_scene(), slot=AnnotationSlot.SEQUENCING))
    return graph.require_publishable()


def render_page() -> str:
    graph = protocol_graph()
    return assemble_page(graph, [head("FISSEQ library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>FISSEQ &mdash; transcriptome sequencing inside fixed cells</h1>',
        '<p class="research-notes"><a href="01_fisseq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1126/science.1250212">Lee et al., <i>Science</i> (2014)</a>. Exact oligos and reaction order are from the authors\' open <a href="https://doi.org/10.1038/nprot.2014.191"><i>Nature Protocols</i> procedure</a>.'),
        AnnotationSlot.BEFORE_WORKFLOW,
        WorkflowSection("Reaction workflow"),
        AnnotationSlot.SEQUENCING,
        OffsetPrimerSection(
            "sequencing_library", SEQUENCING_PRIMERS,
            "SOLiD sequencing-by-ligation uses primer N, then four primers shortened "
            "successively at the 5′ end. The offsets make the five-base ligation stride "
            "interrogate all positions.",
            "in-situ rolony cDNA interrogated by five offset SOLiD ligation primers",
            "FISSEQ uses SOLiD sequencing-by-ligation rather than an R1/I1/I2/R2 run."),
        PlanLimitationSection("sequencing_library"),
        '</div>'])


def selftest() -> None:
    """Check oligos and incubation against the Lee et al. procedure."""
    from checks import Check, run_common
    check = Check()
    check.section("Lee et al. protocol choices")
    check("RT primer as published", RT_PRIMER_WRITTEN,
          "/5phos/TCTCGGGAACGCTGAAGANNNNNN")
    check("RCA primer as published", RCA_PRIMER_WRITTEN, "TCTTCAGCGTTCCCGA*G*A")
    check("CircLigase II incubation", (CIRCLIGASE_C, CIRCLIGASE_MIN), (60, 60))
    run_common(check)
    check.report()
