"""Local source and reaction model for miracl-seq."""
from chemdraw import Construct, Scene, Segment, Workflow, feature
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import role
import illumina as il
import seqprimers as sp
from rare_cell_and_affinity import _finish_pcr, _scene, _workflow
from types import SimpleNamespace

TRANSCRIPT = role("transcript-derived insert")
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x GemCode v2")
UMI = feature("umi", "umi", "random")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "unknown")

from dataclasses import dataclass

@dataclass(frozen=True)
class MiraclLoading:
    coefficient: int

    def __post_init__(self) -> None:
        if self.coefficient < 2:
            raise ValueError("MIRACL-seq is defined by deliberate droplet overloading")

    def partition(self) -> Scene:
        return _scene(f"{self.coefficient}× overloaded GEM",
                      role("one 10x gel bead"), role("multiple RAISIN nuclei"),
                      role("oil-emulsion partition"))

def raisin_particle() -> Scene:
    sc = Scene()
    sc.strand("nucleus", [role("intact nucleus")], label="RAISIN")
    sc.strand("ribosome", [role("rough ER / attached ribosomes"),
                            role("ribosome-bound poly(A) RNA")], label="retained exterior RNA")
    sc.note("ribosome", "CHAPS/salt extraction retains nucleus-associated ribosomes and mRNA")
    return sc

def gemcode_v2_library() -> tuple[Construct, tuple]:
    """The standard 10x v2 endpoint used by droplet RAISIN and MIRACL-seq."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                feature=CELL_BARCODE),
        Segment("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
        Segment("poly(dT)30", "T" * 30),
        Segment("V", "V", placeholder=True), Segment("N", "N", placeholder=True),
        Segment("RAISIN transcript-derived insert", "X" * 36, placeholder=True),
        Segment("dA junction", "A"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="GemCode v2 RAISIN library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid GemCode v2 endpoint: " + "; ".join(problems))
    return lib, primers

def miracl_workflow(coefficient: int = 2) -> Workflow:
    loading = MiraclLoading(coefficient)
    final, _ = gemcode_v2_library()
    barcoded = _scene("barcoded first strand", role("10x Read 1 handle"),
                      role("16-nt cell barcode", molecular_feature=CELL_BARCODE),
                      role("10-nt UMI", molecular_feature=UMI), TRANSCRIPT)
    return _workflow(raisin_particle(), (
        (f"Overload the GemCode v2 channel ({coefficient}× shown)", loading.partition(),
         "The defining experiments tested 2×–8× loading; rare-cell homotypic multiplets remain unlikely."),
        ("Lyse in GEMs and perform barcoded reverse transcription", barcoded,
         "One bead barcode labels cDNA from every nucleus sharing that partition."),
        ("Amplify, shear, attach the 5′ adapter and sample index",
         Scene.duplex(list(final), label="MIRACL-seq library"),
         "This is the standard GemCode v2 molecular endpoint; heterotypic multiplets are removed computationally."),
    ))

final, primers = gemcode_v2_library()
SPEC = SimpleNamespace(
    TITLE="MIRACL-seq — droplet overloading for rare-cell RNA-seq",
    NOTES="01_miracl_seq.html",
    SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.08.003">Drokhlyansky et al., Cell (2020)</a>.',
    SUMMARY="MIRACL-seq deliberately overloads GemCode v2 droplets with unsorted RAISIN nuclei, increasing rare-cell recovery while computationally removing heterotypic multiplets.",
    CAVEAT="MIRACL-seq changes partition loading, not the GemCode v2 oligo chemistry; the final construct therefore uses the standard v2 molecular endpoint.",
    FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers,
    READ_LENGTHS={"Read 1": 26},
    FINAL_CAPTION="GemCode v2 3′ gene-expression library from overloaded RAISIN partitions.",
    SEQUENCING_INTRO="Read 1 records the 16-nt bead barcode and 10-nt UMI; Read 2 interrogates the transcript-derived insert.",
    SEQUENCING_ENDING="",
    SEQUENCING_UNAVAILABLE="",
    workflow=miracl_workflow,
    oligos=lambda: (),
)

def _miracl_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    particles = graph.input("raisin_particles", "unsorted RAISIN nuclei",
                            MaterialKind.BIOLOGICAL, state=raisin_particle(), primary=True)
    droplets = graph.input("overload_condition", "deliberately overloaded GemCode-v2 partition condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_gems", "overloaded GEMs containing one bead and multiple RAISIN nuclei",
                             MaterialKind.MOLECULAR, state=MiraclLoading(2).partition())
    unresolved = graph.product("unresolved_gems", "empty or uninterpretable GEMs",
                               MaterialKind.WASTE)
    graph.react("partition", "Overload RAISIN nuclei into GemCode-v2 partitions",
                StepKind.PARTITION,
                ((particles, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, unresolved), changes_molecule=False)
    graph.discard(unresolved, "does not yield a usable gel-bead barcode partition")
    bead_primer = graph.input("bead_primer", "GemCode-v2 barcode/UMI oligo(dT) bead primer",
                              MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "droplet reverse transcriptase",
                        MaterialKind.CATALYST)
    cdna = graph.product("barcoded_cdna", "partition-barcoded RAISIN transcript cDNA",
                         MaterialKind.MOLECULAR,
                         state=_scene("barcoded first strand", role("16-nt cell barcode",
                             molecular_feature=CELL_BARCODE), role("10-nt UMI",
                             molecular_feature=UMI), TRANSCRIPT))
    graph.react("reverse_transcribe", "Lyse in GEMs and perform barcoded reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cdna,))
    reviewed = ["raisin_particles", "overload_condition", "partition", "bead_primer",
                "reverse_transcriptase", "reverse_transcribe"]
    _finish_pcr(graph, "library", cdna, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="MIRACL-seq GemCode-v2")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _miracl_graph(SPEC)

SPEC.protocol_graph = protocol_graph
