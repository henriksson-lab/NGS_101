"""Construct model for the custom-Nextera itChIP-seq library."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# Verbatim representative oligos from Ai et al., Supplementary Table 1.
SOURCE_T5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
SOURCE_T7_1 = "GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
SOURCE_I5_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC"
SOURCE_I7_N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG"

CONNECTOR_A = "TCGTCGGCAGCGTCTCCACGC"
CONNECTOR_B = "GTCTCGTGGGCTCGGCTGTCCCTGTCC"
T5_BARCODE = "TATAGCCT"
T7_BARCODE = "CGAGTAAT"
ME_A_SPACER = "GCGATCGAGGACGGC"
ME_B_SPACER = "CACCGTCTCCGCCTC"
EXTENDED_ME_A = ME_A_SPACER + nx.ME
EXTENDED_ME_B = ME_B_SPACER + nx.ME
COMMON_ANNEALING = nx.ME_RC
I5_OLIGO_INDEX = "TAGATCGC"
I7_OLIGO_INDEX = "TCGCCTTA"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def t5_segments() -> list[Segment]:
    return [seg("connector A", CONNECTOR_A, "s5"),
            seg("T5 barcode", T5_BARCODE, "cbc",
                feature=feature("cell_t5", "cell_barcode", "combinatorial", group="cell_id", part="T5")),
            seg("ME-A spacer", ME_A_SPACER, "r1"),
            seg("ME", nx.ME, "me")]


def t7_segments() -> list[Segment]:
    return [seg("connector B", CONNECTOR_B, "s7"),
            seg("T7 barcode", T7_BARCODE, "cbc",
                feature=feature("cell_t7", "cell_barcode", "combinatorial", group="cell_id", part="T7")),
            seg("ME-B spacer", ME_B_SPACER, "r2"),
            seg("ME", nx.ME, "me")]


def loaded_end_scene(which: str) -> Scene:
    """Representative transferred strand paired to the phosphorylated ME strand."""
    if which not in ("T5", "T7"):
        raise ValueError("which must be 'T5' or 'T7'")
    transferred = t5_segments() if which == "T5" else t7_segments()
    common = [seg("ME'", COMMON_ANNEALING, "me")]
    sc = Scene()
    sc.strand("transferred", transferred, label=f"{which}-1 oligo")
    sc.anneal("non-transferred", common, to="transferred",
              pair=("ME'", "ME"), label="common annealing primer", mod5="p")
    sc.mark("transferred", f"{which} barcode", "cell barcode")
    return sc


def indexed_p5_segments() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_OLIGO_INDEX, "cbc"),
            seg("s5", nx.S5, "s5")]


def indexed_p7_segments() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7 oligo", I7_OLIGO_INDEX, "cbc"),
            seg("s7", nx.S7, "s7")]


def gap_filled_fragment() -> Construct:
    """A heterotypic T5--T7 product after the tagmentation gaps are filled."""
    right = complement_segments(t7_segments(), suffix="'")
    return Construct([*t5_segments(),
                      seg("immunoprecipitated DNA", "XXXXXXXX...XXXXXXXX",
                          placeholder=True),
                      *right], name="gap-filled T5--T7 fragment")


def final_library() -> Construct:
    """PCR-completed custom-Nextera library, top strand 5' to 3'."""
    right = complement_segments(t7_segments(), suffix="'")
    return Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", I5_OLIGO_INDEX, "cbc",
            feature=feature("sample_i5", "sample_index", "whitelist", whitelist="published Nextera PCR index set")),
        *t5_segments(),
        seg("immunoprecipitated DNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *right,
        seg("i7", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published Nextera PCR index set")),
        seg("P7'", il.P7_RC, "p7"),
    ], name="itChIP-seq custom-Nextera library")


SEQ_PRIMERS = [
    sp.custom("Read 1", "itChIP Read 1", EXTENDED_ME_A,
              "Ai et al. Supplementary Table 1"),
    sp.custom("Index 1 (i7)", "itChIP Index 1", revcomp(EXTENDED_ME_B),
              "Ai et al. Supplementary Table 1"),
    sp.custom("Index 2 (i5)", "grafted P5", il.P5,
              "Ai et al. Supplementary Table 1"),
    sp.custom("Read 2", "itChIP Read 2", EXTENDED_ME_B,
              "Ai et al. Supplementary Table 1"),
]
READ_LENGTHS = {}  # acquired cycle mask remains unverified in the capped source scope


def _validate_model() -> None:
    if CONNECTOR_A[:len(nx.S5)] != nx.S5 or CONNECTOR_B[:len(nx.S7)] != nx.S7:
        raise ValueError("itChIP connectors no longer begin with canonical s5/s7")
    if EXTENDED_ME_A[-len(nx.ME):] != nx.ME or EXTENDED_ME_B[-len(nx.ME):] != nx.ME:
        raise ValueError("extended mosaic ends no longer terminate in canonical ME")
    if "".join(s.top for s in t5_segments()) != SOURCE_T5_1:
        raise ValueError("assembled T5-1 differs from Supplementary Table 1")
    if "".join(s.top for s in t7_segments()) != SOURCE_T7_1:
        raise ValueError("assembled T7-1 differs from Supplementary Table 1")
    if "".join(s.top for s in indexed_p5_segments()) != SOURCE_I5_N501:
        raise ValueError("assembled i5-N501 differs from Supplementary Table 1")
    if "".join(s.top for s in indexed_p7_segments()) != SOURCE_I7_N701:
        raise ValueError("assembled i7-N701 differs from Supplementary Table 1")
    loaded_end_scene("T5").rows()
    loaded_end_scene("T7").rows()
    Scene.duplex(list(gap_filled_fragment())).rows()
    Scene.duplex(list(final_library())).rows()
    problems = sp.verify(final_library(), SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid itChIP sequencing layout: " + "; ".join(problems))


_validate_model()


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41556-019-0383-5"
    graph = ProtocolGraph("itChIP-seq", authority=source)
    chromatin = graph.input("chromatin", "opened chromatin in first-round wells",
                            MaterialKind.BIOLOGICAL, primary=True)
    transposomes = graph.input("barcoded_transposomes", "paired T5/T7 barcoded Tn5 ends",
                               MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_chromatin", "T5/T7-tagged chromatin",
                           MaterialKind.MOLECULAR)
    tagmentation_waste = graph.product("tagmentation_waste",
                                       "candidate non-library tagmentation products", MaterialKind.MOLECULAR)
    graph.react("tagmentation", "Tagment opened chromatin with one T5/T7 pair per well",
                StepKind.TAGMENTATION,
                ((chromatin, InputRole.TEMPLATE), (transposomes, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tagmentation_waste))
    graph.mark_terminal(Terminal(tagmentation_waste.id, Disposition.UNRESOLVED, "Amplification eligibility does not demonstrate physical discard of candidate insertion products."))
    antibody = graph.input("antibody", "ChIP antibody", MaterialKind.SELECTION)
    selected = graph.product("selected_chromatin", "immunoprecipitated tagged chromatin",
                             MaterialKind.MOLECULAR)
    unselected = graph.product("unselected_chromatin", "unselected chromatin",
                               MaterialKind.WASTE)
    graph.react("chip", "Pool, release chromatin and immunoprecipitate the target",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (selected, unselected))
    graph.discard(unselected, "chromatin not retained by immunoprecipitation is washed away")
    polymerase = graph.input("polymerase", "gap-fill and PCR polymerase",
                             MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterotypic T5/T7 fragments",
                           MaterialKind.MOLECULAR, state=gap_filled_fragment())
    graph.react("gap_fill", "Fill the transposase-generated gaps", StepKind.EXTENSION,
                ((selected, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (filled,))
    i5 = graph.input("i5_primer", "indexed i5 PCR primer", MaterialKind.OLIGO)
    i7 = graph.input("i7_primer", "indexed i7 PCR primer", MaterialKind.OLIGO)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify with one i5 and one i7 primer", StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Supplementary Table 1 gives exact oligos; their four custom sequencing roles make a single canonical PCR landing view misleading."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS), read_lengths=READ_LENGTHS,
        limitation="The bounded original supplement confirms custom indexing architecture; actual cycle mask remains unassigned pending main Methods.",
        read_uses={"Read 1": "include", "Index 1 (i7)": "extract",
                   "Index 2 (i5)": "extract", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "chromatin", "barcoded_transposomes", "tn5", "tagmentation", "antibody",
        "chip", "polymerase", "gap_fill", "i5_primer", "i7_primer",
        "library_pcr", "library"), source)
    return graph.require_publishable()
