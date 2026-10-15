"""Molecular model for VASA-seq as defined by Salmen et al. (2022).

Supplementary Table 12 (`MOESM13`) prints the shared RA3/RTP pair and all 46+46
VASA-drop index primers; Supplementary Table 14 (`MOESM15`) names, for every oligo
it does *not* print, the publication it came from: the VASA-plate CEL-seq2 capture
primers to Muraro et al. 2016 (doi:10.1016/j.cels.2016.09.002, barcodes in
GSE176588), the VASA-plate RA3/RTP to Hashimshony et al. 2012
(doi:10.1016/j.celrep.2012.08.003), and the inDrop v3 bead barcoding to Briggs et
al. 2018 (doi:10.1126/science.aar5780).  None of those documents is on disk, so
every region they would supply is an inferred placeholder here.  The scg_lib_structs
drawn page prints reconstructions of all of them; it is not used.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, revcomp
from illumina import P5, P7, TRUSEQ_SMALL_RNA_RA3
import nextera as nx
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


# Supplementary Table 12 (MOESM13) prints "/5rApp/TGGAATTCTCGGGTGCCAAGG/3SpC3/" and
# "GCCTTGGCACCCGAGAATTCCA"; those bases are the canonical TruSeq Small RNA RA3 and its
# matching RT primer, so RA3 is imported rather than retyped.
RA3 = TRUSEQ_SMALL_RNA_RA3
RTP = "GCCTTGGCACCCGAGAATTCCA"
# The 15-nt tail shared by every VD_ILMN8_i7 primer, verbatim from Supplementary Table 12.
DROP_COMMON = "GGGTGTCGGGTGCAG"
INDEX_LEN = 8
UFI_LEN = 6
DROP_BARCODE_LEN = 8
PLATE_BARCODE_LEN = 8
UFI_FEATURE = feature("ufi", "umi", "random")
PLATE_CELL = feature("cell_barcode_plate", "cell_barcode", "whitelist",
                     whitelist="vasa-plate-barcode-set")
DROP_BC1 = feature("cell_barcode_drop_1", "cell_barcode", "combinatorial",
                   group="cell_barcode", part="barcode 1")
DROP_BC2 = feature("cell_barcode_drop_2", "cell_barcode", "combinatorial",
                   group="cell_barcode", part="barcode 2")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UFI": UFI_FEATURE, "UFI'": UFI_FEATURE,
            "cell barcode": PLATE_CELL, "barcode 1'": DROP_BC1,
            "barcode 2'": DROP_BC2, "i5": I5_FEATURE, "i7": I7_FEATURE,
            "i7'": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, token: str, tag: str | None = None) -> Segment:
    """Create an unpublished molecular region with mandatory inferred styling."""
    if not (token.startswith("[") and token.endswith("]")):
        raise ValueError("inferred molecular regions must be bracketed role tokens")
    return seg(name, token, tag, placeholder=True, inferred=True)


def poly_a(name: str = "added poly(A)", display_bases: int = 12) -> Segment:
    """Display the published poly(A)-tailing reaction without asserting tail length."""
    if display_bases < 2:
        raise ValueError("poly(A) display must be legible")
    return seg(name, "A" * display_bases, inferred=True)


def repaired_fragment() -> Construct:
    return Construct([
        seg("RNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        poly_a(),
    ], name="fragmented, repaired and poly(A)-tailed cellular RNA")


def capture_primer(fmt: str) -> Construct:
    if fmt not in {"plate", "drop"}:
        raise ValueError("VASA format must be 'plate' or 'drop'")
    return Construct([
        inferred(f"{fmt} barcode / T7 region", f"[{fmt} capture region]", "cbc"),
        seg("UFI", "N" * UFI_LEN, "umi", placeholder=True),
        seg("cell barcode", "N" * (PLATE_BARCODE_LEN if fmt == "plate" else DROP_BARCODE_LEN),
            "cbc", placeholder=True),
        seg("oligo-dT", "T" * 12, inferred=True),
    ], name=f"VASA-{fmt} capture oligo (schematic)")


def ivt_template(fmt: str) -> Construct:
    p = capture_primer(fmt)
    return Construct([
        *p.segments,
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"VASA-{fmt} IVT template (schematic)")


def schematic_duplex(con: Construct) -> Scene:
    """Pair every sequence-like column; exempt only bracketed role placeholders."""
    unpaired = [s.name + "'" for s in con if s.placeholder and s.top.startswith("[")]
    return Scene.duplex(list(con), unpaired=unpaired)


def ligated_arna() -> Construct:
    return Construct([
        seg("rRNA-depleted aRNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("RA3", RA3, "r2"),
    ], name="rRNA-depleted aRNA after RA3 ligation")


def ra3_ligation_geometry():
    from adapter_ligation import AdapterLigationGeometry
    sc = Scene(); sc.strand("RA3", [seg("5′-rApp ligating end", RA3, "r2"),
                                    inferred("3′ blocking group", "[3SpC3]")],
                            label="published RA3")
    sc.ligation_end("RA3", "5′-rApp ligating end", kind="single-strand",
                    polarity="5-prime", owner="adapter",
                    text="5′-rApp single-strand end")
    return AdapterLigationGeometry.known("single-strand", sc,
        statement="Table 12 prints /5rApp/-RA3-/3SpC3/: its activated 5′ end ligates to aRNA 3′-OH, while the 3′ block prevents concatemerization; this is not a duplex overhang.")


def rtp_segments() -> list[Segment]:
    """RTP assembled so its published 21-nt annealing core is guaranteed from RA3."""
    core = revcomp(RA3)
    if RTP != "G" + core:
        raise ValueError("published RTP no longer equals 5' G plus reverse-complement RA3")
    return [seg("RTP 5' G", RTP[0]), seg("RTP core", core, "r2")]


def drop_i5_primer(index: str) -> Construct:
    if len(index) != INDEX_LEN:
        raise ValueError("VASA-drop i5 index must be 8 nt")
    return Construct([seg("P5", P5, "p5"), seg("i5", index, "cbc"),
                      seg("s5", nx.S5, "s5")], name="VD_ILMN8_i5")


def drop_i7_primer(index: str) -> Construct:
    if len(index) != INDEX_LEN:
        raise ValueError("VASA-drop i7 index must be 8 nt")
    return Construct([seg("P7", P7, "p7"), seg("i7", index, "cbc"),
                      seg("15-nt common sequence", DROP_COMMON)], name="VD_ILMN8_i7")


def plate_library() -> Construct:
    """Final VASA-plate molecule at the resolution supported by the defining paper."""
    con = Construct([
        inferred("plate P5-side PCR arm", "[plate P5-side arm]", "p5"),
        seg("UFI", "N" * UFI_LEN, "umi", placeholder=True),
        seg("cell barcode", "N" * PLATE_BARCODE_LEN, "cbc", placeholder=True),
        seg("poly(T) read-through", "T" * 12),
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("RA3", RA3, "r2"),
        inferred("plate P7-side PCR arm", "[plate P7-side arm]", "p7"),
    ], name="VASA-plate final library")
    _require_order(con, ("UFI", "cell barcode", "poly(T) read-through", "fragment cDNA", "RA3"))
    return con


def drop_library() -> Construct:
    """Final VASA-drop molecule using paper-supported read boundaries.

    The disputed ligation-adapter/bead-linker sequence is deliberately represented
    only by a role placeholder.
    """
    con = Construct([
        *drop_i5_primer("N" * INDEX_LEN).segments,
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        poly_a("poly(A) copy"),
        seg("UFI'", "N" * UFI_LEN, "umi", placeholder=True),
        seg("barcode 2'", "N" * DROP_BARCODE_LEN, "cbc", placeholder=True),
        inferred("bead linker", "[proprietary bead linker]"),
        seg("barcode 1'", "N" * DROP_BARCODE_LEN, "cbc", placeholder=True),
        seg("common sequence'", revcomp(DROP_COMMON)),
        seg("i7'", "N" * INDEX_LEN, "cbc", placeholder=True),
        seg("P7'", revcomp(P7), "p7"),
    ], name="VASA-drop final library")
    _require_order(con, ("fragment cDNA", "UFI'", "barcode 2'", "bead linker",
                         "barcode 1'", "common sequence'", "i7'"))
    return con


def _require_order(con: Construct, ordered: tuple[str, ...]) -> None:
    names = [s.name for s in con]
    try:
        positions = [names.index(name) for name in ordered]
    except ValueError as exc:
        raise ValueError(f"{con.name} lacks a read-layout segment") from exc
    if positions != sorted(positions):
        raise ValueError(f"{con.name} read-layout segments are out of order")


def read_layouts() -> tuple[tuple[str, str, str], ...]:
    """Read boundaries reported by the defining paper, checked against the models."""
    plate, drop = plate_library(), drop_library()
    _require_order(plate, ("UFI", "cell barcode", "poly(T) read-through"))
    _require_order(drop, ("fragment cDNA", "UFI'", "barcode 2'", "barcode 1'",
                          "common sequence'", "i7'"))
    return (
        ("VASA-plate", "Read 1 · 26 cycles", "UFI 6 + cell barcode 8 + poly(T)"),
        ("VASA-plate", "Read 2 · 135 cycles", "fragment cDNA"),
        ("VASA-drop", "Read 1 · 247 cycles", "fragment cDNA"),
        ("VASA-drop", "Index 1 · 31 cycles", "barcode 1 (8) + common sequence (15) + i7 (8)"),
        ("VASA-drop", "Index 2 · 8 cycles", "i5"),
        ("VASA-drop", "Read 2 · 14 cycles", "barcode 2 (8) + UFI (6)"),
    )


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41587-022-01361-8"
    unresolved = source + "; capture oligos and complete terminal PCR arms are not printed"
    graph = ProtocolGraph("VASA-seq", authority=source)
    rna = graph.input("total_rna", "cellular total RNA", MaterialKind.BIOLOGICAL,
                      primary=True)
    heat = graph.input("heat_condition", "heat-fragmentation condition",
                       MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_rna", "heat-fragmented total RNA",
                               MaterialKind.MOLECULAR)
    graph.react("fragmentation", "Heat-fragment total RNA", StepKind.FRAGMENTATION,
                ((rna, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (fragmented,))
    repair = graph.input("repair_enzymes", "RNA end-repair enzymes", MaterialKind.CATALYST)
    repaired = graph.product("repaired_rna", "end-repaired RNA fragments",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "Repair fragmented RNA ends", StepKind.END_REPAIR,
                ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("poly_a_polymerase", "poly(A) polymerase", MaterialKind.CATALYST)
    tailed = graph.product("polya_rna", repaired_fragment().name, MaterialKind.MOLECULAR,
                           state=repaired_fragment())
    graph.react("poly_a_tailing", "Add poly(A) to every repaired RNA fragment",
                StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    plate_input = graph.product("plate_rna", "VASA-plate RNA aliquot", MaterialKind.MOLECULAR,
                                state=repaired_fragment())
    drop_input = graph.product("drop_rna", "VASA-drop RNA aliquot", MaterialKind.MOLECULAR,
                               state=repaired_fragment())
    graph.react("format_split", "Split into VASA-plate and VASA-drop formats", StepKind.SPLIT,
                ((tailed, InputRole.TEMPLATE),), (plate_input, drop_input),
                changes_molecule=False)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    ds_polymerase = graph.input("ds_polymerase", "second-strand DNA polymerase",
                                MaterialKind.CATALYST)
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    rnase_h = graph.input("rnase_h", "thermostable RNase H rRNA-depletion reagent",
                          MaterialKind.SELECTION)
    ra3 = graph.input("ra3", "published 5′-rApp RA3 adapter", MaterialKind.ADAPTER,
                      state=(seg("RA3", RA3, "r2"),))
    ligase = graph.input("rna_ligase", "RNA adapter ligase", MaterialKind.CATALYST)
    rtp = graph.input("rtp", "published RTP reverse-transcription primer",
                      MaterialKind.OLIGO, state=tuple(rtp_segments()))
    pcr_polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    branch_data = (
        ("plate", plate_input, plate_library(),
         "VASA-plate capture oligo and PCR arms are cited but not printed.",
         {"Read 1": 26, "Read 2": 135},
         {"Read 1": "extract", "Read 2": "include"}),
        ("drop", drop_input, drop_library(),
         "The VASA-drop bead linker and complete sequencing-primer landing sites are not published.",
         {"Read 1": 247, "Index 1 (i7)": 31, "Index 2 (i5)": 8, "Read 2": 14},
         {"Read 1": "include", "Index 1 (i7)": "extract",
          "Index 2 (i5)": "demultiplex", "Read 2": "extract"}),
    )
    reviewed = [
        "total_rna", "heat_condition", "fragmentation", "repair_enzymes", "end_repair",
        "poly_a_polymerase", "poly_a_tailing", "format_split", "reverse_transcriptase",
        "ds_polymerase", "t7_polymerase", "rnase_h", "ra3", "rna_ligase", "rtp",
        "pcr_polymerase"]
    for fmt, branch_input, final, limitation, lengths, uses in branch_data:
        capture = graph.input(f"{fmt}_capture_oligo", f"VASA-{fmt} capture oligo",
                              MaterialKind.OLIGO, state=capture_primer(fmt), inferred=True,
                              source=unresolved)
        first = graph.product(f"{fmt}_first_strand", f"VASA-{fmt} barcoded first-strand cDNA",
                              MaterialKind.MOLECULAR, inferred=True, source=unresolved)
        graph.react(f"{fmt}_rt", f"Barcode and reverse-transcribe VASA-{fmt} RNA",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((branch_input, InputRole.TEMPLATE), (capture, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (first,), source=unresolved)
        ivt_template_state = ivt_template(fmt)
        template = graph.product(f"{fmt}_ivt_template", ivt_template_state.name,
                                 MaterialKind.MOLECULAR, state=ivt_template_state,
                                 inferred=True, source=unresolved)
        graph.react(f"{fmt}_second_strand", f"Make the VASA-{fmt} double-stranded IVT template",
                    StepKind.EXTENSION,
                    ((first, InputRole.TEMPLATE), (ds_polymerase, InputRole.CATALYST)),
                    (template,), source=unresolved)
        arna = graph.product(f"{fmt}_arna", f"VASA-{fmt} amplified RNA",
                             MaterialKind.MOLECULAR)
        graph.react(f"{fmt}_ivt", f"T7-amplify the VASA-{fmt} template",
                    StepKind.TRANSCRIPTION,
                    ((template, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,))
        depleted = graph.product(f"{fmt}_depleted", f"rRNA-depleted VASA-{fmt} aRNA",
                                 MaterialKind.MOLECULAR)
        rrna = graph.product(f"{fmt}_rrna_waste", f"removed VASA-{fmt} rRNA-derived aRNA",
                             MaterialKind.WASTE)
        graph.react(f"{fmt}_depletion", f"Deplete rRNA from VASA-{fmt} aRNA",
                    StepKind.SELECTION,
                    ((arna, InputRole.TEMPLATE), (rnase_h, InputRole.SELECTION)),
                    (depleted, rrna))
        graph.discard(rrna, "rRNA-derived amplified RNA is removed before library construction")
        ligated = graph.product(f"{fmt}_ra3_ligated", f"RA3-ligated VASA-{fmt} aRNA",
                                MaterialKind.MOLECULAR, state=ligated_arna())
        graph.react(f"{fmt}_ra3_ligation", f"Ligate RA3 to VASA-{fmt} aRNA",
                    StepKind.LIGATION,
                    ((depleted, InputRole.TEMPLATE), (ra3, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,), geometry=ra3_ligation_geometry())
        second = graph.product(f"{fmt}_second_cdna", f"RTP-primed VASA-{fmt} cDNA",
                               MaterialKind.MOLECULAR)
        graph.react(f"{fmt}_second_rt", f"Reverse-transcribe VASA-{fmt} aRNA from RTP",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((ligated, InputRole.TEMPLATE), (rtp, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (second,))
        primers = graph.input(f"{fmt}_pcr_primers", f"VASA-{fmt} terminal PCR-primer set",
                              MaterialKind.OLIGO, inferred=True, source=unresolved)
        library = graph.product(f"{fmt}_library", final.name, MaterialKind.LIBRARY,
                                state=final, inferred=True, source=unresolved)
        graph.react(f"{fmt}_pcr", f"PCR-complete the VASA-{fmt} library", StepKind.PCR,
                    ((second, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
                     (pcr_polymerase, InputRole.CATALYST)), (library,),
                    geometry=PCRGeometry.unavailable(
                        final, reason=limitation, symmetric_species=True), source=unresolved)
        graph.sequence(library, SequencingPlan(
            limitation=limitation, acquired_roles=tuple(lengths),
            read_lengths=lengths, read_uses=uses))
        reviewed.extend((
            f"{fmt}_capture_oligo", f"{fmt}_rt", f"{fmt}_second_strand",
            f"{fmt}_ivt", f"{fmt}_depletion", f"{fmt}_ra3_ligation",
            f"{fmt}_second_rt", f"{fmt}_pcr_primers", f"{fmt}_pcr",
            f"{fmt}_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
