"""Reusable molecular states for chromatin and epigenetic library protocols."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature, junction_row
from restriction import Digest, MBOI


def illumina_ligation_library(label: str, *, inferred: bool = True):
    """A conventional end-repair/dA/adapter-ligation library.

    Historical papers often name the vendor kit without printing its oligos.  In that
    case inference styling is attached here and cannot be forgotten by page callers.
    """
    return truseq_library(
        [seg(label, "X" * 36, placeholder=True)], label + " library",
        inferred_adapters=inferred)


def ligated_insert_scene(label: str, *, inferred_adapters: bool = True) -> Scene:
    con, _ = illumina_ligation_library(label, inferred=inferred_adapters)
    sc = Scene.duplex(list(con), label="adapter-ligated molecule")
    sc.junction("top", "Read 1 arm", label, "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc


def tagmented_library(label: str):
    return nextera_library([seg(label, "X" * 36, placeholder=True)], label + " library")


def tagmentation_scene(label: str) -> Scene:
    con, _ = tagmented_library(label)
    sc = Scene.duplex(list(con), label="PCR-completed tagmented molecule")
    sc.junction("top", "left mosaic end", label, "Tn5 transfer")
    sc.junction("top", label, "right mosaic end reverse complement", "Tn5 transfer")
    sc.labels("top")
    return sc


def droplet_targeted_tagment_library(label: str, *, assay_barcode_len: int = 8,
                                     assay_name: str = "assay barcode",
                                     connector_seq: str = "",
                                     assay_is_feature: bool = True):
    """Representative 10x-ATAC-backed targeted-tagmentation molecule.

    Several multimodal chromatin methods put a modality barcode between the S5 handle
    and mosaic end before passing nuclei into Chromium ATAC.  This constructor keeps
    cell identity and modality identity distinct and deliberately declares only primer
    sites that remain contiguous after that insertion.
    """
    if assay_barcode_len < 1:
        raise ValueError("targeted tagmentation needs a non-empty assay barcode")
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("10x cell barcode", "B" * 16, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "whitelist",
                            whitelist="10x Chromium ATAC whitelist")),
        seg("S5", nx.S5, "s5"),
        seg(assay_name, "M" * assay_barcode_len, "cbc" if assay_is_feature else None,
            placeholder=True,
            feature=(feature("assay_barcode", "feature_barcode", "whitelist",
                             whitelist="protocol oligo table")
                     if assay_is_feature else None)),
        *([seg("assay-side connector", connector_seq)] if connector_seq else []),
        seg("left mosaic end", nx.ME, "me"),
        seg(label, "X" * 38, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name=label + " targeted-tagmentation library")
    primers = (
        sp.custom("Read 1", "S5-side genomic read", nx.S5,
                  "protocol-specific targeted-tagmentation construct"),
        sp.NEXTERA["I1"],
        sp.custom("Read 2", "S7-side genomic read", nx.S7,
                  "protocol-specific targeted-tagmentation construct"),
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid targeted-tagmentation library: " + "; ".join(problems))
    return lib, primers


def spatial_barcoded_library(label: str, *, barcode_len: int = 8):
    """TruSeq-shell representation of orthogonal in-tissue ligation barcodes."""
    insert = [
        seg("spatial barcode A", "A" * barcode_len, "cbc", placeholder=True,
            feature=feature("spatial_a", "spatial_barcode", "combinatorial",
                            group="spatial_pixel", part="row", whitelist="published A plate")),
        seg("spatial barcode B", "B" * barcode_len, "cbc", placeholder=True,
            feature=feature("spatial_b", "spatial_barcode", "combinatorial",
                            group="spatial_pixel", part="column", whitelist="published B plate")),
        seg(label, "X" * 38, placeholder=True),
    ]
    return truseq_library(insert, label + " spatial library", inferred_adapters=True)


def pbat_priming_scene(round_name: str, tail_name: str, tail_seq: str,
                       *, mod5: str = "") -> Scene:
    template = [seg("bisulfite-converted template", "X" * 26, placeholder=True)]
    primer = [seg(tail_name, tail_seq, "r1"),
              seg("random N4", "N" * 4, placeholder=True)]
    sc = Scene(); sc.strand("template", template, label=round_name + " template")
    sc.anneal("primer", primer, to="template",
              pair=("random N4", "bisulfite-converted template"),
              shift=len(template[0]) - len(primer[-1]),
              unpaired=("random N4",), label=round_name + " random primer", mod5=mod5)
    sc.mark("primer", tail_name, "5′ sequencing-adapter tail")
    sc.arrow("primer", "polymerase extension")
    return sc


def proximity_junction() -> Construct:
    """MboI GATCGATC contact scar with derived fill and biotin coordinates."""
    filled = Digest(MBOI, "X" * 12, "X" * 12).fill_in(biotin_base="A")
    scar = filled.junction
    top_bio = filled.junction_biotin_top[0]
    boundary = MBOI.bottom_cut
    return Construct([
        seg("locus A", "X" * 20, placeholder=True),
        seg("filled end A before biotin", scar[:top_bio], "me"),
        seg("top-strand biotin-dA", scar[top_bio:top_bio + 1], "w1"),
        seg("filled end A after biotin", scar[top_bio + 1:boundary], "me"),
        seg("filled restriction end B", scar[boundary:], "me"),
        seg("locus B", "X" * 20, placeholder=True),
    ], name="proximity-ligation product")


def proximity_rows() -> list[Row]:
    con = proximity_junction()
    return [*Scene.duplex(list(con), label="contact product").rows(),
            junction_row(con, "filled end A after biotin", "filled restriction end B",
                         "proximity ligation")]


def antibody_enzyme_rows(enzyme: str, action: str) -> list[Row]:
    return [
        Row(chunks=[("chromatin ---- [target protein or histone mark] ---- chromatin", None, False)]),
        Row(chunks=[("                         ^ primary antibody", "cbc", False)]),
        Row(chunks=[(f"                         | Protein A–{enzyme}", "me", False)]),
        Row(chunks=[(f"                         └─ {action}", None, False)]),
    ]


def dimelo_mark_rows() -> list[Row]:
    return antibody_enzyme_rows("Hia5", "SAM-dependent m6A marks deposited on nearby DNA") + [
        Row(chunks=[("single long molecule: ---A--m6A---A----m6A------A---", "w1", False)]),
    ]


def strand_selection_rows() -> list[Row]:
    return [
        Row(chunks=[("parental template       5′ -------- parental -------- 3′", None, False)]),
        Row(chunks=[("new BrdU strand         3′ --BrdU--BrdU--BrdU------ 5′", "w1", False)]),
        Row(chunks=[("                                      UV + Hoechst ↓", None, False)]),
        Row(chunks=[("retained library strand 5′ -------- parental -------- 3′", None, False)]),
    ]
