"""Shared exact chemistry for Spatial-CUT&Tag and Spatial-ATAC-seq.

The two assays differ in how Tn5 is targeted, but use the same deterministic
barcoding and sequencing-library architecture.  The oligos are from Farzad et
al., Nature Protocols 2024, Supplementary Tables 2 and 3.  The read-facing
orientation is derived here from those as-ordered oligos, so protocol pages do
not independently transcribe or reverse them.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Scene, Workflow, feature,
                      revcomp)


# As ordered, 5' -> 3'.  Asterisks in NEW_P5_SMART_WRITTEN are phosphorothioate
# bonds and therefore remain display metadata rather than bases in the model.
TN5_ME_REV_WRITTEN = "/5Phos/CTGTCTCTTATACACATCT"
TN5_ME_B = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
LIGATION_LINKER_1 = "AGTCGTACGCCGATGCTGTCTCTTATACACATCT"
LIGATION_LINKER_2 = "CGAATGCTCTGGCCTCTCAAGCACGTGGAT"
CUSTOM_READ_1 = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
NEW_P5_SMART = il.P5 + CUSTOM_READ_1
NEW_P5_SMART_WRITTEN = NEW_P5_SMART[:-3] + "T*A*C"
AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTG"

# Representative members of the two 50-oligo spatial plates.
BARCODE_A_1 = "AACGTGAT"
BARCODE_B_1 = "AACGTGAT"
BARCODE_A_OLIGO_1 = "/5Phos/AGGCCAGAGCATTCG" + BARCODE_A_1 + nx.ME
BARCODE_B_OLIGO_1 = CUSTOM_READ_1 + BARCODE_B_1 + "ATCCACGTGCTTGAG"

# The sequenced strand contains the reverse complements of the ordered splints.
# Split them at the experimentally ligated bonds.  Concatenation is checked at
# import time, making the exact junctions part of construction rather than prose.
LINKER_2_READ = revcomp(LIGATION_LINKER_2)
LINKER_2_B_SIDE = BARCODE_B_OLIGO_1[-15:]
LINKER_2_A_SIDE = BARCODE_A_OLIGO_1.removeprefix("/5Phos/")[:15]
LINKER_1_READ = revcomp(LIGATION_LINKER_1)
LINKER_1_A_SIDE = BARCODE_A_OLIGO_1.removeprefix("/5Phos/")[-len(nx.ME):]
LINKER_1_GENOMIC_SIDE = LINKER_1_READ[len(nx.ME):]
if LINKER_2_B_SIDE + LINKER_2_A_SIDE != LINKER_2_READ:
    raise ValueError("spatial linker 2 oligos do not reconstruct the read strand")
if LINKER_1_A_SIDE + LINKER_1_GENOMIC_SIDE != LINKER_1_READ:
    raise ValueError("spatial linker 1 oligos do not reconstruct the read strand")


SPATIAL_A = feature(
    "spatial_a", "spatial_barcode", "whitelist", group="spatial_coordinate",
    part="A", whitelist="Farzad et al. Supplementary Table 3, first barcodes")
SPATIAL_B = feature(
    "spatial_b", "spatial_barcode", "whitelist", group="spatial_coordinate",
    part="B", whitelist="Farzad et al. Supplementary Table 3, second barcodes")
SAMPLE_I7 = feature("sample_index_i7", "sample_index", "whitelist",
                    whitelist="Farzad et al. Supplementary Table 2, Ad2 panel")


def _right_arm():
    return [
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
    ]


def tagged_fragment(insert_name: str) -> Construct:
    """Representative linker-1/S7 product after Tn5 strand transfer."""
    return Construct([
        seg("ligation linker 1", LINKER_1_READ, "me"),
        seg(insert_name, "X" * 46, placeholder=True),
        *_right_arm(),
    ], name="linker-bearing tagmented fragment")


def a_barcoded_fragment(insert_name: str) -> Construct:
    """Product after barcode A is ligated at the linker-1 junction."""
    return Construct([
        seg("ligation linker 2 A-side", LINKER_2_A_SIDE),
        seg("spatial barcode A", "A" * 8, "cbc", placeholder=True,
            feature=SPATIAL_A),
        seg("ligation linker 1 A-side", LINKER_1_A_SIDE, "me"),
        seg("ligation linker 1 genomic-side", LINKER_1_GENOMIC_SIDE, "me"),
        seg(insert_name, "X" * 46, placeholder=True),
        *_right_arm(),
    ], name="A-barcoded tagmented fragment")


def ab_barcoded_fragment(insert_name: str) -> Construct:
    """Product after barcode B is ligated inside linker 2."""
    return Construct([
        seg("custom Read 1 site", CUSTOM_READ_1, "r1"),
        seg("spatial barcode B", "B" * 8, "cbc", placeholder=True,
            feature=SPATIAL_B),
        seg("ligation linker 2 B-side", LINKER_2_B_SIDE),
        *list(a_barcoded_fragment(insert_name)),
    ], name="A/B-barcoded tagmented fragment")


def library(insert_name: str):
    """PCR-completed, single-index library with every primer site validated."""
    core = list(ab_barcoded_fragment(insert_name))
    # The intermediate already carries the S7 end; PCR adds the P5 and indexed P7
    # flow-cell arms without changing the internal barcode/linker order.
    lib = Construct([
        seg("P5", il.P5, "p5"),
        *core,
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=SAMPLE_I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name=insert_name + " spatial epigenome library")
    primers = (
        sp.custom("Read 1", "Spatial epigenome custom Read 1", CUSTOM_READ_1,
                  "Farzad et al. 2024 Supplementary Tables 1-2"),
        sp.NEXTERA["I1"],
        sp.NEXTERA["R2"],
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid spatial epigenome library: " + "; ".join(problems))
    return lib, primers


def _scene(con: Construct, label: str, junctions=()):
    sc = Scene.duplex(list(con), label=label)
    for left, right, text in junctions:
        sc.junction("top", left, right, text)
    sc.labels("top")
    return sc


def tagged_rows(insert_name: str):
    return _scene(tagged_fragment(insert_name), "tagmented fragment").rows()


def a_barcoded_rows(insert_name: str):
    return _scene(a_barcoded_fragment(insert_name), "after barcode A", (
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def ab_barcoded_rows(insert_name: str):
    return _scene(ab_barcoded_fragment(insert_name), "after barcode B", (
        ("ligation linker 2 B-side", "ligation linker 2 A-side", "ligation"),
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def final_rows(lib: Construct):
    return _scene(lib, "library", (
        ("ligation linker 2 B-side", "ligation linker 2 A-side", "ligation"),
        ("ligation linker 1 A-side", "ligation linker 1 genomic-side", "ligation"),
    )).rows()


def workflow(insert_name: str, *, first_action: str, first_note: str) -> Workflow:
    """Shared deterministic-barcoding workflow with assay-specific targeting."""
    start = Scene.duplex(
        [seg("fixed tissue chromatin", "X" * 46, placeholder=True)],
        label="fixed tissue chromatin")
    wf = Workflow(MolecularState("Fixed tissue chromatin", tuple(start.rows())))
    wf.react(first_action, tagged_rows(insert_name), name="Linker-bearing fragment",
             note=first_note)
    wf.react("Flow barcode A through the first microchannel array",
             a_barcoded_rows(insert_name), name="A-barcoded fragment",
             note="The first 8-nt whitelist identifies one tissue coordinate; ligation occurs inside linker 1.")
    wf.react("Rotate the device and flow barcode B through perpendicular channels",
             ab_barcoded_rows(insert_name), name="A/B-barcoded fragment",
             note="The second 8-nt whitelist identifies the orthogonal coordinate; ligation occurs inside linker 2.")
    lib, _ = library(insert_name)
    wf.react("Reverse crosslinks and PCR-complete P5, indexed P7 and flow-cell arms",
             final_rows(lib), name="Sequencing library",
             note="The ordered B+A pair identifies the channel intersection. The construct uses one i7 sample index.")
    return wf


def oligos():
    """Copyable representative oligos, exactly as written in the protocol tables."""
    return (
        f"<b>Tn5MErev</b> — 5′-{TN5_ME_REV_WRITTEN}-3′",
        f"<b>Tn5ME-B</b> — 5′-{TN5_ME_B}-3′",
        f"<b>Ligation linker 1</b> — 5′-{LIGATION_LINKER_1}-3′",
        f"<b>Ligation linker 2</b> — 5′-{LIGATION_LINKER_2}-3′",
        f"<b>Tn5ME-A_BC1_1</b> — 5′-{BARCODE_A_OLIGO_1}-3′",
        f"<b>BC2_R1_1</b> — 5′-{BARCODE_B_OLIGO_1}-3′",
        f"<b>New-P5-SMART PCR hybrid oligo</b> — 5′-{NEW_P5_SMART_WRITTEN}-3′",
        f"<b>Ad2.1_TAAGGCGA</b> — 5′-{AD2_1}-3′",
        f"<b>Custom Read 1 primer</b> — 5′-{CUSTOM_READ_1}-3′",
    )
