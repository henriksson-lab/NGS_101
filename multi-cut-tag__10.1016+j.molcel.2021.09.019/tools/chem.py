"""Source-backed bulk Multi-CUT&Tag molecular model.

The defining paper also adapts the chemistry to the 10x scATAC workflow. The exact
antibody-coded transposon and bulk sequencing geometry are printed in Table S1, whereas
the vendor cell-barcode oligos are not; this module therefore draws the published bulk
construct exactly and keeps the single-cell branch explicit but separate.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Row, Scene, Workflow,
                      duplex_rows, feature, revcomp)
from tagmentation import loaded_tn5_adapter

TITLE = "Multi-CUT&Tag — antibody-coded targeted transposition"
NOTES = "01_multi-cut-tag.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.molcel.2021.09.019">'
          'Gopalan et al., <i>Molecular Cell</i> (2021)</a>.')
SUMMARY = ("Each antibody is coupled to pA–Tn5 loaded with a distinct eight-base adapter "
           "code. Both ends of every chromatin fragment consequently report which "
           "antibody-targeted transposomes made it.")
CAVEAT = ("The exact construct below is the paper’s bulk Multi-CUT&Tag branch. The "
          "single-cell branch passes the same target-coded fragments into the 10x Chromium "
          "Single Cell ATAC workflow; its proprietary gel-bead oligo sequence is not printed "
          "in the paper and is therefore not invented here.")

# Table S1. The four entries in each family differ only in the antibody barcode.
READ1 = "TCGTCGGCAGCGTCTCCACGC"
READ2 = "GTCTCGTGGGCTCGGCTGTCCCTGTCC"
INDEX1 = "GGACAGGGACAGCCGAGCCCACGAGAC"
INDEX2 = "GCGTGGAGACGCTGCCGACGA"
CONNECTOR_A = "GCGATCGAGGACGGC"
CONNECTOR_B = "CACCGTCTCCGCCTC"
BARCODES_A = ("TATAGCCT", "ATAGAGGC", "CCTATCCT", "GGCTCTGA")
BARCODES_B = ("CGAGTAAT", "TCTCCGGA", "AATGAGCG", "GGAATCTC")
TN5_ME_REV = nx.ME_RC


def _target_feature(side: str):
    return feature(f"antibody_barcode_{side.lower()}", "feature_barcode", "fixed",
                   group="antibody_identity", part=side,
                   note="one of four Table S1 sequences")


def adapter_a(barcode: str = BARCODES_A[0]):
    if barcode not in BARCODES_A:
        raise ValueError("adapter A needs a published antibody barcode")
    return [seg("Read 1 priming site", READ1, "r1"),
            seg("A-side antibody barcode", barcode, "cbc",
                feature=_target_feature("A")),
            seg("connector A", CONNECTOR_A, "me"),
            seg("A mosaic end", nx.ME, "me")]


def adapter_b(barcode: str = BARCODES_B[1]):
    if barcode not in BARCODES_B:
        raise ValueError("adapter B needs a published antibody barcode")
    return [seg("Read 2 priming site", READ2, "r2"),
            seg("B-side antibody barcode", barcode, "cbc",
                feature=_target_feature("B")),
            seg("connector B", CONNECTOR_B, "me"),
            seg("B mosaic end", nx.ME, "me")]


def loaded_adapter_scene(parts, label: str) -> Scene:
    """A source-listed barcoded top oligo annealed to phosphorylated Tn5MErev."""
    sc = loaded_tn5_adapter(
        parts, mosaic_end=parts[-1].name, label=label,
        bottom_name="Tn5MErev", bottom_label="5′-phosphorylated Tn5MErev",
        bottom_sequence=TN5_ME_REV, mod5="p")
    barcode = next(s.name for s in parts if s.feature)
    sc.mark("Tn5 top", barcode, "8-nt antibody identity")
    return sc


def transferred_fragment() -> Construct:
    """Representative heterologous A/B product after gap fill.

    Different barcode choices at the two ends are deliberately independent: equal codes
    are homogeneous reads, while unequal codes are the paper's mixed reads.
    """
    left = adapter_a(BARCODES_A[0])
    right = adapter_b(BARCODES_B[1])
    right_rc = [seg(s.name + " reverse complement", revcomp(s.top), s.tag,
                    feature=s.feature) for s in reversed(right)]
    return Construct([*left, seg("target-proximal chromatin", "X" * 36,
                                 placeholder=True), *right_rc],
                     name="gap-filled antibody-coded fragment")


TRANSFERRED = transferred_fragment()


def _final_library() -> Construct:
    return Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "J" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i5", "sample_index", "unknown")),
        *list(TRANSFERRED),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("sample_index_i7", "sample_index", "unknown")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="bulk Multi-CUT&Tag library")


FINAL_LIBRARY = _final_library()
SEQ_PRIMERS = (
    sp.custom("Read 1", "Multi-CUT&Tag Read1", READ1, "Gopalan et al. Table S1"),
    sp.custom("Index 1 (i7)", "Multi-CUT&Tag Index1", INDEX1,
              "Gopalan et al. Table S1"),
    sp.custom("Index 2 (i5)", "Multi-CUT&Tag Index2", INDEX2,
              "Gopalan et al. Table S1"),
    sp.custom("Read 2", "Multi-CUT&Tag Read2", READ2, "Gopalan et al. Table S1"),
)
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS)
if _problems:
    raise ValueError("Multi-CUT&Tag final library: " + "; ".join(_problems))

FINAL_CAPTION = ("Representative bulk A/B library. The independently chosen A- and B-side "
                 "codes may identify the same antibody (homogeneous) or two antibodies "
                 "occupying the same fragment (mixed).")
SEQUENCING_INTRO = ("The four custom Table S1 primers land by exact sequence. Each genomic "
                    "read first reports an 8-nt antibody barcode, 15-nt connector and 19-nt "
                    "mosaic end: the 42 bases removed before genomic alignment.")
READ_LENGTHS = {"Read 1": 72, "Index 1 (i7)": 8,
                "Index 2 (i5)": 8, "Read 2": 72}
READ_USES = {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "extract_then_include"}


def oligos():
    rows = []
    for i, bc in enumerate(BARCODES_A, 1):
        rows.append(f"<b>P5_i5_{i}_Universal_Connector_A</b> — 5′-{READ1}{bc}{CONNECTOR_A}{nx.ME}-3′")
    for i, bc in enumerate(BARCODES_B, 1):
        rows.append(f"<b>P7_i7_{i}_Universal_Connector_B</b> — 5′-{READ2}{bc}{CONNECTOR_B}{nx.ME}-3′")
    rows.extend((
        f"<b>Tn5MErev</b> — 5′-/5Phos/{TN5_ME_REV}-3′",
        f"<b>Read1 / Read2</b> — 5′-{READ1}-3′ / 5′-{READ2}-3′",
        f"<b>Index1 / Index2</b> — 5′-{INDEX1}-3′ / 5′-{INDEX2}-3′",
    ))
    return rows


def workflow() -> Workflow:
    start = Scene.duplex([seg("chromatin", "X" * 42, placeholder=True)],
                         label="native chromatin")
    wf = Workflow(MolecularState("Antibody-bound native chromatin", tuple(start.rows())))
    wf.react("Anneal the two target-coded adapter families and load pA–Tn5",
             [*loaded_adapter_scene(adapter_a(), "adapter A").rows(),
              *loaded_adapter_scene(adapter_b(), "adapter B").rows()],
             name="Antibody-specific A/B transposomes",
             note="Each antibody receives one matched A/B barcode pair; Tn5MErev is 5′-phosphorylated.")
    wf.react("Tether the loaded pA–Tn5 complexes through their antibodies",
             [Row(chunks=[("chromatin target ← antibody ← pA–Tn5 [exact A/B duplexes above]", None, False)])],
             name="Washed, target-bound transposomes")
    wf.react("Activate tagmentation and fill the 9-bp staggered gaps",
             duplex_rows(TRANSFERRED, label="target-coded chromatin fragment"),
             name="Gap-filled antibody-coded fragments",
             note="The representative has different end codes; same-code products are equally valid.")
    wf.react("Amplify with indexed i5/i7 primers",
             duplex_rows(FINAL_LIBRARY, label="bulk sequencing library"),
             name="Bulk Multi-CUT&Tag library",
             note="PCR adds P5, P7 and the two 8-cycle sample indexes around the exact target-coded insert.")
    wf.react("For the single-cell branch, instead enter the 10x scATAC GEM workflow",
             [Row(chunks=[("same target-coded fragment → 10x gel-bead cell barcode → sample-index PCR", None, False)])],
             name="Single-cell branch (vendor oligo sequence unavailable)",
             note="The paper uses 16 Index 2 cycles for the single-cell run; the proprietary gel-bead construct is not printed.")
    return wf
