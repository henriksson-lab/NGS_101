"""Reusable molecular views for multimodal and spatial barcode-transfer protocols."""
from __future__ import annotations

from chemdraw import Row, Scene, Segment, feature


def barcode_rounds(rounds: int, *, label: str = "cell barcode") -> list[Row]:
    """Draw split-pool barcode ligations, with every covalent junction explicit."""
    if rounds < 1:
        raise ValueError("at least one barcode round is required")
    text = "molecule"
    for n in range(1, rounds + 1):
        text += f" ** [{label} {n}]"
    return [Row(chunks=[(text, "cbc", False)]),
            Row(chunks=[("** = ligation at the named barcode boundary", "me", False)])]


def modality_split(*names: str) -> list[Row]:
    """Show one tagged cell yielding separately amplified sequencing libraries."""
    if len(names) < 2:
        raise ValueError("a modality split needs at least two libraries")
    rows = [Row(chunks=[("one cell identity", "cbc", False)])]
    for i, name in enumerate(names):
        rows.append(Row(chunks=[(("├─ " if i < len(names) - 1 else "└─ ") + name,
                                 None, False)]))
    return rows


def adjacent_probe_scene() -> Scene:
    """Adjacent probe halves on fixed RNA; Scene enforces the duplex geometry."""
    target = [Segment("target RNA left", "X" * 18, placeholder=True),
              Segment("target RNA right", "X" * 18, placeholder=True)]
    left = [Segment("left probe barcode", "L" * 12, "cbc", placeholder=True,
                    feature=feature("probe_barcode_left", "feature_barcode", "unknown",
                                    group="probe_barcode", part="left")),
            Segment("left probe target arm", "X" * 18, placeholder=True)]
    right = [Segment("right probe target arm", "X" * 18, placeholder=True),
             Segment("right probe barcode", "R" * 12, "cbc", placeholder=True,
                     feature=feature("probe_barcode_right", "feature_barcode", "unknown",
                                     group="probe_barcode", part="right"))]
    sc = Scene(); sc.strand("RNA", target, label="fixed target RNA")
    sc.anneal("left probe", left, to="RNA", pair=("left probe target arm", "target RNA left"),
              label="left probe", unpaired=("left probe barcode",))
    sc.anneal("right probe", right, to="RNA", pair=("right probe target arm", "target RNA right"),
              label="right probe", unpaired=("right probe barcode",))
    sc.arrow("left probe", "adjacent probe-pair ligation")
    sc.labels("left probe"); sc.labels("right probe")
    return sc


def surface_capture_rows(surface: str, barcode: str, molecule: str) -> list[Row]:
    return [
        Row(chunks=[(f"{surface}:  ●—[{barcode}]—[UMI]—poly(dT)", "cbc", False)]),
        Row(chunks=[(f"                                      || {molecule} poly(A)", None, False)]),
        Row(chunks=[("                                      ← reverse transcription", "me", False)]),
    ]
