"""Reusable TruSeq/Nextera molecular states and bounded workflow operations."""
from __future__ import annotations

from dataclasses import dataclass

import nextera as nx
import illumina as il
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      complement_segments)
from endprep import ThreePrimeOverhangAdapter

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

def _duplex(name: str, *, bases: str = "N" * 36) -> Scene:
    return Scene.duplex([Segment(name, bases, placeholder=True)], label=name)

def _single(name: str, parts: list[Segment], *, label: str) -> Scene:
    sc = Scene()
    sc.strand(name, parts, label=label)
    sc.labels(name)
    return sc

def _workflow(initial_name: str, initial: Scene,
              steps: tuple[tuple[str, str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, name, state, note in steps:
        w.react(action, state.rows(), name=name, note=note)
    return w


def nextera_final_rows(library: Construct) -> list:
    sc = Scene.duplex(list(library), label=library.name)
    insert = next(s.name for s in library if "insert" in s.name or "DNA" in s.name)
    sc.junction("top", "left mosaic end", insert, "transposition boundary")
    sc.junction("top", insert, "right mosaic end reverse complement",
                "transposition boundary")
    return sc.rows()

def _tagged_intermediate(label: str, *, gaps: bool = False) -> Scene:
    """Heterologous S5/S7 Tn5 product, drawn as the caller's source supports.

    ``gaps=True`` draws the transposition product before gap filling: each transferred
    adapter (entry point + mosaic end) is joined to a genomic 5' end and Tn5's 9-bp
    stagger leaves the nine genomic bases beside each join unpaired.  The
    non-transferred mosaic-end strand is not drawn, because vendor documents do not
    print it.  The legacy default (``gaps=False``) is a closed-duplex simplification
    that omits the gaps; a page using it must say so.
    """
    if not gaps:
        parts = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
                 Segment(label, "N" * 36, placeholder=True),
                 Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
                 Segment("S7 reverse complement", nx.S7_RC, "s7")]
        sc = Scene.duplex(parts, label="tagmented DNA")
        sc.junction("top", "left mosaic end", label, "transposition boundary")
        sc.junction("top", label, "right mosaic end reverse complement",
                    "transposition boundary")
        return sc
    gap = nx.TAGMENTATION_GAP
    core = Segment(label, "N" * 27, placeholder=True)
    left_gap = Segment("left 9-nt end", "N" * gap, placeholder=True)
    right_gap = Segment("right 9-nt end", "N" * gap, placeholder=True)
    top = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
           left_gap, core]
    bottom = [Segment("S7", nx.S7, "s7"), Segment("right mosaic end", nx.ME, "me"),
              right_gap, *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="tagmented DNA")
    sc.anneal("bottom", bottom, to="top", pair=(label + "'", label), label="tagmented DNA")
    sc.junction("top", "left mosaic end", left_gap.name, "transposition boundary")
    sc.junction("bottom", "right mosaic end", right_gap.name, "transposition boundary")
    sc.mark("top", left_gap.name, "unpaired: 9-nt gap opposite")
    sc.mark("bottom", right_gap.name, "unpaired: 9-nt gap opposite")
    sc.labels("top")
    return sc

@dataclass(frozen=True)
class Definition:
    title: str
    notes: str
    source: str
    summary: str
    caveat: str
    final_caption: str
    sequencing_intro: str
    final_kind: str
    insert_name: str
    workflow_builder: object
