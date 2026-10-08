"""Molecular model for the original Hsieh et al. Micro-C protocol."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments
from end_repair import fill_five_prime_overhang
from single_cell_hic import unresolved_illumina_library


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# A representative MNase-generated end, not a recognition sequence.  Both the repair
# strand and positions at which biotin-dA/dC are incorporated are computed from it.
REPAIR = fill_five_prime_overhang("TAGC", substituted_bases=("A", "C"))
SIZE_SELECTION_BP = (250, 350)
PCR_CYCLES = (12, 15)


def mnase_scene():
    top = [seg("left linker", "X" * 10, placeholder=True),
           seg("nucleosomal DNA", "N" * 34, "insert", placeholder=True),
           seg("right linker", "X" * 10, placeholder=True)]
    sc = Scene.duplex(top, label="crosslinked chromatin")
    sc.junction("top", "left linker", "nucleosomal DNA", "MNase", ch="^")
    sc.junction("top", "nucleosomal DNA", "right linker", "MNase", ch="^")
    sc.mark("top", "nucleosomal DNA", "protected mononucleosome")
    return sc


def resected_end_scene():
    top = [seg(f"overhang {i + 1}", b, "me")
           for i, b in enumerate(REPAIR.overhang_5p)]
    body = seg("nucleosomal DNA", "X" * 28, placeholder=True)
    sc = Scene(); sc.strand("top", [*top, body], label="resected end")
    sc.anneal("bottom", complement_segments([body]), to="top",
              pair=("nucleosomal DNA'", "nucleosomal DNA"), label="resected end")
    sc.mark("top", "overhang 1", "5′ single-stranded end", through="overhang 4")
    return sc


def repaired_end_scene():
    top = [seg(f"overhang {i + 1}", b, "me")
           for i, b in enumerate(REPAIR.overhang_5p)]
    body = seg("nucleosomal DNA", "X" * 28, placeholder=True)
    fill = [seg(f"fill {i + 1}", b, "w1" if b in REPAIR.substituted_bases else "me")
            for i, b in enumerate(REPAIR.fill_5p)]
    bottom = [*complement_segments([body]), *fill]
    sc = Scene(); sc.strand("top", [*top, body], label="blunt repaired end")
    sc.anneal("bottom", bottom, to="top", pair=("fill 1", "overhang 4"),
              label="blunt repaired end")
    for i in REPAIR.substituted_positions:
        sc.mark("bottom", f"fill {i + 1}", f"new biotin-d{REPAIR.fill_5p[i]}")
    return sc


def contact_junction():
    """Representative blunt ligation; unlike restriction Hi-C it has no fixed motif."""
    return Construct([
        seg("nucleosome A", "X" * 28, placeholder=True),
        seg("repaired end A", REPAIR.overhang_5p, "w1"),
        seg("repaired end B", "ACGT", "w1"),
        seg("nucleosome B", "X" * 28, placeholder=True),
    ], name="Micro-C proximity-ligation product")


def selected_insert():
    return Construct([
        seg("nucleosome A read end", "X" * 32, placeholder=True),
        *list(contact_junction().slice("repaired end A", "repaired end B")),
        seg("nucleosome B read end", "X" * 32, placeholder=True),
    ], name="biotin-selected Micro-C insert")


def final_library():
    return unresolved_illumina_library(selected_insert(),
        "Micro-C sequencing library", indexed=None)
