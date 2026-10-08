"""Topology model for PacBio SMRTbell prep kit 3.0 libraries."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, revcomp
from dumbbell import Dumbbell
from endprep import dA_tailed_scene, repair_and_dA_tail


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


INSERT = Construct([seg("genomic insert", "X" * 48, placeholder=True)],
                   name="15–20 kb genomic insert")
END_PREP = repair_and_dA_tail(INSERT)
# Full current adapter bases are proprietary. N is molecularly unknown DNA, not a guessed
# sequence; the terminal T overhang is represented separately by Dumbbell.
HAIRPIN_PLACEHOLDER = "N" * 18
SMRTBELL = Dumbbell(INSERT, HAIRPIN_PLACEHOLDER, HAIRPIN_PLACEHOLDER,
                    name="SMRTbell library", insert_overhang="A", adapter_overhang="T")
END_PREP_TIMES = ((37, 30), (65, 5))
LIGATION = (20, 30)
NUCLEASE = (37, 15)


def end_prep_scene():
    return dA_tailed_scene(END_PREP)


def polymerase_cycle() -> Construct:
    con = Construct([
        seg("insert strand 1", INSERT.top(), placeholder=True),
        seg("right hairpin / primer site", HAIRPIN_PLACEHOLDER, "r1", placeholder=True),
        seg("insert strand 2 (reverse complement)", revcomp(INSERT.top()),
            placeholder=True),
        seg("left hairpin / primer site", HAIRPIN_PLACEHOLDER, "r1", placeholder=True),
    ], name="one SMRTbell polymerase circuit")
    if con.top() != SMRTBELL.pass_sequence:
        raise AssertionError("polymerase circuit must be the covalently closed dumbbell path")
    return con
