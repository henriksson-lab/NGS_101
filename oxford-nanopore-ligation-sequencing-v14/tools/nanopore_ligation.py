"""Molecular model for Oxford Nanopore SQK-LSK114 ligation sequencing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment
from endprep import (ThreePrimeOverhangAdapter, dA_tailed_scene, ligate_both_ends,
                     repair_and_dA_tail)
from nanopore_library import LA_BOTTOM, LA_TOP, motor_sequencing_plan
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


LIGATION_ADAPTER = ThreePrimeOverhangAdapter(
    "SQK-LSK114 Ligation Adapter disclosed duplex", LA_TOP, LA_BOTTOM
)
INSERT = Construct([seg("genomic DNA", "X" * 42, placeholder=True)], name="DNA fragment")
END_PREP = repair_and_dA_tail(INSERT)
SEALED_JUNCTIONS = ligate_both_ends(END_PREP, LIGATION_ADAPTER)
END_PREP_TIMES = ((20, 5), (65, 5))
ADAPTER_END = reported_adapter_end(
    "SQK-LSK114 Ligation Adapter", kind="cohesive", polarity="3-prime",
    terminal_sequence="T", body_sequence=LA_TOP,
    statement="The disclosed LA duplex has a 27-bp core and one 3-prime dT overhang.")
SOURCE = ('Commercial protocol: Oxford Nanopore '
          '<a href="https://nanoporetech.com/resources/technique/whole-genome-sequencing/document/genomic-dna-by-ligation-sqk-lsk114">GDE_9161_v114 Rev AC</a>.')
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing is motor-driven; the full leader, motor "
                          "identity and exact motor-bearing strand are proprietary.")


def end_prep_scene() -> Scene:
    return dA_tailed_scene(END_PREP)


def adapter_scene() -> Scene:
    top = [seg("27-bp adapter core", LIGATION_ADAPTER.core, "r1"),
           seg("3-prime dT", LIGATION_ADAPTER.overhang_3, "r2")]
    bottom = [seg("adapter bottom", LIGATION_ADAPTER.bottom, "r1")]
    sc = Scene(); sc.strand("adapter top", top, label="LA top")
    sc.anneal("adapter bottom", bottom, to="adapter top",
              pair=("adapter bottom", "27-bp adapter core"), label="LA bottom")
    sc.mark("adapter top", "3-prime dT", "pairs the DNA dA")
    return sc


def final_library() -> Construct:
    left = SEALED_JUNCTIONS.left
    right = SEALED_JUNCTIONS.right
    return Construct([
        seg("left motor-loaded leader", "[motor-loaded leader]", "r2", placeholder=True),
        seg("left adapter core", LIGATION_ADAPTER.core, "r1"),
        seg("left T:A junction", left.adapter_overhang, bottom=left.insert_overhang),
        seg("genomic DNA", "X" * 42, placeholder=True),
        seg("right A:T junction", right.insert_overhang, bottom=right.adapter_overhang),
        seg("right adapter core", LIGATION_ADAPTER.bottom, "r1"),
        seg("right motor-loaded leader", "[motor-loaded leader]", "r2", placeholder=True),
    ], name="SQK-LSK114 adapted duplex")


def final_scene() -> Scene:
    lib = final_library()
    sc = Scene.duplex(list(lib), label="adapted DNA",
                      unpaired=("left motor-loaded leader'", "right motor-loaded leader'"))
    sc.junction("top", "left T:A junction", "genomic DNA", "sealed ligation")
    sc.junction("top", "genomic DNA", "right A:T junction", "sealed ligation")
    return sc


def protocol_graph():
    g = ProtocolGraph("Oxford Nanopore ligation sequencing V14", authority=SOURCE)
    dna = g.input("dna", "genomic DNA fragments", MaterialKind.MOLECULAR,
                  state=INSERT, primary=True)
    repair = g.input("repair", "damage-repair and end-prep enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired phosphorylated dA-tailed DNA",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Repair damage, polish ends and add one 3-prime dA",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    adapter = g.input("adapter", "motor-loaded dT Ligation Adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (LA core disclosed; full motor leader proprietary)")
    ligase = g.input("ligase", "Salt-T4 ligase", MaterialKind.CATALYST)
    final = final_library()
    library = g.product("library", "motor-adapted ligation library",
                        MaterialKind.LIBRARY, state=final, inferred=True,
                        source=SOURCE + " (motor leader proprietary)")
    g.react("ligate", "Ligate the motor-loaded dT adapter to both DNA ends",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (library,), geometry=ADAPTER_END)
    g.sequence(library, motor_sequencing_plan(
        final, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair", "end_prep", "adapter", "ligase",
                              "ligate", "library"), SOURCE)
    return g.require_publishable()
