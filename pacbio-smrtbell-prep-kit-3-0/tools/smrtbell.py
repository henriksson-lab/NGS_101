"""Topology model for PacBio SMRTbell prep kit 3.0 libraries."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment
from dumbbell import (Dumbbell, smrtbell_circuit, smrtbell_hairpin_geometry,
                      smrtbell_sequencing_plan)
from endprep import dA_tailed_scene, repair_and_dA_tail
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind


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
FINAL_LIBRARY = None
SEQUENCING_UNAVAILABLE = ("The kit guide names the sequencing primer and SMRTbell "
                          "adapter but does not disclose their current bases.")
HAIRPIN_END = smrtbell_hairpin_geometry(
    SMRTBELL,
    statement="SMRTbell prep kit 3.0 ligates T-overhang hairpins to dA-tailed inserts.")


def end_prep_scene():
    return dA_tailed_scene(END_PREP)


def polymerase_cycle() -> Construct:
    return smrtbell_circuit(SMRTBELL)


FINAL_LIBRARY = polymerase_cycle()


def protocol_graph():
    authority = ('Commercial protocol: PacBio WGS/metagenome procedure '
                 '<a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-whole-genome-and-metagenome-libraries-using-SMRTbell-prep-kit-3.0.pdf">102-166-600 Rev09</a>.')
    g = ProtocolGraph("PacBio SMRTbell prep kit 3.0", authority=authority)
    dna = g.input("dna", "sheared 15–20 kb genomic DNA", MaterialKind.MOLECULAR,
                  state=INSERT, primary=True)
    repair = g.input("repair", "damage-repair, end-repair and dA-tailing enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed genomic inserts",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Repair damage and ends, phosphorylate and dA-tail",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    hairpins = g.input("hairpins", "T-overhang SMRTbell hairpin adapters",
                       MaterialKind.ADAPTER, inferred=True,
                       source=authority + " (adapter bases proprietary)")
    ligase = g.input("ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    closed = g.product("closed", "closed SMRTbell molecules", MaterialKind.MOLECULAR,
                       state=FINAL_LIBRARY, inferred=True,
                       source=authority + " (hairpin bases proprietary)")
    g.react("ligate", "Ligate a hairpin to each insert end", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (closed,), geometry=HAIRPIN_END)
    nuclease = g.input("nuclease", "free-end nuclease", MaterialKind.CATALYST)
    library = g.product("library", "nuclease-resistant closed SMRTbell library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=authority + " (hairpin bases proprietary)")
    linear = g.product("linear", "molecules retaining free DNA ends", MaterialKind.WASTE)
    g.react("remove_linear", "Digest molecules retaining a free end", StepKind.DIGESTION,
            ((closed, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (library, linear), changes_molecule=False)
    g.discard(linear, "not a completely closed SMRTbell")
    g.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair", "end_prep", "hairpins", "ligase",
                              "ligate", "nuclease", "remove_linear", "library"), authority)
    return g.require_publishable()
