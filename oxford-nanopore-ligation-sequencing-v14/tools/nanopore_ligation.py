"""Molecular model for Oxford Nanopore SQK-LSK114 ligation sequencing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment
from endprep import (EndPreparedDuplex, ThreePrimeOverhangAdapter, dA_tailed_scene,
                     ligate_both_ends, cohesive_ligation_geometry)
from protocol_model import role
from protocol_graph import SequencingPlan, Terminal, Disposition

# ONT CHTD_500 v1 revAV, section15; source-owned declarations.
LA_TOP = "CCTGTACTTCGTTCAGTTACGTATTGCT"
LA_BOTTOM = "GCAATACGTAACTGAACGAAGTACAGG"
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


LIGATION_ADAPTER = ThreePrimeOverhangAdapter(
    "SQK-LSK114 Ligation Adapter disclosed duplex", LA_TOP, LA_BOTTOM
)
INSERT = Construct([seg("genomic DNA", "X" * 42, placeholder=True)], name="DNA fragment")
END_PREP = EndPreparedDuplex(INSERT, five_prime_phosphate=True, three_prime_overhang="A")
SEALED_JUNCTIONS = ligate_both_ends(END_PREP, LIGATION_ADAPTER)
END_PREP_TIMES = ((20, 5), (65, 5))
ADAPTER_END = cohesive_ligation_geometry(
    END_PREP, LIGATION_ADAPTER,
    statement="CHTD_500 revAV sections3/15: opposed insert3′A and LA3′T, "
              "with the disclosed LA strands paired computationally.")
# Adapter DNA assignment is conditional on the separately dated reference.
for _strand in ADAPTER_END.scenes[1].strands.values():
    for _part in _strand.segs:
        _part.inferred = True
TITLE = "Oxford Nanopore ligation sequencing V14 — SQK-LSK114"
NOTES = "01_nanopore-ligation-v14.html"
SUMMARY = "Repair and end-prepare genomic DNA, recover it on AMPure beads, ligate Ligation Adapter, recover the adapted library, and sequence with a nanopore motor."
CAVEAT = ("The workflow is GDE_9161_v114_revAC_24Sep2025; the disclosed DNA strands are "
          "a separately dated reference from CHTD_500_v1_revAV_24July2026. These DNA strands and their end arrangement are a later-reference conditional "
          "assignment, not proof of the2025 reagent lot or its complete motor-bearing assembly. The graph draws the "
          "disclosed DNA portion and a separate associated motor/leader boundary; motor "
          "attachment and strand-entry polarity remain unassigned in this scope.")
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "A motor-associated adapter captures DNA at a nanopore; no synthesis sequencing primer is used. Strand-entry polarity and read length remain molecule-dependent and unassigned here. Flow Cell Tether is supplied during flow-cell priming; no covalent DNA/tether connection is assigned."
SOURCE = ('Commercial protocol: Oxford Nanopore '
          '<a href="https://nanoporetech.com/resources/technique/whole-genome-sequencing/document/genomic-dna-by-ligation-sqk-lsk114">GDE_9161_v114 Rev AC</a>; later DNA-component reference: <a href="https://nanoporetech.com/document/chemistry-technical-document">CHTD_500 v1 revAV, section15</a>.')
MOTOR_LIMITATION = ("Nanopore sequencing is motor-driven; the full leader, motor "
                          "identity and exact motor-bearing strand remain unverified "
                          "in the inspected SQK-LSK114 sources.")


def end_prep_scene() -> Scene:
    return dA_tailed_scene(END_PREP)


def adapter_scene() -> Scene:
    top = [seg("27-bp adapter core", LIGATION_ADAPTER.core, "r1"),
           seg("3-prime dT", LIGATION_ADAPTER.overhang_3, "r2")]
    bottom = [seg("adapter bottom", LIGATION_ADAPTER.bottom, "r1")]
    for part in (*top,*bottom):
        part.inferred = True  # later sequence reference, not proof of2025lot identity
    sc = Scene(); sc.strand("adapter top", top, label="INFERRED later-reference LA top")
    sc.anneal("adapter bottom", bottom, to="adapter top",
              pair=("adapter bottom", "27-bp adapter core"), label="LA bottom")
    sc.mark("adapter top", "3-prime dT", "pairs the DNA dA")
    return sc


def final_library() -> Construct:
    left = SEALED_JUNCTIONS.left
    right = SEALED_JUNCTIONS.right
    return Construct([
        seg("left adapter core", LIGATION_ADAPTER.core, "r1", inferred=True),
        seg("left T:A junction", left.adapter_overhang, bottom=left.insert_overhang, inferred=True),
        seg("genomic DNA", "X" * 42, placeholder=True),
        seg("right A:T junction", right.insert_overhang, bottom=right.adapter_overhang, inferred=True),
        seg("right adapter core", LIGATION_ADAPTER.bottom, "r1", inferred=True),
    ], name="SQK-LSK114 adapted duplex")


def final_scene() -> Scene:
    lib = final_library()
    sc = Scene.duplex(list(lib), label="disclosed DNA portion")
    sc.strand("motor boundary", [role("unassigned motor/leader association", inferred=True)],
              label="associated complex; no covalent DNA placement")
    sc.junction("top", "left T:A junction", "genomic DNA", "sealed ligation")
    sc.junction("top", "genomic DNA", "right A:T junction", "sealed ligation")
    return sc


def protocol_graph():
    g = ProtocolGraph("Oxford Nanopore ligation sequencing V14", authority=SOURCE)
    dna = g.input("dna", "genomic DNA fragments", MaterialKind.MOLECULAR,
                  state=INSERT, primary=True)
    repair = g.input("repair", "NEBNext FFPE DNA Repair Mix (Buffer v2) and Ultra II "
                     "End-prep Enzyme Mix", MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired phosphorylated dA-tailed DNA",
                       MaterialKind.MOLECULAR, state=end_prep_scene())
    g.react("end_prep", "Repair damage, polish ends and add one 3-prime dA "
            "(20 °C 5 min, 65 °C 5 min)", StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    beads_1 = g.input("end_prep_beads", "60 µl AMPure XP Beads, 80% ethanol washes, "
                      "61 µl water elution", MaterialKind.SELECTION)
    tailed_clean = g.product("tailed_clean", "bead-purified dA-tailed DNA",
                             MaterialKind.MOLECULAR, state=end_prep_scene())
    end_prep_waste = g.product("end_prep_waste", "end-prep supernatant and ethanol washes",
                               MaterialKind.WASTE)
    g.react("end_prep_cleanup", "Clean up the end-prep reaction on AMPure XP beads",
            StepKind.CLEANUP,
            ((tailed, InputRole.TEMPLATE), (beads_1, InputRole.SELECTION)),
            (tailed_clean, end_prep_waste))
    g.discard(end_prep_waste, "unbound end-prep reaction components removed with the supernatant")
    qc_rule = g.input("qc_rule", "1 µl from61 µl eluate", MaterialKind.CONDITION)
    qc = g.product("qubit_aliquot", "1 µl for Qubit quantification", MaterialKind.MOLECULAR,
                   state=end_prep_scene())
    ligation_input = g.product("ligation_input", "60 µl end-prepped DNA for ligation",
                              MaterialKind.MOLECULAR, state=end_prep_scene())
    g.react("sample_qc", "Take 1 µl of the61 µl eluate for Qubit", StepKind.SAMPLING,
            ((tailed_clean, InputRole.TEMPLATE), (qc_rule, InputRole.CONDITION)), (qc, ligation_input), changes_molecule=False)
    g.mark_terminal(Terminal(qc.id, Disposition.UNRESOLVED,
                            "Qubit measurement aliquot; subsequent physical fate unassigned."))
    adapter = g.input("adapter", "motor-associated dT Ligation Adapter",
                      MaterialKind.ADAPTER, state=adapter_scene(), inferred=True,
                      source=SOURCE + " (disclosed DNA portion; full motor association unassigned)")
    ligase = g.input("ligase", "Ligation Buffer (LNB) and Salt-T4 DNA Ligase",
                     MaterialKind.CATALYST)
    residual = g.product("ligation_residual", "unligated/partially ligated DNA and other ligation products",
                         MaterialKind.MOLECULAR, state=tailed_clean.state, inferred=True)
    final = final_library()
    library = g.product("library", "motor-adapted ligation library",
                        MaterialKind.LIBRARY, state=final_scene(), inferred=True,
                        source=SOURCE + " (motor association unassigned)")
    g.react("ligate", "Ligate the motor-loaded dT adapter to both DNA ends "
            "(10 min, room temperature)", StepKind.LIGATION,
            ((ligation_input, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (library,residual), geometry=ADAPTER_END)
    g.mark_terminal(Terminal(residual.id, Disposition.UNRESOLVED,
                    "Ligation does not guarantee complete adaptation; later bead partition is component-specific and unassigned."))
    beads_2 = g.input("library_beads", "40 µl AMPure XP Beads, 250 µl LFB or SFB washes "
                      "(not ethanol), 15 µl Elution Buffer", MaterialKind.SELECTION)
    library_clean = g.product("library_clean", "bead-purified motor-adapted library",
                              MaterialKind.LIBRARY, state=final_scene(), inferred=True,
                              source=SOURCE + " (motor association unassigned)")
    library_waste = g.product("library_waste", "ligation supernatant (unligated Ligation "
                              "Adapter) and LFB/SFB washes", MaterialKind.WASTE)
    g.react("library_cleanup", "Clean up the adapted library; LFB enriches ≥3 kb, "
            "SFB keeps all sizes", StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (beads_2, InputRole.SELECTION)),
            (library_clean, library_waste))
    g.discard(library_waste, "unbound material removed with the bead supernatant and washes")
    g.sequence(library_clean, SequencingPlan(
        acquired_roles=("Nanopore strand",), read_uses={"Nanopore strand": "include"},
        limitation=MOTOR_LIMITATION))
    g.declare_reviewed_scope(("dna", "repair", "end_prep", "end_prep_beads",
                              "end_prep_cleanup", "qc_rule", "sample_qc", "adapter", "ligase", "ligate",
                              "library_beads", "library_cleanup", "library_clean"), SOURCE)
    return g.require_publishable()


def oligos():
    return [f"LA top (CHTD revAV): 5′-{LA_TOP}-3′",
            f"LA bottom (CHTD revAV): 5′-{LA_BOTTOM}-3′"]
