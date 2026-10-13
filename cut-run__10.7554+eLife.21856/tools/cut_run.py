"""Original Skene & Henikoff CUT&RUN molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import duplex_fragment, seg, targeting_rows
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "CUT&RUN — antibody-tethered MNase chromatin profiling"
NOTES = "01_cut-run.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.7554/eLife.21856">Skene & Henikoff, <i>eLife</i> (2017)</a>.'
SUMMARY = "An antibody recruits Protein A–MNase to native chromatin. Calcium activates local cleavage; fragments cut on both sides of the bound particle diffuse into the supernatant and receive conventional Illumina adapters afterward."
CAVEAT = "CUT&RUN does not install adapters during targeted cleavage. The defining paper used the KAPA DNA polymerase library-preparation chemistry, but KAPA supplied no adapters in that kit and the paper does not identify the separately chosen adapter/index oligos. Their sequences, index count and sequencing-primer set are therefore left unavailable instead of being replaced by generic TruSeq."

FINAL_LIBRARY = Construct([
    seg("unreported left adapter", "[adapter]", placeholder=True, inferred=True),
    seg("released target-proximal DNA", "X" * 38, placeholder=True),
    seg("unreported right adapter / index", "[adapter/index]", placeholder=True,
        inferred=True, feature=feature("sample_index_region", "sample_index", "unknown",
                                      note="index count and sub-boundary are not reported")),
], name="CUT&RUN library (adapter identity unavailable)")
SEQ_PRIMERS = ()
FINAL_CAPTION = "INFERRED — the released MNase fragment is end-repaired, dA-tailed, adapter-ligated and PCR-amplified. Dotted role tokens mark adapter information that the defining paper does not report."
SEQUENCING_ENDING = "The paper reports 25-cycle paired-end HiSeq 2500 sequencing, but does not name the adapter oligos or sequencing primers."
SEQUENCING_UNAVAILABLE = "adapter and sequencing-primer identities are absent from the defining paper"

def sections():
    missing=AdapterLigationGeometry.unavailable("Skene and Henikoff specify KAPA end repair/dA-tailing but do not identify the separately supplied adapter/index oligos; the complementary adapter overhang is therefore unavailable.")
    return [
        ("Immobilize nuclei and tether pA–MNase",
         targeting_rows("MNase", "inactive nuclease positioned at the epitope"),
         "Concanavalin-A beads retain permeabilized nuclei during antibody and pA–MNase binding and washing."),
        ("Activate cleavage and release the selected fragment",
         [Row(chunks=[("chromatin ---- ^ [antibody-bound particle] ^ ---- chromatin", None, False)]),
          Row(chunks=[("              Ca²⁺-activated cuts; two-cut fragment enters supernatant", "me", False)])],
         "Calcium activates MNase at 0 °C. Chelators stop digestion; fragments cleaved on both sides of a particle are recovered from the soluble fraction."),
        ("Recover DNA for conventional library preparation",
         duplex_fragment("released MNase fragment").rows()+missing.rows(),
         "Proteinase treatment and extraction yield an unadapted duplex. Adapter ligation is a separate downstream operation, unlike CUT&Tag."),
    ]

def workflow():
    s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2])
    g=AdapterLigationGeometry.unavailable("Skene and Henikoff specify KAPA end repair/dA-tailing but do not identify the separately supplied adapter/index oligos; the complementary adapter overhang is therefore unavailable.")
    w.ligate_adapter("Attach the unreported library adapters",g,s[2][1],note=s[2][2]);return w


def protocol_graph() -> ProtocolGraph:
    """Material graph through targeted cleavage and source-limited library prep."""
    source = "Skene and Henikoff 2017, doi:10.7554/eLife.21856"
    graph = ProtocolGraph(TITLE, authority=source)

    def add_input(key, label, kind, *, state=None, primary=False):
        value = Material(key, label, kind, state=state, primary=primary, source=source)
        graph.add_input(value)
        return value

    nuclei = add_input("nuclei", "permeabilized native-chromatin nuclei",
                       MaterialKind.BIOLOGICAL, state=sections()[0][1], primary=True)
    cona = add_input("cona_beads", "concanavalin-A magnetic beads",
                     MaterialKind.SELECTION)
    antibody = add_input("primary_antibody", "target-specific primary antibody",
                         MaterialKind.SELECTION)
    pamnase = add_input("pa_mnase", "Protein A–MNase", MaterialKind.CATALYST)
    tethered = Material("tethered_nuclei", "washed antibody-tethered pA–MNase nuclei",
                        MaterialKind.MOLECULAR, state=sections()[0][1], source=source)
    wash_waste = Material("unbound_targeting_reagents", "unbound targeting reagents",
                          MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "target_chromatin", "Immobilize nuclei and tether pA–MNase",
        StepKind.AFFINITY_CAPTURE,
        (Input(nuclei.id, InputRole.TEMPLATE), Input(cona.id, InputRole.SELECTION),
         Input(antibody.id, InputRole.SELECTION),
         Input(pamnase.id, InputRole.CATALYST)),
        (tethered.id, wash_waste.id), source=source), (tethered, wash_waste))
    graph.mark_terminal(Terminal(
        wash_waste.id, Disposition.DISCARD,
        "removed during washes before calcium activation"))

    calcium = add_input("calcium", "calcium activation condition", MaterialKind.CONDITION)
    released = Material("released_fragment", "released target-proximal MNase fragment",
                        MaterialKind.MOLECULAR, state=duplex_fragment(
                            "released MNase fragment"), source=source)
    residual = Material("residual_nuclei", "bead-bound residual nuclei and chromatin",
                        MaterialKind.WASTE, source=source)
    graph.add_step(Step(
        "targeted_cleavage", "Activate MNase and release selected fragments",
        StepKind.DIGESTION,
        (Input(tethered.id, InputRole.TEMPLATE),
         Input(pamnase.id, InputRole.CATALYST),
         Input(calcium.id, InputRole.CONDITION)),
        (released.id, residual.id), source=source), (released, residual))
    graph.mark_terminal(Terminal(
        residual.id, Disposition.DISCARD,
        "the soluble two-cut fragments are recovered; bead-bound nuclei are not sequenced"))

    proteinase = add_input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    recovered = Material("recovered_dna", "deproteinized CUT&RUN DNA",
                         MaterialKind.MOLECULAR,
                         state=duplex_fragment("released MNase fragment"), source=source)
    graph.add_step(Step(
        "recover_dna", "Digest protein and recover released DNA",
        StepKind.DEPROTEINIZATION,
        (Input(released.id, InputRole.TEMPLATE),
         Input(proteinase.id, InputRole.CATALYST)),
        (recovered.id,), source=source), (recovered,))

    repair = add_input("end_repair_enzymes", "KAPA end-repair enzymes",
                       MaterialKind.CATALYST)
    repaired = Material("end_repaired_dna", "end-repaired CUT&RUN DNA",
                        MaterialKind.MOLECULAR,
                        state=duplex_fragment("end-repaired MNase fragment"), source=source)
    graph.add_step(Step(
        "end_repair", "End-repair recovered DNA", StepKind.END_REPAIR,
        (Input(recovered.id, InputRole.TEMPLATE), Input(repair.id, InputRole.CATALYST)),
        (repaired.id,), source=source), (repaired,))

    tailing = add_input("a_tailing_enzyme", "KAPA dA-tailing enzyme",
                        MaterialKind.CATALYST)
    da_scene = Scene.duplex([
        Segment("left 3′ dA", "A"),
        Segment("released MNase fragment", "X" * 34, placeholder=True),
        Segment("right 3′ dA", "A"),
    ], label="dA-tailed CUT&RUN DNA")
    tailed = Material("da_tailed_dna", "dA-tailed CUT&RUN DNA",
                      MaterialKind.MOLECULAR, state=da_scene, source=source)
    graph.add_step(Step(
        "da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
        (Input(repaired.id, InputRole.TEMPLATE),
         Input(tailing.id, InputRole.CATALYST)),
        (tailed.id,), source=source), (tailed,))

    adapters = add_input("unreported_adapters", "separately supplied adapter/index oligos",
                         MaterialKind.ADAPTER, state=FINAL_LIBRARY,)
    ligase = add_input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = Material("adapter_ligated_dna", "adapter-ligated CUT&RUN DNA",
                       MaterialKind.MOLECULAR, state=FINAL_LIBRARY, source=source,
                       inferred=True)
    missing = AdapterLigationGeometry.unavailable(
        "The paper specifies KAPA end repair/dA-tailing but does not identify the separately supplied adapter/index oligos.")
    graph.add_step(Step(
        "adapter_ligation", "Attach the unreported library adapters",
        StepKind.LIGATION,
        (Input(tailed.id, InputRole.TEMPLATE),
         Input(adapters.id, InputRole.STRUCTURAL),
         Input(ligase.id, InputRole.CATALYST)),
        (ligated.id,), geometry=missing, source=source), (ligated,))

    forward = add_input("unreported_forward_primer", "unreported PCR primer 1",
                        MaterialKind.OLIGO)
    reverse = add_input("unreported_reverse_primer", "unreported PCR primer 2",
                        MaterialKind.OLIGO)
    library = Material("cut_run_library", "CUT&RUN sequencing library",
                       MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source,
                       inferred=True)
    graph.add_step(Step(
        "library_pcr", "PCR-amplify the ligated library", StepKind.PCR,
        (Input(ligated.id, InputRole.TEMPLATE),
         Input(forward.id, InputRole.STRUCTURAL),
         Input(reverse.id, InputRole.STRUCTURAL)),
        (library.id,), geometry=PCRGeometry.unavailable(
            FINAL_LIBRARY,
            reason="The defining paper does not identify the separately supplied PCR primer oligos."),
        source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        limitation="The defining paper reports paired-end HiSeq sequencing but not the adapter or sequencing-primer identities."))

    graph.declare_reviewed_scope((
        "nuclei", "cona_beads", "primary_antibody", "pa_mnase",
        "target_chromatin", "calcium", "targeted_cleavage", "proteinase_k",
        "recover_dna", "end_repair_enzymes", "end_repair", "a_tailing_enzyme",
        "da_tail", "unreported_adapters", "dna_ligase", "adapter_ligation",
        "unreported_forward_primer", "unreported_reverse_primer", "library_pcr",
        "cut_run_library",
    ), source)
    return graph.require_publishable()
