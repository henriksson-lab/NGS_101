"""Original Skene & Henikoff CUT&RUN molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
from batch_ngs import duplex_fragment, seg, targeting_rows, truseq_library
from chemdraw import Construct, MolecularState, Row, Scene, Workflow
from endprep import (cohesive_ligation_geometry, dA_tailed_scene,
                     ligation_product_scene, repair_and_dA_tail)
from illumina_vendor_workflows import TRUSEQ_Y_END
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "CUT&RUN — antibody-tethered MNase chromatin profiling"
NOTES = "01_cut-run.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.7554/eLife.21856">Skene & Henikoff, <i>eLife</i> (2017)</a>.'
SUMMARY = "An antibody recruits Protein A–MNase to native chromatin. Calcium activates local cleavage; fragments cut on both sides of the bound particle diffuse into the supernatant and receive conventional Illumina adapters afterward."
CAVEAT = ("CUT&RUN does not install adapters during targeted cleavage. Its Methods point "
          "to the same library preparation used by Kasinathan et al.; that primary "
          "supplement explicitly identifies Illumina TruSeq adapters. The selected "
          "sample-index member is not reported, so its eight bases remain a placeholder, "
          "but the adapter family, ligation end and standard sequencing-primer sites are known.")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("released target-proximal DNA", "X" * 38, placeholder=True),
], "single-index TruSeq CUT&RUN library", dual_index=False)
FINAL_CAPTION = ("PCR-completed single-index TruSeq library. The i7 position is exact, "
                 "but the particular eight-base sample index selected for this experiment "
                 "was not reported.")
SEQUENCING_INTRO = ("Standard TruSeq primers acquire paired reads from the two ends of "
                    "the released CUT&RUN fragment and one i7 sample-index read.")


END_PREP = repair_and_dA_tail(Construct([
    seg("released target-proximal DNA", "X" * 38, placeholder=True),
], name="released CUT&RUN fragment"))


def truseq_ligation_geometry():
    return cohesive_ligation_geometry(
        END_PREP, TRUSEQ_Y_END,
        statement=("The CUT&RUN paper reuses the Kasinathan et al. library procedure; "
                   "its supplement names Illumina TruSeq adapters, whose public "
                   "3-prime-dT end is complementary to the dA-tailed fragment."))


def ligated_library_scene() -> Scene:
    return ligation_product_scene(
        FINAL_LIBRARY, "released target-proximal DNA",
        label="TruSeq-adapter-ligated CUT&RUN molecule",
        adapter_group="TruSeq adapter")

def sections():
    return [
        ("Immobilize nuclei and tether pA–MNase",
         targeting_rows("MNase", "inactive nuclease positioned at the epitope"),
         "Concanavalin-A beads retain permeabilized nuclei during antibody and pA–MNase binding and washing."),
        ("Activate cleavage and release the selected fragment",
         [Row(chunks=[("chromatin ---- ^ [antibody-bound particle] ^ ---- chromatin", None, False)]),
          Row(chunks=[("              Ca²⁺-activated cuts; two-cut fragment enters supernatant", "me", False)])],
         "Calcium activates MNase at 0 °C. Chelators stop digestion; fragments cleaved on both sides of a particle are recovered from the soluble fraction."),
        ("Recover DNA for conventional library preparation",
         duplex_fragment("released MNase fragment").rows(),
         "Proteinase treatment and extraction yield an unadapted duplex. Adapter ligation is a separate downstream operation, unlike CUT&Tag."),
    ]

def workflow():
    s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2])
    w.react(s[2][0], s[2][1], note=s[2][2])
    w.react("End-repair and add one 3-prime dA to each end",
            dA_tailed_scene(END_PREP, "dA-tailed CUT&RUN DNA").rows(),
            note="KAPA enzymes prepare each released fragment for cohesive adapter ligation.")
    w.ligate_adapter(
        "Ligate Illumina TruSeq adapters", truseq_ligation_geometry(),
        ligated_library_scene().rows(),
        note="The cited primary library workflow identifies the adapter family; the two adapter copies share one reagent identity.")
    w.pcr(
        "PCR-enrich the TruSeq library",
        PCRGeometry.regional(
            FINAL_LIBRARY,
            PCRPrimer("TruSeq universal PCR primer", il.TRUSEQ_P5_FULL),
            reverse_pool="selected TruSeq indexed PCR primer",
            reverse_region="Index 1 / Read 2 arm",
            statement=("The workflow uses TruSeq library amplification, but the paper "
                       "does not identify the selected eight-base index member.")),
        ligated_library_scene().rows(),
        note="The universal left primer lands exactly; the indexed right-primer member is source-limited to its TruSeq arm and pool.")
    return w


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
    da_scene = dA_tailed_scene(END_PREP, "dA-tailed CUT&RUN DNA")
    tailed = Material("da_tailed_dna", "dA-tailed CUT&RUN DNA",
                      MaterialKind.MOLECULAR, state=da_scene, source=source)
    graph.add_step(Step(
        "da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
        (Input(repaired.id, InputRole.TEMPLATE),
         Input(tailing.id, InputRole.CATALYST)),
        (tailed.id,), source=source), (tailed,))

    adapters = add_input("truseq_adapters", "Illumina TruSeq adapter oligos",
                         MaterialKind.ADAPTER, state=TRUSEQ_Y_END.scene())
    ligase = add_input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = Material("adapter_ligated_dna", "TruSeq-adapter-ligated CUT&RUN DNA",
                       MaterialKind.MOLECULAR, state=ligated_library_scene(), source=source)
    graph.add_step(Step(
        "adapter_ligation", "Ligate Illumina TruSeq adapters",
        StepKind.LIGATION,
        (Input(tailed.id, InputRole.TEMPLATE),
         Input(adapters.id, InputRole.STRUCTURAL),
         Input(ligase.id, InputRole.CATALYST)),
        (ligated.id,), geometry=truseq_ligation_geometry(), source=source), (ligated,))

    forward = add_input("truseq_universal_primer", "TruSeq universal PCR primer",
                        MaterialKind.OLIGO)
    reverse = add_input("truseq_index_primer_pool", "selected TruSeq indexed PCR primer",
                        MaterialKind.OLIGO)
    library = Material("cut_run_library", "CUT&RUN sequencing library",
                       MaterialKind.LIBRARY, state=FINAL_LIBRARY, source=source)
    graph.add_step(Step(
        "library_pcr", "PCR-amplify the ligated library", StepKind.PCR,
        (Input(ligated.id, InputRole.TEMPLATE),
         Input(forward.id, InputRole.STRUCTURAL),
         Input(reverse.id, InputRole.STRUCTURAL)),
        (library.id,), geometry=PCRGeometry.regional(
            FINAL_LIBRARY,
            PCRPrimer("TruSeq universal PCR primer", il.TRUSEQ_P5_FULL),
            reverse_pool="selected TruSeq indexed PCR primer",
            reverse_region="Index 1 / Read 2 arm",
            statement=("The cited procedure identifies TruSeq library amplification; "
                       "the selected eight-base index member is not reported.")),
        source=source), (library,))
    graph.mark_terminal(Terminal(library.id, Disposition.SEQUENCING), SequencingPlan(
        FINAL_LIBRARY, tuple(SEQ_PRIMERS), tuple(p.role for p in SEQ_PRIMERS),
        read_lengths={"Read 1": 25, "Read 2": 25}))

    graph.declare_reviewed_scope((
        "nuclei", "cona_beads", "primary_antibody", "pa_mnase",
        "target_chromatin", "calcium", "targeted_cleavage", "proteinase_k",
        "recover_dna", "end_repair_enzymes", "end_repair", "a_tailing_enzyme",
        "da_tail", "truseq_adapters", "dna_ligase", "adapter_ligation",
        "truseq_universal_primer", "truseq_index_primer_pool", "library_pcr",
        "cut_run_library",
    ), source)
    return graph.require_publishable()
