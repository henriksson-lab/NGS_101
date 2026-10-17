"""Protocol-owned source, oligos, workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import UNASSIGNED_LENGTH, nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import ThreePrimeOverhangAdapter, ligation_product_scene, react_dA_tail_then_ligate, repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from accuracy_protocols import _workflow

BOTSEQS_TYPICAL_DILUTION_FOLD = 100_000

BOTSEQS_INSERT_BP = 350

BOTSEQS_READ_BP = 100

# PNAS SI Appendix (Hoang et al. 2016), SI Materials and Methods, "Dilution of
# Y-adapter-ligated molecules" and "PCR amplification of diluted Y-adapter-ligated
# molecules" (sapp.pdf.txt lines 4-28). The same primers and program appear in US
# 12,359,250 B2 Example 1, read earlier as a lead. '*' = one 3' phosphorothioate bond.
SI_PCR = ("Hoang et al. PNAS (2016) SI Appendix, SI Materials and Methods: 'PCR amplification "
          "of diluted Y-adapter-ligated molecules'; cycles per library from Dataset S1 Table S1")
TS_PCR_OLIGO1_WRITTEN = "AATGATACGGCGACCACCGAG*A"
TS_PCR_OLIGO2_WRITTEN = "CAAGCAGAAGACGGCATACGA*G"
TS_PCR_OLIGO1 = TS_PCR_OLIGO1_WRITTEN.replace("*", "")
TS_PCR_OLIGO2 = TS_PCR_OLIGO2_WRITTEN.replace("*", "")
# Dataset S1 Table S1: 35 of 43 libraries 10^5-fold/18 cycles, six 10^7-fold/20 cycles,
# two low-input PMS2-/- libraries 10-fold/20 cycles.
PCR_CYCLES = (18, 20)

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

TRUSEQ_DA_LIGATION = reported_adapter_end(
    "TruSeq Y adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    statement="The documented TruSeq Y adapter presents a single 3-prime dT to a dA-tailed insert.")


def _botseqs_workflow(library: Construct) -> Workflow:
    dna = molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),))
    repaired = Construct([role("genomic insert")], name="end-repaired genomic insert")
    prepared = repair_and_dA_tail(repaired)
    ligated = ligation_product_scene(library, "genomic insert",
                                     label="TruSeq PCR-free molecules",
                                     adapter_group="Illumina Y adapter")
    ligated_plain = ligation_product_scene(
        library, "genomic insert", label="bottlenecked adapter-ligated library",
        adapter_group="Illumina Y adapter", highlight_components=False)
    families = molecule("redundantly sequenced bottleneck", (
        role("Watson read family"), role("same two genomic endpoints"), role("Crick read family")))
    w = Workflow(MolecularState("Sheared genomic DNA", tuple(dna.rows())))
    w.react("End-repair the sheared fragments",
            Scene.duplex(list(repaired), unpaired=("genomic insert'",)).rows(),
            name="end-repaired DNA")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each end",
        ligation_action="Ligate Illumina Y adapters",
        ligation_note="Adapter orientation distinguishes Watson and Crick descendants.")
    w.react("Dilute to about one haploid-genome equivalent", ligated_plain.rows(),
            name="bottlenecked adapter-ligated library",
            note="For 500 ng input the typical dilution was 100,000-fold; it is imposed before PCR.")
    w.react("PCR-amplify and paired-end sequence redundantly", families.rows(),
            name="Watson and Crick read families",
            note="The two genomic endpoints are the endogenous family key; no synthetic UMI is used. "
                 "P5/P7-terminal TS-PCR primers, Q5 Hot Start, 18 or 20 cycles and a 1.0x AMPure "
                 "cleanup (SI Appendix).")
    return w


def _botseqs_final_rows(library: Construct):
    """Draw the actual adapter-ligation boundaries on the canonical TruSeq molecule."""
    sc = ligation_product_scene(library, "genomic insert", label=library.name,
                                adapter_group="Illumina Y adapter",
                                highlight_components=False)
    sc.labels("top")
    return sc.rows()


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([role("genomic insert")], "BotSeqS TruSeq PCR-free library", dual_index=False,
    # The TruSeq PCR-Free index set (and so its index length) is not stated.
    i7_nt=UNASSIGNED_LENGTH)
TITLE = 'BotSeqS — bottleneck sequencing system'
NOTES = '01_botseqs.html'
SOURCE = '<a href="https://doi.org/10.1073/pnas.1607794113">Hoang et al., <i>PNAS</i> (2016)</a>.'
SUMMARY = 'A TruSeq PCR-free library is diluted after adapter ligation and before PCR; native fragment endpoints group redundant Watson and Crick read families without a synthetic UMI.'
CAVEAT = ('The post-dilution PCR primers (TS-PCR Oligo1/Oligo2, one 3-prime phosphorothioate), '
          'Q5 Hot Start polymerase, cycling, 1.0x AMPure cleanup and the pBlueScript-carrier serial '
          'dilution are from the PNAS SI Appendix; per-library dilution (10-, 10^5- or 10^7-fold) '
          'and cycle number (18 or 20) are from Dataset S1 Table S1. The TruSeq PCR-Free adapter '
          'index set and whether an index read was acquired are not stated. '
          'INFERRED: the 350-bp size selection is placed after end repair following the 2017 '
          'TruSeq DNA PCR-Free guide; the 2016 kit revision used was not inspected.')
READ_LENGTHS = {'Read 1': 100, 'Read 2': 100}
FINAL_CAPTION = 'Ordinary TruSeq paired-end molecule; error correction comes from the post-ligation bottleneck and shared genomic endpoints, not an added barcode.'
SEQUENCING_INTRO = 'Paired reads establish both endpoints and P5/P7 orientation for Watson- and Crick-family consensus.'
def oligos():
    return [f"TS-PCR Oligo1 (SI Appendix): 5′-{TS_PCR_OLIGO1_WRITTEN}-3′",
            f"TS-PCR Oligo2 (SI Appendix): 5′-{TS_PCR_OLIGO2_WRITTEN}-3′"]
def final_rows():
    return _botseqs_final_rows(FINAL_LIBRARY)
def workflow():
    return _botseqs_workflow(FINAL_LIBRARY)

def protocol_graph():
    library = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = "Hoang et al., PNAS (2016); doi:10.1073/pnas.1607794113"
    graph = ProtocolGraph("BotSeqS", authority=source)
    genomic = graph.input(
        "genomic_dna", "Qiagen AllPrep tissue genomic DNA, 34 ng–1 µg in 55 µL TE",
        MaterialKind.BIOLOGICAL,
        state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    shearing = graph.input("shearing_condition",
                           "BioRuptor, high intensity, 15 s on / 90 s off, seven cycles at 3 °C",
                           MaterialKind.CONDITION)
    fragments = graph.product(
        "sheared_fragments", "sheared genomic DNA with native endpoints",
        MaterialKind.MOLECULAR,
        state=molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),)))
    graph.react("shear", "Shear genomic DNA", StepKind.FRAGMENTATION,
                ((genomic, InputRole.TEMPLATE), (shearing, InputRole.CONDITION)),
                (fragments,))
    repair = graph.input("end_repair_mix",
                         "TruSeq DNA PCR-Free end-repair reagents (standard low-input protocol)",
                         MaterialKind.CATALYST)
    repaired_state = Construct([role("genomic insert")], name="end-repaired genomic insert")
    repaired = graph.product("repaired_fragments", "end-repaired genomic fragments",
                             MaterialKind.MOLECULAR, state=repaired_state)
    graph.react("end_repair", "End-repair the sheared fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    size_beads = graph.input(
        "size_selection", "TruSeq DNA PCR-Free bead selection for 350-bp inserts",
        MaterialKind.SELECTION, inferred=True,
        source=source + "; Methods: 'selection for 350-bp insert sizes'. Placement after end "
               "repair follows the October 2017 TruSeq DNA PCR-Free guide (1000000039279 v00), "
               "not the uninspected 2016 kit revision")
    selected = graph.product("size_selected", "350-bp-insert-selected end-repaired fragments",
                             MaterialKind.MOLECULAR, state=repaired_state, inferred=True)
    off_size = graph.product("off_size_fragments", "fragments outside the selected size window",
                             MaterialKind.WASTE, inferred=True)
    graph.react("select_size", "Select the 350-bp insert size", StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (size_beads, InputRole.SELECTION)),
                (selected, off_size), changes_molecule=False)
    graph.discard(off_size, "bead supernatants discarded during size selection")
    tailing = graph.input("a_tailing_enzyme", "TruSeq DNA PCR-Free A-tailing mix",
                          MaterialKind.CATALYST)
    tailed_state = repair_and_dA_tail(repaired_state)
    tailed = graph.product("tailed_fragments", "dA-tailed genomic fragments",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react("a_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((selected, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,))
    adapters = graph.input("truseq_y_adapters",
                           "TruSeq DNA PCR-Free indexed Y adapters (index set not stated)",
                           MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "TruSeq DNA PCR-Free ligation reagents",
                         MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "TruSeq adapter-ligated genomic library pool",
                            MaterialKind.MOLECULAR, state=library)
    graph.react("adapter_ligation", "Ligate Illumina Y adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=TRUSEQ_DA_LIGATION)

    dilution = graph.input(
        "dilution_bottleneck",
        "dilution of Y-adapter-ligated DNA (20 µL): serial 2 µL into 18 µL transfers in TE with "
        "1 ng/µL pBlueScript carrier, 2 µL used as PCR template; typically 10^5-fold for 500 ng "
        "input, chosen from 10^3–10^7 dilutions test-sequenced 2x50 on MiSeq (SI Fig. S1); "
        "Dataset S1 lists 10-, 10^5- and 10^7-fold final dilutions",
        MaterialKind.CONDITION,
        source="Hoang et al. PNAS (2016) SI Appendix 'Dilution of Y-adapter-ligated molecules', "
               "'MiSeq run and analysis', Fig. S1; Dataset S1 Table S1")
    bottleneck = graph.product("bottleneck", "bottlenecked adapter-ligated molecules",
                               MaterialKind.MOLECULAR, state=library)
    excluded = graph.product("excluded_molecules", "adapter-ligated molecules excluded by dilution",
                             MaterialKind.WASTE)
    graph.react("bottleneck_sampling", "Dilute to about one haploid-genome equivalent",
                StepKind.SAMPLING,
                ((adapted, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)),
                (bottleneck, excluded), changes_molecule=False)
    graph.discard(excluded, "falls outside the deliberately sampled molecular bottleneck")

    forward = graph.input("pcr_forward",
                          f"TS-PCR Oligo1 5′-{TS_PCR_OLIGO1_WRITTEN}-3′, HPLC-purified, 0.5 µM",
                          MaterialKind.OLIGO, state=TS_PCR_OLIGO1, source=SI_PCR)
    reverse = graph.input("pcr_reverse",
                          f"TS-PCR Oligo2 5′-{TS_PCR_OLIGO2_WRITTEN}-3′, HPLC-purified, 0.5 µM",
                          MaterialKind.OLIGO, state=TS_PCR_OLIGO2, source=SI_PCR)
    polymerase = graph.input("pcr_polymerase",
                             "Q5 Hot Start High-Fidelity 2x Master Mix (NEB), 50 µL; 98 °C 30 s, "
                             "18 or 20 x (98 °C 10 s, 69 °C 30 s, 72 °C 30 s), 72 °C 2 min",
                             MaterialKind.CATALYST, source=SI_PCR)
    amplified = graph.product("amplified_library", "redundantly amplified BotSeqS PCR product",
                              MaterialKind.MOLECULAR, state=library)
    graph.react("library_pcr", "PCR-amplify the bottlenecked library", StepKind.PCR,
                ((bottleneck, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), source=SI_PCR, geometry=PCRGeometry.exact(
                    library, PCRPrimer("TS-PCR Oligo1", TS_PCR_OLIGO1),
                    PCRPrimer("TS-PCR Oligo2", TS_PCR_OLIGO2),
                    statement="Primer sequences printed in the PNAS SI Appendix; sites computed on "
                              "the canonical TruSeq P5/P7 termini. The sample index is carried by the "
                              "ligated PCR-free adapter, not added by this PCR."))
    beads = graph.input("ampure", "AMPure XP beads, 1.0x bead-to-sample ratio",
                        MaterialKind.CONDITION, source=SI_PCR)
    final = graph.product("sequencing_library", "purified BotSeqS library",
                          MaterialKind.LIBRARY, state=library)
    pcr_waste = graph.product("pcr_cleanup_waste", "1.0x AMPure supernatant (primers, short products)",
                              MaterialKind.WASTE)
    graph.react("pcr_cleanup", "Purify the PCR product with 1.0x AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
                (final, pcr_waste), source=SI_PCR)
    graph.discard(pcr_waste, "AMPure supernatant removed")
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), ("Read 1", "Read 2"),
        {"Read 1": BOTSEQS_READ_BP, "Read 2": BOTSEQS_READ_BP}))
    graph.declare_reviewed_scope((
        "genomic_dna", "shearing_condition", "shear", "end_repair_mix", "end_repair",
        "size_selection", "select_size", "a_tailing_enzyme", "a_tail", "truseq_y_adapters", "adapter_ligase",
        "adapter_ligation", "dilution_bottleneck", "bottleneck_sampling",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "amplified_library",
        "ampure", "pcr_cleanup", "pcr_cleanup_waste", "sequencing_library",
    ), source)
    return graph.require_publishable()
