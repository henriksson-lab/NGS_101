"""Protocol-owned source, oligos, workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import ThreePrimeOverhangAdapter, ligation_product_scene, react_dA_tail_then_ligate, repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal

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
# Dataset S1 Table S1: 36 of 44 libraries 10^5-fold/18 cycles, six 10^7-fold/20 cycles,
# two low-input PMS2-/- libraries 10-fold/20 cycles.
PCR_CYCLES = (18, 20)

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

TRUSEQ_DA_LIGATION = AdapterLigationGeometry.unavailable(
    "The original specifies TruSeq DNA PCR-Free Y adapters, but the study-specific2016 "
    "kit/index-set revision and physical adapter ends remain unassigned. Canonical A/T "
    "compatibility diagrams are inferred hypotheses, not a source-exact reagent reconstruction.")


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
    families = molecule("amplified bottlenecked library", tuple(library))
    w = Workflow(MolecularState("Sheared genomic DNA", tuple(dna.rows())))
    w.react("End-repair the sheared fragments",
            Scene.duplex(list(repaired), unpaired=("genomic insert'",)).rows(),
            name="end-repaired DNA")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each end",
        ligation_action="INFERRED — canonical TruSeq adapter compatibility",
        ligation_note="Adapter orientation distinguishes Watson and Crick descendants.")
    w.react("Dilute to about one haploid-genome equivalent", ligated_plain.rows(),
            name="bottlenecked adapter-ligated library",
            note="For 500 ng input the typical dilution was 100,000-fold; it is imposed before PCR.")
    w.react("PCR-amplify and paired-end sequence redundantly", families.rows(),
            name="amplified BotSeqS DNA (families are computational groups)",
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


# Local canonical compatibility model for the source-named TruSeq PCR-Free kit.
# Source proves the terminal22-base TS-PCR primers; the rest of the adapter and
# single-index arrangement is a hypothesis until the study kit revision is assigned.
# No synthetic UMI is installed. Native endpoint coordinates are an analysis key.
def final_library():
    con = Construct([
        Segment("P5 terminal TS-PCR1 region", TS_PCR_OLIGO1, "p5"),
        Segment("remaining P5/Read1 arm", il.TRUSEQ_P5_FULL[len(TS_PCR_OLIGO1):], "r1", inferred=True),
        role("genomic insert"), Segment("dA junction", "A", inferred=True),
        Segment("Read2 adapter arm", il.INDEX1_PRIMER, "r2", inferred=True),
        role("i7 index: member/length unassigned", molecular_feature=feature("sample_index_i7", "sample_index", "unknown"), inferred=True),
        Segment("P7 internal region", il.P7_RC[:-len(TS_PCR_OLIGO2)], "p7", inferred=True),
        Segment("terminal TS-PCR2 complement", revcomp(TS_PCR_OLIGO2), "p7"),
    ], name="BotSeqS canonical TruSeq compatibility model (inferred adapter arrangement)")
    primers=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    problems=sp.verify(con,primers,required_roles=tuple(p.role for p in primers))
    if problems:raise ValueError("BotSeqS candidate readout: "+"; ".join(problems))
    return con,primers

FINAL_LIBRARY, SEQ_PRIMERS = final_library()
TITLE = 'BotSeqS — bottleneck sequencing system'
NOTES = '01_botseqs.html'
SOURCE = '<a href="https://doi.org/10.1073/pnas.1607794113">Hoang et al., <i>PNAS</i> (2016)</a>.'
SUMMARY = 'A TruSeq PCR-free library is diluted after adapter ligation and before PCR; native fragment endpoints group redundant Watson and Crick read families without a synthetic UMI.'
CAVEAT = ('INFERRED — the canonical single-index adapter arrangement is a compatibility model; the source-named2016 kit/index-set revision and physical ends remain unassigned. The post-dilution PCR primers (TS-PCR Oligo1/Oligo2, one 3-prime phosphorothioate), '
          'Q5 Hot Start polymerase, cycling, 1.0x AMPure cleanup and the pBlueScript-carrier serial '
          'dilution are from the PNAS SI Appendix; per-library dilution (10-, 10^5- or 10^7-fold) '
          'and cycle number (18 or 20) are from Dataset S1 Table S1. The TruSeq PCR-Free adapter '
          'index set and whether an index read was acquired are not stated. '
          'INFERRED: the 350-bp size selection is placed after end repair following the 2017 '
          'TruSeq DNA PCR-Free guide; the 2016 kit revision used was not inspected.')
READ_LENGTHS = {'Read 1': 100, 'Read 2': 100}
FINAL_CAPTION = 'INFERRED — canonical TruSeq single-index compatibility model; error correction comes from the post-ligation bottleneck and shared genomic endpoints, not an added barcode.'
SEQUENCING_INTRO = 'Paired reads establish both endpoints and P5/P7 orientation for Watson- and Crick-family consensus.'
def oligos():
    return [f"TS-PCR Oligo1 (SI Appendix): 5′-{TS_PCR_OLIGO1_WRITTEN}-3′",
            f"TS-PCR Oligo2 (SI Appendix): 5′-{TS_PCR_OLIGO2_WRITTEN}-3′"]
def final_rows():
    return _botseqs_final_rows(FINAL_LIBRARY)
def workflow():
    return _botseqs_workflow(FINAL_LIBRARY)

def adapter_ligated_precursor():
    """Actual pre-PCR population without assuming the candidate final duplex."""
    sc = Scene()
    sc.strand("Y-adapter precursor population",
        [role("Y-adapter-ligated genomic molecule; native strands and ends unresolved", inferred=True)],
        label="source-named PCR-Free Y-adapter population; physical architecture unassigned")
    return sc


def protocol_graph():
    library = FINAL_LIBRARY
    precursor = adapter_ligated_precursor()
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
                           "TruSeq DNA PCR-Free Y adapters (kit/index revision unassigned)",
                           MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "TruSeq DNA PCR-Free ligation reagents",
                         MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "adapter-ligated pool (canonical architecture inferred)",
                            MaterialKind.MOLECULAR, state=precursor, inferred=True, source=source+" source-named kit; exact2016 adapter architecture not assigned")
    graph.react("adapter_ligation", "INFERRED architecture — ligate source-named kit Y adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=TRUSEQ_DA_LIGATION)

    split_condition=graph.input("reference_split_condition",
        "SI WGS route:10 µL of the20 µL adapter pool for matched WGS when prepared from this pool; "
        "the remaining source material supplies the bottleneck", MaterialKind.CONDITION, source=SI_PCR)
    bot_aliquot=graph.product("botseqs_aliquot","adapter-ligated material assigned to BotSeqS",MaterialKind.MOLECULAR,state=adapted.state)
    reference=graph.product("reference_aliquot","optional matched-WGS aliquot (source subset)",MaterialKind.MOLECULAR,state=adapted.state)
    graph.react("reference_split","Retain an optional matched-WGS aliquot before bottleneck dilution",StepKind.SPLIT,
        ((adapted,InputRole.TEMPLATE),(split_condition,InputRole.CONDITION)),(bot_aliquot,reference),changes_molecule=False,
        source="PNAS SI Whole-genome sequencing; main article describes same DNA or same library references")
    wgs_forward=graph.input("wgs_forward","forward component of TruSeq PCR Primer Cocktail (member unassigned)",MaterialKind.OLIGO)
    wgs_reverse=graph.input("wgs_reverse","reverse component of TruSeq PCR Primer Cocktail (member unassigned)",MaterialKind.OLIGO)
    wgs_mix=graph.input("wgs_polymerase","TruSeq PCR Master Mix per manufacturer instructions",MaterialKind.CATALYST)
    wgs_amp=graph.product("wgs_amplified","matched-WGS PCR product",MaterialKind.MOLECULAR,state=library)
    graph.react("wgs_pcr","Amplify the retained matched-WGS aliquot",StepKind.PCR,
        ((reference,InputRole.TEMPLATE),(wgs_forward,InputRole.STRUCTURAL),(wgs_reverse,InputRole.STRUCTURAL),
         (wgs_mix,InputRole.CATALYST)),(wgs_amp,),geometry=PCRGeometry.unavailable(reference.state,
            reason="SI names the TruSeq PCR Primer Cocktail/Master Mix; selected vendor primer members "
                   "and native study-specific adapter geometry are unassigned."))
    wgs_beads=graph.input("wgs_beads","1.0× AMPure XP for WGS PCR cleanup",MaterialKind.SELECTION)
    wgs_library=graph.product("wgs_library","matched-WGS reference library (where prepared)",MaterialKind.LIBRARY,state=library)
    wgs_waste=graph.product("wgs_waste","WGS cleanup supernatant",MaterialKind.WASTE)
    graph.react("wgs_cleanup","Purify the matched-WGS PCR with1.0× AMPure",StepKind.CLEANUP,
        ((wgs_amp,InputRole.TEMPLATE),(wgs_beads,InputRole.SELECTION)),(wgs_library,wgs_waste),changes_molecule=False)
    graph.discard(wgs_waste,"physically removed AMPure supernatant")
    graph.sequence(wgs_library,SequencingPlan(acquired_roles=("Read 1","Read 2"),
        read_lengths={"Read 1":100,"Read 2":100},limitation="SI WGS route: HiSeq2500 PE2×100,>30×. "
        "The matched reference was available for32/34 individuals; two used Sanger. "
        "Exact primer/index geometry is unassigned; no index read is guessed."))

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
                               MaterialKind.MOLECULAR, state=bot_aliquot.state)
    excluded = graph.product("excluded_molecules", "unsampled adapter stock and serial-dilution material",
                             MaterialKind.MOLECULAR)
    graph.react("bottleneck_sampling", "Dilute to about one haploid-genome equivalent",
                StepKind.SAMPLING,
                ((bot_aliquot, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)),
                (bottleneck, excluded), changes_molecule=False)
    graph.mark_terminal(Terminal(excluded.id,Disposition.UNRESOLVED,
        "These molecules do not enter the sampled PCR aliquot. Remaining source/serial stock "
        "is not thereby a physical discard; matched WGS retention is a separate source branch."))

    forward = graph.input("pcr_forward",
                          f"TS-PCR Oligo1 5′-{TS_PCR_OLIGO1_WRITTEN}-3′, HPLC-purified, 0.5 µM",
                          MaterialKind.OLIGO, state=molecule("TS-PCR Oligo1",(Segment("primer",TS_PCR_OLIGO1,note="one terminal3′ phosphorothioate"),)), source=SI_PCR)
    reverse = graph.input("pcr_reverse",
                          f"TS-PCR Oligo2 5′-{TS_PCR_OLIGO2_WRITTEN}-3′, HPLC-purified, 0.5 µM",
                          MaterialKind.OLIGO, state=molecule("TS-PCR Oligo2",(Segment("primer",TS_PCR_OLIGO2,note="one terminal3′ phosphorothioate"),)), source=SI_PCR)
    polymerase = graph.input("pcr_polymerase",
                             "Q5 Hot Start High-Fidelity 2x Master Mix (NEB), 50 µL; 98 °C 30 s, "
                             "18 or 20 x (98 °C 10 s, 69 °C 30 s, 72 °C 30 s), 72 °C 2 min",
                             MaterialKind.CATALYST, source=SI_PCR)
    amplified = graph.product("amplified_library", "redundantly amplified BotSeqS PCR product",
                              MaterialKind.MOLECULAR, state=library)
    graph.react("library_pcr", "PCR-amplify the bottlenecked library", StepKind.PCR,
                ((bottleneck, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), source=SI_PCR, geometry=PCRGeometry.unavailable(
                    bottleneck.state,
                    reason="The incoming PCR-free Y-adapter population is the PCR substrate. "
                    "SI specifies both modified primers and program; exact native2016 adapter "
                    "ends and first-cycle geometry are unassigned. The inferred canonical final "
                    "library is a compatibility model, not a substitute duplex template."))
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
        "adapter_ligation", "reference_split_condition", "reference_split", "wgs_forward", "wgs_reverse", "wgs_polymerase", "wgs_pcr", "wgs_beads", "wgs_cleanup", "wgs_library", "dilution_bottleneck", "bottleneck_sampling",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr", "amplified_library",
        "ampure", "pcr_cleanup", "pcr_cleanup_waste", "sequencing_library",
    ), source)
    return graph.require_publishable()
