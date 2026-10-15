"""Protocol-owned source, oligos, workflow and reaction graph."""
import illumina as il
import nextera as nx
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg, truseq_library
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
            note="The two genomic endpoints are the endogenous family key; no synthetic UMI is used.")
    return w


def _botseqs_final_rows(library: Construct):
    """Draw the actual adapter-ligation boundaries on the canonical TruSeq molecule."""
    sc = ligation_product_scene(library, "genomic insert", label=library.name,
                                adapter_group="Illumina Y adapter",
                                highlight_components=False)
    sc.labels("top")
    return sc.rows()


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([role("genomic insert")], "BotSeqS TruSeq PCR-free library", dual_index=False)
TITLE = 'BotSeqS — bottleneck sequencing system'
NOTES = '01_botseqs.html'
SOURCE = '<a href="https://doi.org/10.1073/pnas.1607794113">Hoang et al., <i>PNAS</i> (2016)</a>.'
SUMMARY = 'A TruSeq PCR-free library is diluted after adapter ligation and before PCR; native fragment endpoints group redundant Watson and Crick read families without a synthetic UMI.'
CAVEAT = ''
READ_LENGTHS = {'Read 1': 100, 'Read 2': 100}
FINAL_CAPTION = 'Ordinary TruSeq paired-end molecule; error correction comes from the post-ligation bottleneck and shared genomic endpoints, not an added barcode.'
SEQUENCING_INTRO = 'Paired reads establish both endpoints and P5/P7 orientation for Watson- and Crick-family consensus.'
def oligos():
    return ()
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
        "genomic_dna", "genomic DNA", MaterialKind.BIOLOGICAL,
        state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True)
    shearing = graph.input("shearing_condition", "mechanical DNA-shearing condition",
                           MaterialKind.CONDITION)
    fragments = graph.product(
        "sheared_fragments", "sheared genomic DNA with native endpoints",
        MaterialKind.MOLECULAR,
        state=molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),)))
    graph.react("shear", "Shear genomic DNA", StepKind.FRAGMENTATION,
                ((genomic, InputRole.TEMPLATE), (shearing, InputRole.CONDITION)),
                (fragments,))
    repair = graph.input("end_repair_mix", "end-repair enzyme mix", MaterialKind.CATALYST)
    repaired_state = Construct([role("genomic insert")], name="end-repaired genomic insert")
    repaired = graph.product("repaired_fragments", "end-repaired genomic fragments",
                             MaterialKind.MOLECULAR, state=repaired_state)
    graph.react("end_repair", "End-repair the sheared fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed_state = repair_and_dA_tail(repaired_state)
    tailed = graph.product("tailed_fragments", "dA-tailed genomic fragments",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react("a_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,))
    adapters = graph.input("truseq_y_adapters", "TruSeq Y adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "DNA adapter ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "TruSeq adapter-ligated genomic library pool",
                            MaterialKind.MOLECULAR, state=library)
    graph.react("adapter_ligation", "Ligate Illumina Y adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=TRUSEQ_DA_LIGATION)

    dilution = graph.input("dilution_bottleneck", "approximately one-genome-equivalent dilution",
                           MaterialKind.CONDITION)
    bottleneck = graph.product("bottleneck", "bottlenecked adapter-ligated molecules",
                               MaterialKind.MOLECULAR, state=library)
    excluded = graph.product("excluded_molecules", "adapter-ligated molecules excluded by dilution",
                             MaterialKind.WASTE)
    graph.react("bottleneck_sampling", "Dilute to about one haploid-genome equivalent",
                StepKind.SAMPLING,
                ((adapted, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)),
                (bottleneck, excluded), changes_molecule=False)
    graph.discard(excluded, "falls outside the deliberately sampled molecular bottleneck")

    forward = graph.input("pcr_forward", "TruSeq library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("pcr_reverse", "TruSeq indexed library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "redundantly amplified BotSeqS library",
                          MaterialKind.LIBRARY, state=library)
    graph.react("library_pcr", "PCR-amplify the bottlenecked library", StepKind.PCR,
                ((bottleneck, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (final,), geometry=PCRGeometry.unavailable(
                    library,
                    reason="The paper uses ordinary TruSeq amplification; the chosen index oligo varies by library."))
    graph.sequence(final, SequencingPlan(
        library, tuple(primers), tuple(p.role for p in primers),
        {"Read 1": BOTSEQS_READ_BP, "Read 2": BOTSEQS_READ_BP}))
    graph.declare_reviewed_scope((
        "genomic_dna", "shearing_condition", "shear", "end_repair_mix", "end_repair",
        "a_tailing_enzyme", "a_tail", "truseq_y_adapters", "adapter_ligase",
        "adapter_ligation", "dilution_bottleneck", "bottleneck_sampling",
        "pcr_forward", "pcr_reverse", "pcr_polymerase", "library_pcr",
        "sequencing_library",
    ), source)
    return graph.require_publishable()
