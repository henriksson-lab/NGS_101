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

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

TRUSEQ_DA_LIGATION = reported_adapter_end(
    "TruSeq Y adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    statement="The documented TruSeq Y adapter presents a single 3-prime dT to a dA-tailed insert.")


def _cirseq_workflow(library: Construct) -> Workflow:
    rna = molecule("purified viral RNA", (role("viral RNA"),))
    fragment = molecule("size-selected RNA", (role("85–100-nt RNA fragment", length_bp=92),))
    circ = closed_circle("circular RNA template", (role("RNA fragment"),),
                         closure_label="RNA intramolecular ligation")
    tandem = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    ds_insert = Construct([role("tandem-repeat cDNA")],
                          name="double-stranded tandem-repeat cDNA")
    ds = Scene.duplex(list(ds_insert), label=ds_insert.name,
                      unpaired=("tandem-repeat cDNA'",))
    prepared = repair_and_dA_tail(ds_insert)
    ligated = ligation_product_scene(library, "tandem-repeat cDNA",
                                     label="TruSeq-ligated CirSeq library",
                                     adapter_group="TruSeq indexed adapter")
    w = Workflow(MolecularState("Starting material", tuple(rna.rows())))
    w.react("Fragment with Zn2+ and size-select", fragment.rows(), name="Fragment and size-select",
            note="Fragments are kept short enough for at least three copies in one read.")
    w.react("5′-phosphorylate and circularize RNA", circ.rows(), name="Circularize RNA",
            note="T4 RNA ligase 1 closes the native 3′→5′ junction.")
    w.react("Random-prime rolling-circle reverse transcription", tandem.rows(),
            name="Rolling-circle reverse transcription",
            note="Physically linked repeats derive independently from the same RNA template.")
    w.react("Synthesize the second strand and blunt the cDNA", ds.rows(),
            name="double-stranded tandem-repeat cDNA")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each cDNA end",
        ligation_action="Ligate TruSeq indexed adapters",
        ligation_note="Size selection retains tandem inserts long enough for within-read consensus.")
    return w


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([role("tandem-repeat cDNA")], "CirSeq library", dual_index=False, i7_nt=6)
TITLE = 'CirSeq — circular resequencing of RNA virus populations'
NOTES = '01_cirseq.html'
SOURCE = '<a href="https://doi.org/10.1038/nprot.2014.118">Acevedo and Andino, <i>Nature Protocols</i> (2014)</a>.'
SUMMARY = 'Short viral RNA fragments are circularized and copied repeatedly by rolling-circle reverse transcription, placing at least three linked template copies inside one Illumina read.'
CAVEAT = 'The named FC-121-4001 kit contains twelve six-base Set A indexes, but the protocol does not identify the index member used for a particular library. The defining readout is single-end; the opposite adapter retains a Read 2 primer site even though no paired read is reported.'
READ_LENGTHS = {'Read 1': 300}
FINAL_CAPTION = 'TruSeq adapters flank a double-stranded cDNA insert containing tandem copies of one circular RNA fragment.'
SEQUENCING_INTRO = "The protocol used single-end 300-cycle MiSeq reads. The adapter has an index and an opposite-end primer site, but the paper does not specify an index read or a paired Read 2 run."
def oligos():
    return ()
def final_rows():
    return ligation_product_scene(FINAL_LIBRARY, "tandem-repeat cDNA", label=FINAL_LIBRARY.name, adapter_group="TruSeq indexed adapter", highlight_components=False).rows()
def workflow():
    return _cirseq_workflow(FINAL_LIBRARY)

def protocol_graph():
    library = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    source = "Acevedo and Andino, Nature Protocols (2014); doi:10.1038/nprot.2014.118"
    graph = ProtocolGraph("CirSeq", authority=source)
    viral_rna = graph.input(
        "viral_rna", "purified viral RNA", MaterialKind.BIOLOGICAL,
        state=molecule("purified viral RNA", (role("viral RNA"),)), primary=True)
    zinc = graph.input("zinc_fragmentation", "Zn2+ fragmentation condition",
                       MaterialKind.CONDITION)
    fragments = graph.product(
        "rna_fragments", "fragmented viral RNA", MaterialKind.MOLECULAR,
        state=molecule("fragmented viral RNA", (role("short viral RNA fragments"),)))
    graph.react("fragment", "Fragment viral RNA with Zn2+", StepKind.FRAGMENTATION,
                ((viral_rna, InputRole.TEMPLATE), (zinc, InputRole.CONDITION)),
                (fragments,))
    size_rule = graph.input("size_selection", "85–100-nt size-selection rule",
                            MaterialKind.SELECTION)
    selected = graph.product(
        "selected_rna", "85–100-nt viral RNA fragments", MaterialKind.MOLECULAR,
        state=molecule("size-selected RNA", (role("85–100-nt RNA fragment", length_bp=92),)))
    rejected = graph.product("rejected_sizes", "RNA outside the 85–100-nt window",
                             MaterialKind.WASTE)
    graph.react("select_size", "Select 85–100-nt RNA fragments", StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not fit at least three linked copies within the long read")

    kinase = graph.input("polynucleotide_kinase", "5-prime polynucleotide kinase",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated_rna", "5-prime-phosphorylated RNA fragments",
                                   MaterialKind.MOLECULAR, state=selected.state)
    graph.react("phosphorylate", "5-prime-phosphorylate selected RNA",
                StepKind.PHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 1", MaterialKind.CATALYST)
    circle_state = closed_circle("circular RNA template", (role("RNA fragment"),),
                                 closure_label="RNA intramolecular ligation")
    circles = graph.product("rna_circles", "circularized viral RNA fragments",
                            MaterialKind.MOLECULAR, state=circle_state)
    graph.react("circularize", "Circularize RNA fragments", StepKind.CIRCULARIZATION,
                ((phosphorylated, InputRole.TEMPLATE), (rna_ligase, InputRole.CATALYST)),
                (circles,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes phosphorylated RNA fragments with T4 RNA ligase 1."))

    random_rt = graph.input("random_rt_primer", "random reverse-transcription primer",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "rolling-circle reverse transcriptase",
                                       MaterialKind.CATALYST)
    tandem_state = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    tandem = graph.product("tandem_cdna", "rolling-circle tandem-repeat cDNA",
                           MaterialKind.MOLECULAR, state=tandem_state)
    graph.react("rolling_circle_rt", "Random-prime rolling-circle reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((circles, InputRole.TEMPLATE), (random_rt, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (tandem,))
    second_primer = graph.input("second_strand_primer", "random second-strand primer",
                                MaterialKind.OLIGO)
    dna_polymerase = graph.input("second_strand_polymerase", "second-strand DNA polymerase",
                                 MaterialKind.CATALYST)
    ds_construct = Construct([role("tandem-repeat cDNA")],
                             name="double-stranded tandem-repeat cDNA")
    double_stranded = graph.product("double_stranded_cdna", ds_construct.name,
                                    MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((tandem, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (dna_polymerase, InputRole.CATALYST)), (double_stranded,))
    repair = graph.input("end_repair_mix", "blunt-end repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired_cdna", "blunt tandem-repeat cDNA",
                             MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("end_repair", "Blunt the double-stranded cDNA", StepKind.END_REPAIR,
                ((double_stranded, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed_cdna", "dA-tailed tandem-repeat cDNA",
                           MaterialKind.MOLECULAR,
                           state=repair_and_dA_tail(ds_construct))
    graph.react("a_tail", "Add one 3-prime dA to each cDNA end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("truseq_adapters", "TruSeq indexed adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "DNA adapter ligase", MaterialKind.CATALYST)
    final = graph.product("sequencing_library", "TruSeq-ligated CirSeq library",
                          MaterialKind.LIBRARY, state=library)
    graph.react("adapter_ligation", "Ligate TruSeq indexed adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (final,), geometry=TRUSEQ_DA_LIGATION)
    graph.sequence(final, SequencingPlan(
        library, (primers[0],), ("Read 1",),
        READ_LENGTHS, {"Read 1": "include"}))
    graph.declare_reviewed_scope((
        "viral_rna", "zinc_fragmentation", "fragment", "size_selection", "select_size",
        "polynucleotide_kinase", "phosphorylate", "rna_ligase", "circularize",
        "random_rt_primer", "reverse_transcriptase", "rolling_circle_rt",
        "second_strand_primer", "second_strand_polymerase", "second_strand",
        "end_repair_mix", "end_repair", "a_tailing_enzyme", "a_tail",
        "truseq_adapters", "adapter_ligase", "adapter_ligation", "sequencing_library",
    ), source)
    return graph.require_publishable()
