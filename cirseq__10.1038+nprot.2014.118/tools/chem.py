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
            note="Step 18: T4 polynucleotide kinase and T4 RNA ligase 1 act in one ATP-containing mix (37 °C, 30 min); the ligase closes the native 3′→5′ junction.")
    w.react("Random-prime rolling-circle reverse transcription", tandem.rows(),
            name="Rolling-circle reverse transcription",
            note="Random hexamers and SuperScript III; RNase H is added after 2 min at 42 °C. Physically linked repeats derive from the same RNA template.")
    w.react("Synthesize the second strand and blunt the cDNA", ds.rows(),
            name="double-stranded tandem-repeat cDNA",
            note="NEBNext second-strand synthesis module, then NEBNext end repair; the second-strand priming mechanism is not described by the protocol.")
    react_dA_tail_then_ligate(
        w, prepared, TRUSEQ_Y_END, ligated,
        dA_action="Add one 3-prime dA to each cDNA end",
        ligation_action="Ligate TruSeq indexed adapters",
        ligation_note="NEBNext Quick T4 DNA ligase with the TruSeq indexed adapter.")
    final = final_rows()
    w.react("Excise 450–600-nt adapter-ligated DNA from a 10% urea-PAGE gel", final,
            name="size-selected ligated library",
            note="Step 43: at least 435 bases = 135 adapter bases plus at least 300 bases of tandem-repeat insert.")
    w.react("Amplify 15 cycles with Phusion and the TruSeq PCR primer cocktail; purify on 7.5% PAGE", final,
            name="amplified CirSeq library",
            note="Steps 46–53. The primer-cocktail sequences are not printed in the protocol.")
    return w


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([role("tandem-repeat cDNA")], "CirSeq library", dual_index=False, i7_nt=6)
TITLE = 'CirSeq — circular resequencing of RNA virus populations'
NOTES = '01_cirseq.html'
SOURCE = '<a href="https://doi.org/10.1038/nprot.2014.118">Acevedo and Andino, <i>Nature Protocols</i> (2014)</a>.'
SUMMARY = 'Short viral RNA fragments are circularized and copied repeatedly by rolling-circle reverse transcription, placing at least three linked template copies inside one Illumina read.'
CAVEAT = 'The named FC-121-4001 kit contains twelve six-base Set A indexes, but the protocol does not identify the index member used for a particular library. The protocol recommends 300-nt reads (MiSeq v2 300-cycle kit, or HiSeq Rapid single-read kits); the opposite adapter retains a Read 2 primer site even though no paired read is described. The second-strand priming mechanism inside the NEBNext module and the PCR primer-cocktail sequences are not printed in the protocol, so those steps are drawn without primer geometry.'
READ_LENGTHS = {'Read 1': 300}
FINAL_CAPTION = 'TruSeq adapters flank a double-stranded cDNA insert containing tandem copies of one circular RNA fragment.'
SEQUENCING_INTRO = "The protocol recommends single 300-nt reads (MiSeq v2 300-cycle kit). The adapter has an index and an opposite-end primer site, but the paper does not specify an index read or a paired Read 2 run."
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

    kinase = graph.input("polynucleotide_kinase", "T4 polynucleotide kinase (Step 18 mix)",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product("phosphorylated_rna", "5-prime-phosphorylated RNA fragments",
                                   MaterialKind.MOLECULAR, state=selected.state)
    graph.react("phosphorylate", "5-prime-phosphorylate selected RNA (same Step 18 mix as ligation)",
                StepKind.PHOSPHORYLATION,
                ((selected, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 1 with 1 mM ATP", MaterialKind.CATALYST)
    circle_state = closed_circle("circular RNA template", (role("RNA fragment"),),
                                 closure_label="RNA intramolecular ligation")
    circles = graph.product("rna_circles", "circularized viral RNA fragments",
                            MaterialKind.MOLECULAR, state=circle_state)
    graph.react("circularize", "Circularize RNA fragments", StepKind.CIRCULARIZATION,
                ((phosphorylated, InputRole.TEMPLATE), (rna_ligase, InputRole.CATALYST)),
                (circles,), geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The protocol circularizes phosphorylated RNA fragments with T4 RNA ligase 1."))

    random_rt = graph.input("random_rt_primer", "random hexamers (Invitrogen N8080127)",
                            MaterialKind.OLIGO)
    reverse_transcriptase = graph.input("reverse_transcriptase", "SuperScript III reverse transcriptase",
                                       MaterialKind.CATALYST)
    rnase_h = graph.input("rnase_h", "RNase H (added after 2 min at 42 °C)", MaterialKind.CATALYST)
    tandem_state = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    tandem = graph.product("tandem_cdna", "rolling-circle tandem-repeat cDNA",
                           MaterialKind.MOLECULAR, state=tandem_state)
    graph.react("rolling_circle_rt", "Random-prime rolling-circle reverse transcription",
                StepKind.REVERSE_TRANSCRIPTION,
                ((circles, InputRole.TEMPLATE), (random_rt, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST),
                 (rnase_h, InputRole.CATALYST)), (tandem,))
    dna_polymerase = graph.input("second_strand_polymerase",
                                 "NEBNext second-strand synthesis enzyme mix (E6111)",
                                 MaterialKind.CATALYST)
    ds_construct = Construct([role("tandem-repeat cDNA")],
                             name="double-stranded tandem-repeat cDNA")
    double_stranded = graph.product("double_stranded_cdna", ds_construct.name,
                                    MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("second_strand", "Synthesize the second cDNA strand", StepKind.EXTENSION,
                ((tandem, InputRole.TEMPLATE),
                 (dna_polymerase, InputRole.CATALYST)), (double_stranded,),
                reason="Step 25 uses the NEBNext second-strand module at 16 °C for 2.5 h; the protocol does not describe how the second strand is primed.")
    repair = graph.input("end_repair_mix", "NEBNext end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired_cdna", "blunt tandem-repeat cDNA",
                             MaterialKind.MOLECULAR, state=ds_construct)
    graph.react("end_repair", "Blunt the double-stranded cDNA", StepKind.END_REPAIR,
                ((double_stranded, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("a_tailing_enzyme", "Klenow fragment (3′→5′ exo−)", MaterialKind.CATALYST)
    tailed = graph.product("tailed_cdna", "dA-tailed tandem-repeat cDNA",
                           MaterialKind.MOLECULAR,
                           state=repair_and_dA_tail(ds_construct))
    graph.react("a_tail", "Add one 3-prime dA to each cDNA end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("truseq_adapters", "TruSeq indexed adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "NEBNext Quick T4 DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_library", "TruSeq-ligated CirSeq library",
                            MaterialKind.MOLECULAR, state=library)
    graph.react("adapter_ligation", "Ligate TruSeq indexed adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=TRUSEQ_DA_LIGATION)
    ligated_rule = graph.input("ligated_size_selection",
                               "10% urea-PAGE excision of 450–600-nt ligated DNA",
                               MaterialKind.SELECTION)
    ligated_selected = graph.product("selected_ligated_library",
                                     "450–600-nt adapter-ligated library",
                                     MaterialKind.MOLECULAR, state=library)
    ligated_rejected = graph.product("rejected_ligated_sizes",
                                     "ligation products outside the excised 450–600-nt band",
                                     MaterialKind.WASTE)
    graph.react("select_ligated", "Excise 450–600-nt adapter-ligated DNA (Steps 36–45)",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE), (ligated_rule, InputRole.SELECTION)),
                (ligated_selected, ligated_rejected))
    graph.discard(ligated_rejected, "left in the urea-PAGE gel; excised fragments must be at least 435 bases (135 adapter bases plus at least 300 insert bases)")
    pcr_primers = graph.input("pcr_primer_cocktail", "TruSeq PCR primer cocktail (FC-121-4001), P5-side member",
                              MaterialKind.OLIGO, inferred=True)
    pcr_primers_p7 = graph.input("pcr_primer_cocktail_p7",
                                 "TruSeq PCR primer cocktail (FC-121-4001), P7-side member",
                                 MaterialKind.OLIGO, inferred=True)
    phusion = graph.input("phusion", "Phusion high-fidelity DNA polymerase", MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "amplified CirSeq library",
                              MaterialKind.MOLECULAR, state=library)
    graph.react("library_pcr", "Amplify 15 cycles with Phusion (Steps 46–47)", StepKind.PCR,
                ((ligated_selected, InputRole.TEMPLATE), (pcr_primers, InputRole.STRUCTURAL),
                 (pcr_primers_p7, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    library, reason="The protocol names the FC-121-4001 PCR primer cocktail (0.5 µM each primer if ordered separately) but does not print its sequences; splitting it into P5- and P7-side members is a model inference and the vendor original was not inspected in this pass."))
    page_rule = graph.input("native_page_purification", "7.5% native PAGE excision",
                            MaterialKind.SELECTION)
    final = graph.product("sequencing_library", "PAGE-purified CirSeq library",
                          MaterialKind.LIBRARY, state=library)
    page_rejected = graph.product("rejected_amplified_sizes",
                                  "amplified DNA outside the excised band",
                                  MaterialKind.WASTE)
    graph.react("purify_amplified", "Gel-purify the amplified library (Steps 48–53)",
                StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (page_rule, InputRole.SELECTION)),
                (final, page_rejected))
    graph.discard(page_rejected, "left in the 7.5% PAGE gel (Step 52 repeats Steps 42–44)")
    graph.sequence(final, SequencingPlan(
        library, (primers[0],), ("Read 1",),
        READ_LENGTHS, {"Read 1": "include"}))
    graph.declare_reviewed_scope((
        "viral_rna", "zinc_fragmentation", "fragment", "size_selection", "select_size",
        "polynucleotide_kinase", "phosphorylate", "rna_ligase", "circularize",
        "random_rt_primer", "reverse_transcriptase", "rnase_h", "rolling_circle_rt",
        "second_strand_polymerase", "second_strand",
        "end_repair_mix", "end_repair", "a_tailing_enzyme", "a_tail",
        "truseq_adapters", "adapter_ligase", "adapter_ligation", "ligated_library",
        "ligated_size_selection", "selected_ligated_library", "rejected_ligated_sizes",
        "select_ligated", "pcr_primer_cocktail", "pcr_primer_cocktail_p7", "phusion", "amplified_library",
        "library_pcr", "native_page_purification", "rejected_amplified_sizes",
        "purify_amplified", "sequencing_library",
    ), source)
    return graph.require_publishable()
