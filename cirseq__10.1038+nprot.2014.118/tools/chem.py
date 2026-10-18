"""Protocol-owned source, oligos, workflow and reaction graph."""
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from endprep import repair_and_dA_tail
from protocol_model import closed_circle, molecule, role
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

# Local endpoint records only source-supported roles. The complete historical
# adapters are not selected by a canonical platform constructor.
FINAL_LIBRARY = Construct([
    role("historical TruSeq left adapter; exact species unassigned"),
    role("tandem-repeat cDNA"),
    Segment("source dA junction", "A"),
    role("historical TruSeq right insert-facing arm; unassigned"),
    role("sample index; molecular extent and member unassigned", tag="cbc",
         molecular_feature=feature("cirseq_sample_index", "sample_index", "unknown")),
    role("historical TruSeq right terminal arm; unassigned"),
], name="CirSeq source-limited adapter endpoint")
SEQ_PRIMERS = ()
TRUSEQ_DA_LIGATION = AdapterLigationGeometry.unavailable(
    "Step 31 dA-tails the cDNA and Step 35 ligates TruSeq indexed adapters. The later 2015 FC-121-4001 guide confirms a 3′-T family end, but complete historical adapter identity/revision has not been assigned. The native graph does not import a canonical Y-stem or later primer shell.")
TITLE = 'CirSeq — circular resequencing of RNA virus populations'
NOTES = '01_cirseq.html'
SOURCE = '<a href="https://doi.org/10.1038/nprot.2014.118">Acevedo and Andino, <i>Nature Protocols</i> (2014)</a>.'
SUMMARY = 'Short viral RNA fragments are circularized and copied repeatedly by rolling-circle reverse transcription, placing at least three linked template copies inside one Illumina read.'
CAVEAT = ('The original protocol names TruSeq indexed adapters and PCR primer cocktail FC-121-4001, but this model has not reconstructed the complete historical adapter/primer species. '
          'Step 43 states 135 combined adapter bases; the former generic canonical shell did not match that source value and is no longer used. '
          'The 2015 kit guide is a comparator: it lists six reported index bases plus a seventh adapter base and a 3′-T ligating end, without establishing full 2014 lot equivalence. '
          'The native graph retains the original RNA 85–100-nt window and the recommended 300-nt single-read profile, with primer sites unassigned.')
READ_LENGTHS = {'Read 1': 300}
SEQUENCING_UNAVAILABLE = 'The selected source recommendation uses a 300-nt single read; the current model has not reconstructed the historical TruSeq primer/adapter sites.'
SEQUENCING_ENDING = SEQUENCING_UNAVAILABLE
FINAL_CAPTION = 'TruSeq adapters flank a double-stranded cDNA insert containing tandem copies of one circular RNA fragment.'
SEQUENCING_INTRO = "The selected source recommendation is a 300-nt single read on MiSeq v2. The exact adapter/primer reconstruction remains open; an indexed-adapter name does not establish an acquired index cycle plan or a paired Read 2 run."
def oligos():
    return ()

def ligated_candidate_scene():
    scene = Scene()
    scene.strand("candidate", list(FINAL_LIBRARY), label="ligated candidate strand; opposite structure not reconstructed")
    scene.junction("candidate", "historical TruSeq left adapter; exact species unassigned", "tandem-repeat cDNA", "adapter ligation")
    scene.junction("candidate", "source dA junction", "historical TruSeq right insert-facing arm; unassigned", "adapter ligation")
    return scene

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
        state=molecule("size-selected RNA", (role("RNA fragment selected within 85–100nt; variable length"),)))
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
                            MaterialKind.OLIGO, state=Construct([Segment("random hexamer", "N" * 6, placeholder=True, length_bp=6)]))
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
                reason="Step 25 uses the NEBNext second-strand module at 16 °C for 2.5 h; its priming mechanism remains unassigned pending applicable module-original review.")
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
    adapters = graph.input("truseq_adapters", "historical TruSeq indexed adapters (FC-121-4001); complete species unassigned", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "NEBNext Quick T4 DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_library", "TruSeq-ligated CirSeq library",
                            MaterialKind.MOLECULAR, state=ligated_candidate_scene())
    graph.react("adapter_ligation", "Ligate TruSeq indexed adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=TRUSEQ_DA_LIGATION)
    ligated_rule = graph.input("ligated_size_selection",
                               "10% urea-PAGE excision of 450–600-nt ligated DNA",
                               MaterialKind.SELECTION)
    ligated_selected = graph.product("selected_ligated_library",
                                     "450–600-nt adapter-ligated library",
                                     MaterialKind.MOLECULAR, state=ligated.state)
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
                    ligated_selected.state, reason="The original protocol names the FC-121-4001 cocktail (0.5 µM each primer if ordered separately); the model has not assigned the historical primer species or copied adapter ends on this actual ligation product. Splitting the cocktail into opposite-end roles remains inferred; the 2015 comparator does not supply the missing historical model state."))
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
    graph.sequence(final, SequencingPlan(library=library, acquired_roles=("Read 1",),
        read_lengths=READ_LENGTHS, read_uses={"Read 1": "include"},
        limitation="The selected recommended MiSeq profile uses single 300-nt reads to contain tandem copies; the exact historical sequencing primer/site is not reconstructed. No paired Read 2 or index acquisition is assigned from a canonical adapter shell."))
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
