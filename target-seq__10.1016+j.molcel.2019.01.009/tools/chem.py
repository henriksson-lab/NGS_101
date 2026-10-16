from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, branches, nextera_library, single, sp, token
from chemdraw import Segment, revcomp

from rna_cell_assays import chain

def target_workflow():
    start = branches("one FACS-indexed cell", (
        ("genomic molecules", (token("genomic DNA"),)),
        ("RNA molecules", (token("poly(A) RNA"),))))
    preamp = branches("one-tube targeted RT-PCR", (
        ("gDNA amplicons", (token("targeted genomic loci"),)),
        ("cDNA amplicons", (token("targeted coding loci"),)),
        ("transcriptome", (token("Smart-seq+ cDNA"),))))
    return chain("single cell", start, (
        ("Lyse mildly, reverse-transcribe, and multiplex-pre-amplify", preamp,
         "Genomic and transcript-specific primers amplify mutation loci alongside whole-transcriptome cDNA."),
        ("Split the common pre-amplification product", branches("matched cell", (
            ("whole-transcriptome branch", (token("Smart-seq+ cDNA"),)),
            ("targeted genotype branch", (token("gDNA and cDNA target amplicons"),)))),
         "Both branches originate from the same physical well and FACS record."),
        ("Tagment the transcriptome; perform two nested PCRs on genotype targets",
         branches("two sequencing libraries", (
             ("Nextera transcriptome", (token("Nextera arm"), token("transcript cDNA"), token("Nextera arm"))),
             ("CS1/CS2 genotype", (token("CS1 site"), token("mutation locus"), token("CS2 site"))))),
         "PCR1 adds CS1/CS2 tails; PCR2 adds the cell-specific Illumina-compatible index."),
    ))


TARGET_R1 = "ACACTGACGACATGGTTCTACA"
TARGET_R2 = "TACGGTAGCAGAGACTTGGTCT"

def target_genotype_library() -> tuple[Construct, tuple[sp.SeqPrimer, ...]]:
    lib = Construct([token("P5 and cell index"), Segment("CS1", TARGET_R1, "r1"),
                     token("targeted gDNA or cDNA locus", length=40),
                     Segment("CS2 reverse complement", revcomp(TARGET_R2), "r2"),
                     token("P7")], name="TARGET-seq genotyping library")
    primers = (sp.custom("Read 1", "TARGET-seq CS1-seq", TARGET_R1,
                         "Rodriguez-Meira et al. 2019 Key Resources"),
               sp.custom("Read 2", "TARGET-seq CS2-seq", TARGET_R2,
                         "Rodriguez-Meira et al. 2019 Key Resources"))
    problems = sp.verify(lib, primers, required_roles=("Read 1", "Read 2"))
    if problems:
        raise ValueError("TARGET-seq primer geometry: " + "; ".join(problems))
    return lib, primers


DEFINITION = Definition("TARGET-seq — targeted genotyping and transcriptomics", "10.1016/j.molcel.2019.01.009",
        "Rodriguez-Meira et al., <i>Molecular Cell</i> (2019)", "01_target-seq.html",
        "Targeted gDNA and coding-DNA loci are co-amplified with Smart-seq+ whole-transcriptome cDNA from one FACS-indexed cell, then split into genotype and expression libraries.",
        "The genotype library uses published CS1/CS2 custom Read 1/Read 2 primers (151 cycles each) plus a 10-cycle cell index; its matched transcriptome branch is Nextera XT.")

d = DEFINITION
final, primers = target_genotype_library()
reads = {'Read 1': 151, 'Read 2': 151}
OLIGOS = ("<b>Smart-seq+ TSO-LNA</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′ (2019 Key Resources, JATS line 479)",
                   "<b>ISPCR</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′ (2019 Key Resources, JATS line 481)",
                   "<b>CS1-seq custom Read 1 primer</b> — 5′-A+CA+CTG+ACGACATGGTTCTACA-3′",
                   "<b>CS2-seq custom Read 2 primer</b> — 5′-T+AC+GGT+AGCAGAGACTTGGTCT-3′",
                   "The published ordering notation uses + before an LNA base; the model uses the corresponding unmodified binding sequences.")
workflow_fn = target_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def target_seq_graph(definition, final: Construct,
                     primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """TARGET-seq graph with matched transcriptome and genotype libraries."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input("single_cell", "one FACS-indexed cell", MaterialKind.BIOLOGICAL,
                       state=target_workflow().initial, primary=True)
    rt_primers = graph.input("rt_primer_pool", "targeted and whole-transcriptome RT primers",
                             MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    copied = graph.product("copied_cell", "cell lysate with cDNA and genomic targets",
                           MaterialKind.MOLECULAR,
                           state=branches("copied single-cell material", (
                               ("genomic", (token("genomic DNA"),)),
                               ("cDNA", (token("target and transcriptome cDNA"),)))))
    graph.react("reverse_transcription", "Reverse-transcribe cellular RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cell, InputRole.TEMPLATE), (rt_primers, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    preamp_forward = graph.input("preamp_forward_pool", "multiplex preamplification forward pool",
                                 MaterialKind.OLIGO)
    preamp_reverse = graph.input("preamp_reverse_pool", "multiplex preamplification reverse pool",
                                 MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "preamplification/PCR polymerase",
                             MaterialKind.CATALYST)
    common = graph.product("common_preamp", "one-tube targeted and transcriptome preamplification",
                           MaterialKind.MOLECULAR,
                           state=branches("one-tube targeted RT-PCR", (
                               ("gDNA", (token("targeted genomic loci"),)),
                               ("cDNA", (token("targeted coding loci"),)),
                               ("transcriptome", (token("Smart-seq+ cDNA"),)))))
    graph.react("multiplex_preamp", "Multiplex-pre-amplify genomic and cDNA targets",
                StepKind.PCR,
                ((copied, InputRole.TEMPLATE),
                 (preamp_forward, InputRole.STRUCTURAL),
                 (preamp_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (common,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex pool contains locus-specific pairs rather than one immutable amplicon."))
    transcriptome = graph.product("transcriptome_branch", "Smart-seq+ transcriptome branch",
                                  MaterialKind.MOLECULAR,
                                  state=single("transcriptome branch", (token("Smart-seq+ cDNA"),)))
    genotype = graph.product("genotype_branch", "targeted gDNA/cDNA amplicon branch",
                             MaterialKind.MOLECULAR,
                             state=single("genotype branch", (
                                 token("gDNA and cDNA target amplicons"),)))
    graph.react("split_preamp", "Split the common preamplification product",
                StepKind.SPLIT, ((common, InputRole.TEMPLATE),),
                (transcriptome, genotype))

    loaded = graph.input("loaded_nextera_adapters", "loaded Nextera XT adapters",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_transcriptome", "tagmented transcriptome cDNA",
                           MaterialKind.MOLECULAR, state=transcriptome.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment_transcriptome", "Tagment the transcriptome branch",
                StepKind.TAGMENTATION,
                ((transcriptome, InputRole.TEMPLATE),
                 (loaded, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    gap_fill = graph.input("gap_fill_polymerase", "Nextera gap-fill polymerase",
                           MaterialKind.CATALYST)
    filled = graph.product("gap_filled_transcriptome", "gap-filled tagmented transcriptome",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (gap_fill, InputRole.CATALYST)), (filled,))
    nx_forward = graph.input("nextera_forward_primer", "Nextera XT indexed forward primer",
                             MaterialKind.OLIGO)
    nx_reverse = graph.input("nextera_reverse_primer", "Nextera XT indexed reverse primer",
                             MaterialKind.OLIGO)
    transcriptome_final, transcriptome_primers = nextera_library(
        [token("transcript cDNA", length=44)], "TARGET-seq transcriptome library")
    transcriptome_library = graph.product(
        "transcriptome_library", "TARGET-seq transcriptome library",
        MaterialKind.LIBRARY, state=transcriptome_final)
    graph.react("transcriptome_pcr", "Index and amplify the transcriptome library",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (nx_forward, InputRole.STRUCTURAL),
                 (nx_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (transcriptome_library,),
                geometry=PCRGeometry.unavailable(
                    transcriptome_final,
                    reason="The Nextera XT index choice varies by cell; the endpoint geometry is canonical."))
    graph.sequence(transcriptome_library, SequencingPlan(
        transcriptome_final, tuple(transcriptome_primers),
        tuple(p.role for p in transcriptome_primers)))

    nested1_forward = graph.input("nested1_forward_pool", "first nested target-PCR forward pool",
                                  MaterialKind.OLIGO)
    nested1_reverse = graph.input("nested1_reverse_pool", "first nested target-PCR reverse pool",
                                  MaterialKind.OLIGO)
    nested = graph.product("nested_target_amplicons", "CS1/CS2-tailed target amplicons",
                           MaterialKind.MOLECULAR, state=genotype.state)
    graph.react("nested_pcr1", "First nested PCR adds CS1/CS2 tails", StepKind.PCR,
                ((genotype, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex target pool contains locus-specific nested primer pairs."))
    cs1 = graph.input("cs1_index_primer", "cell-indexed CS1 PCR primer",
                      MaterialKind.OLIGO)
    cs2 = graph.input("cs2_index_primer", "cell-indexed CS2 PCR primer",
                      MaterialKind.OLIGO)
    genotype_library = graph.product("genotype_library", "TARGET-seq genotype library",
                                     MaterialKind.LIBRARY, state=final)
    graph.react("nested_pcr2", "Second nested PCR adds the cell-specific index",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE), (cs1, InputRole.STRUCTURAL),
                 (cs2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (genotype_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The indexed CS1/CS2 amplification family varies by cell."))
    graph.sequence(genotype_library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        read_lengths={"Read 1": 151, "Read 2": 151}))
    graph.declare_reviewed_scope((
        "single_cell", "rt_primer_pool", "reverse_transcriptase", "reverse_transcription",
        "preamp_forward_pool", "preamp_reverse_pool", "pcr_polymerase",
        "multiplex_preamp", "split_preamp", "loaded_nextera_adapters", "tn5",
        "tagment_transcriptome", "gap_fill_polymerase", "fill_tagmentation_gaps",
        "nextera_forward_primer", "nextera_reverse_primer",
        "transcriptome_pcr", "transcriptome_library", "nested1_forward_pool",
        "nested1_reverse_pool", "nested_pcr1", "cs1_index_primer", "cs2_index_primer",
        "nested_pcr2", "genotype_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'target-seq'
    return target_seq_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
