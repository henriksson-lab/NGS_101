from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, branches, nextera_library, single, sp, token
from chemdraw import Segment, revcomp

from rna_cell_assays import chain
from protocol_graph import SourceClaim, SourceDiscrepancy

def target_workflow():
    start = branches("one FACS-indexed cell", (
        ("genomic molecules", (token("genomic DNA"),)),
        ("RNA molecules", (token("poly(A) RNA"),))))
    preamp = branches("one-tube targeted RT-PCR", (
        ("gDNA amplicons", (token("targeted genomic loci"),)),
        ("cDNA amplicons", (token("targeted coding loci"),)),
        ("transcriptome", (token("ISPCR site"), token("Smart-seq+ cDNA"), token("ISPCR site")))))
    return chain("single cell", start, (
        ("Lyse with protease; reverse-transcribe with template switching; multiplex-pre-amplify", preamp,
         "Table S1b: Triton X-100/Qiagen protease lysis, 72 °C 15 min; SMARTScribe with oligo(dT)-ISPCR, mRNA target primers and TSO-LNA; "
         "SeqAmp PCR with ISPCR plus gDNA and cDNA target primers (20 cycles cell lines, 22 HSPCs)."),
        ("Split the unpurified cDNA-amplicon mix", branches("matched cell", (
            ("whole-transcriptome branch", (token("Smart-seq+ cDNA"),)),
            ("targeted genotype branch", (token("gDNA and cDNA target amplicons"),)))),
         "15 of 25 µL go to 0.6:1 AMPure XP purification for Nextera XT; a 1.5 µL aliquot is the genotyping input."),
        ("Tagment the transcriptome; perform two PCRs on genotype targets",
         branches("two sequencing libraries", (
             ("Nextera transcriptome", (token("Nextera arm"), token("transcript cDNA"), token("Nextera arm"))),
             ("CS1/CS2 genotype", (token("PE1 arm"), token("CS1 site"), token("mutation locus"),
                                   token("CS2 site"), token("PE2 arm"))))),
         "PCR1 (KAPA 2G Robust, 20 cycles) adds CS1/CS2 tails; PCR2 (FastStart HF, 10 cycles) adds PE1/PE2 arms with a 10-nt "
         "Fluidigm Access Array single-direction barcode, whose arm is not stated in the inspected TARGET-seq sources."),
    ))


TARGET_R1 = "ACACTGACGACATGGTTCTACA"
TARGET_R2 = "TACGGTAGCAGAGACTTGGTCT"

def target_genotype_library() -> tuple[Construct, tuple[sp.SeqPrimer, ...]]:
    # PCR2 attaches PE1/PE2 arms carrying one 10-nt Access Array barcode. The inspected
    # TARGET-seq sources do not say which arm carries it, so neither arm token claims it.
    lib = Construct([token("PE1 arm (Access Array primer)"), Segment("CS1", TARGET_R1, "r1"),
                     token("targeted gDNA or cDNA locus", length=40),
                     Segment("CS2 reverse complement", revcomp(TARGET_R2), "r2"),
                     token("PE2 arm (Access Array primer)")], name="TARGET-seq genotyping library")
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
        "The genotype library uses published CS1/CS2 custom Read 1/Read 2 primers (151 cycles each) plus a 10-cycle index read primed by a CS1rc-seq/CS2rc-seq mix. "
        "The 10-nt Fluidigm Access Array single-direction barcode lies in one PE arm; the inspected TARGET-seq sources do not say which, so the drawn arms do not place it. "
        "Its matched transcriptome branch is Nextera XT, read as 75-base single-end NextSeq reads.")

d = DEFINITION
final, primers = target_genotype_library()
reads = {'Read 1': 151, 'Read 2': 151}
OLIGOS = ("<b>OligodT-ISPCR</b> — 5′-AAGCAGTGGTATCAACGCAGAGTAC(T)<sub>30</sub>VN-3′ (printed lower-case, 2019 Key Resources, JATS line 477)",
           "<b>Smart-seq+ TSO-LNA</b> — 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrG+G-3′ (2019 Key Resources, JATS line 479)",
                   "<b>ISPCR</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′ (2019 Key Resources, JATS line 481)",
                   "<b>CS1-seq custom Read 1 primer</b> — 5′-A+CA+CTG+ACGACATGGTTCTACA-3′",
                   "<b>CS2-seq custom Read 2 primer</b> — 5′-T+AC+GGT+AGCAGAGACTTGGTCT-3′",
                   "<b>CS1rc-seq / CS2rc-seq index-read primers</b> — 5′-T+GT+AG+AACCATGTCGTCAGTGT-3′ / 5′-A+GAC+CA+AGTCTCTGCTACCGTA-3′ (JATS lines 496–500; used together for the index read)",
                   "Target-specific RT/pre-amplification and CS1/CS2-tailed PCR1 primers are in Table S2A/S2B; PCR2 primers are the Fluidigm Access Array Barcode Library (sequences not printed in the paper).",
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
    lysis = graph.input("lysis_buffer",
                        "4.2 µL lysis: Triton X-100, RNase inhibitor, dNTPs, oligo(dT)-ISPCR, Qiagen protease (Table S1b)",
                        MaterialKind.CONDITION)
    ercc = graph.input("ercc", "ERCC spike-in RNA", MaterialKind.MOLECULAR)
    heat = graph.input("protease_inactivation", "72 °C 15 min protease heat inactivation (Table S1b)",
                       MaterialKind.CONDITION)
    lysate = graph.product("lysate", "protease-treated single-cell lysate",
                           MaterialKind.MOLECULAR, state=target_workflow().initial)
    graph.react("lyse_and_inactivate", "Lyse with protease, then heat-inactivate at 72 °C",
                StepKind.THERMAL_TREATMENT,
                ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION),
                 (ercc, InputRole.SPIKE_IN), (heat, InputRole.CONDITION)), (lysate,))
    oligo_dt = graph.input("oligo_dt_ispcr", "oligo(dT)-ISPCR (in lysis buffer)", MaterialKind.OLIGO)
    mrna_primers = graph.input("mrna_target_primers", "mRNA target-specific RT primers (Table S2A)",
                               MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "SMARTScribe reverse transcriptase (10.2 U/µL)",
                        MaterialKind.CATALYST)
    first_strand = graph.product("first_strand", "first-strand cDNA on cellular RNA",
                                 MaterialKind.MOLECULAR,
                                 state=branches("first-strand single-cell material", (
                                     ("genomic", (token("genomic DNA"),)),
                                     ("cDNA", (token("oligo(dT) and target-primed first-strand cDNA"),)))))
    graph.react("reverse_transcription", "Reverse-transcribe: 42 °C 90 min, then 10 × (50 °C 2 min, 42 °C 2 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((lysate, InputRole.TEMPLATE), (oligo_dt, InputRole.STRUCTURAL),
                 (mrna_primers, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (first_strand,))
    tso = graph.input("tso_lna", "TSO-LNA 5′-…ACATrGrG+G-3′ (1 µM)", MaterialKind.OLIGO)
    copied = graph.product("copied_cell", "cell lysate with template-switched cDNA and genomic targets",
                           MaterialKind.MOLECULAR,
                           state=branches("copied single-cell material", (
                               ("genomic", (token("genomic DNA"),)),
                               ("cDNA", (token("target and transcriptome cDNA"),)))))
    graph.react("template_switching", "Template-switch onto TSO-LNA during the RT program",
                StepKind.TEMPLATE_SWITCHING,
                ((first_strand, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    ispcr = graph.input("ispcr_primer", "ISPCR primer (50 nM)", MaterialKind.OLIGO)
    gdna_primers = graph.input("gdna_target_primers", "gDNA target-specific primers (Table S2A)",
                               MaterialKind.OLIGO)
    cdna_primers = graph.input("cdna_target_primers", "cDNA target-specific primers (Table S2A)",
                               MaterialKind.OLIGO)
    seqamp = graph.input("seqamp", "SeqAmp DNA polymerase", MaterialKind.CATALYST)
    common = graph.product("common_preamp", "one-tube targeted and transcriptome preamplification",
                           MaterialKind.MOLECULAR,
                           state=branches("one-tube targeted RT-PCR", (
                               ("gDNA", (token("targeted genomic loci"),)),
                               ("cDNA", (token("targeted coding loci"),)),
                               ("transcriptome", (token("Smart-seq+ cDNA"),)))))
    graph.react("multiplex_preamp",
                "SeqAmp PCR, 67 °C annealing: 20 cycles (cell lines) or 22 (HSPCs)",
                StepKind.PCR,
                ((copied, InputRole.TEMPLATE), (ispcr, InputRole.STRUCTURAL),
                 (gdna_primers, InputRole.STRUCTURAL),
                 (cdna_primers, InputRole.STRUCTURAL),
                 (seqamp, InputRole.CATALYST)), (common,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex pool contains locus-specific pairs rather than one immutable amplicon."))
    transcriptome = graph.product("transcriptome_branch", "15 µL Smart-seq+ transcriptome aliquot",
                                  MaterialKind.MOLECULAR,
                                  state=single("transcriptome branch", (token("Smart-seq+ cDNA"),)))
    genotype = graph.product("genotype_branch", "1.5 µL unpurified cDNA-amplicon aliquot",
                             MaterialKind.MOLECULAR,
                             state=single("genotype branch", (
                                 token("gDNA and cDNA target amplicons"),)))
    graph.react("split_preamp", "Split the unpurified cDNA-amplicon mix",
                StepKind.SPLIT, ((common, InputRole.TEMPLATE),),
                (transcriptome, genotype))

    beads = graph.input("ampure_06x", "AMPure XP beads, 0.6:1", MaterialKind.SELECTION)
    purified = graph.product("purified_transcriptome", "bead-purified cDNA in 8 µL EB",
                             MaterialKind.MOLECULAR, state=transcriptome.state)
    bead_waste = graph.product("ampure_06x_supernatant", "0.6:1 supernatant: primers and short products",
                               MaterialKind.WASTE)
    graph.react("purify_transcriptome", "Dilute with 11 µL water and purify with 0.6:1 AMPure XP",
                StepKind.CLEANUP,
                ((transcriptome, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (purified, bead_waste))
    graph.discard(bead_waste, "left in the bead supernatant")
    loaded = graph.input("loaded_nextera_adapters", "loaded Nextera XT adapters",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Nextera XT transposase (one-quarter volume reactions)", MaterialKind.CATALYST)
    tagged = graph.product("tagged_transcriptome", "tagmented transcriptome cDNA",
                           MaterialKind.MOLECULAR, state=transcriptome.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment_transcriptome", "Tagment 700 pg cDNA with ATM, 55 °C 6 min; neutralize with NT",
                StepKind.TAGMENTATION,
                ((purified, InputRole.TEMPLATE),
                 (loaded, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    nextera_pcr = graph.input("nextera_pcr_mix", "Nextera XT kit PCR master mix",
                              MaterialKind.CATALYST)
    filled = graph.product("gap_filled_transcriptome", "gap-filled tagmented transcriptome",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (nextera_pcr, InputRole.CATALYST)), (filled,))
    nx_forward = graph.input("nextera_forward_primer", "Nextera XT Index Kit v2 (Set A/C) i5 primer",
                             MaterialKind.OLIGO)
    nx_reverse = graph.input("nextera_reverse_primer", "Nextera XT Index Kit v2 (Set A/C) i7 primer",
                             MaterialKind.OLIGO)
    transcriptome_final, transcriptome_primers = nextera_library(
        [token("transcript cDNA", length=44)], "TARGET-seq transcriptome library")
    amplified_transcriptome = graph.product(
        "amplified_transcriptome", "indexed Nextera XT transcriptome PCR product",
        MaterialKind.MOLECULAR, state=transcriptome_final)
    graph.react("transcriptome_pcr", "NPM indexing PCR: 72 °C 3 min, then 14 cycles (55 °C annealing)",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (nx_forward, InputRole.STRUCTURAL),
                 (nx_reverse, InputRole.STRUCTURAL),
                 (nextera_pcr, InputRole.CATALYST)), (amplified_transcriptome,),
                geometry=PCRGeometry.unavailable(
                    transcriptome_final,
                    reason="The Nextera XT index choice varies by cell; the endpoint geometry is canonical."))
    final_beads = graph.input("ampure_06x_library", "AMPure XP beads, 0.6:1 (post-PCR)", MaterialKind.SELECTION)
    transcriptome_library = graph.product(
        "transcriptome_library", "TARGET-seq transcriptome library",
        MaterialKind.LIBRARY, state=transcriptome_final)
    library_waste = graph.product("ampure_06x_library_supernatant", "post-PCR 0.6:1 supernatant",
                                  MaterialKind.WASTE)
    graph.react("purify_transcriptome_library", "Dilute 1:1 and purify the indexed library with 0.6:1 AMPure XP",
                StepKind.CLEANUP,
                ((amplified_transcriptome, InputRole.TEMPLATE), (final_beads, InputRole.SELECTION)),
                (transcriptome_library, library_waste))
    graph.discard(library_waste, "left in the bead supernatant")
    graph.annotate(SourceDiscrepancy(
        "transcriptome_library", "Published Nextera XT index kits differ",
        (SourceClaim("Key Resources table (PMC6436961 JATS lines 443–445)",
                     "Nextera XT Index Kit v2 Set A (FC-131-2001) and Set C (FC-131-2003)"),
         SourceClaim("Supplemental Experimental Procedures, full-length protocol step 25 (mmc1 line 7150)",
                     "Nextera XT Index Kit, FC-131-1001")),
        "The model draws generic 8-nt i5/i7 windows and assigns no index set.",
        "Both statements are in the published article and its supplement.",
        "Which index kit was used for which sequenced plate is not stated."))
    graph.sequence(transcriptome_library, SequencingPlan(
        transcriptome_final, tuple(transcriptome_primers),
        tuple(p.role for p in transcriptome_primers if p.role != "Read 2"),
        read_lengths={"Read 1": 75}))

    nested1_forward = graph.input("pcr1_forward_pool", "CS1-tailed target-specific PCR1 primers (Table S2B)",
                                  MaterialKind.OLIGO)
    nested1_reverse = graph.input("pcr1_reverse_pool", "CS2-tailed target-specific PCR1 primers (Table S2B)",
                                  MaterialKind.OLIGO)
    kapa = graph.input("kapa_2g_robust", "KAPA 2G Robust HS ReadyMix", MaterialKind.CATALYST)
    nested = graph.product("nested_target_amplicons", "CS1/CS2-tailed target amplicons",
                           MaterialKind.MOLECULAR, state=genotype.state)
    graph.react("nested_pcr1", "PCR1 adds CS1/CS2 tails: 20 cycles, 60 °C annealing", StepKind.PCR,
                ((genotype, InputRole.TEMPLATE),
                 (nested1_forward, InputRole.STRUCTURAL),
                 (nested1_reverse, InputRole.STRUCTURAL),
                 (kapa, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The multiplex target pool contains locus-specific CS1/CS2-tailed primer pairs."))
    barcode_pe1 = graph.input("access_array_pe1",
                              "Access Array single-direction PE1-side primer (barcode side unstated)",
                              MaterialKind.OLIGO)
    barcode_pe2 = graph.input("access_array_pe2",
                              "Access Array single-direction PE2-side primer (barcode side unstated)",
                              MaterialKind.OLIGO)
    faststart = graph.input("faststart_hf", "FastStart High Fidelity PCR System", MaterialKind.CATALYST)
    indexed = graph.product("indexed_genotype_amplicons", "one cell's barcoded genotype amplicons",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("nested_pcr2", "PCR2 adds PE1/PE2 arms and a 10-nt cell barcode: 10 cycles, 60 °C",
                StepKind.PCR,
                ((nested, InputRole.TEMPLATE), (barcode_pe1, InputRole.STRUCTURAL),
                 (barcode_pe2, InputRole.STRUCTURAL),
                 (faststart, InputRole.CATALYST)), (indexed,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The Access Array barcode primers are not printed in the TARGET-seq sources."))
    pool_beads = graph.input("ampure_08x", "AMPure XP beads, 0.8:1", MaterialKind.SELECTION)
    genotype_library = graph.product("genotype_library", "TARGET-seq genotype library",
                                     MaterialKind.LIBRARY, state=final)
    pool_waste = graph.product("ampure_08x_supernatant", "0.8:1 supernatant: primers and short products",
                               MaterialKind.WASTE)
    graph.react("pool_and_purify", "Pool barcoded cells (Mosquito HTS) and purify with 0.8:1 AMPure XP",
                StepKind.CLEANUP,
                ((indexed, InputRole.TEMPLATE), (pool_beads, InputRole.SELECTION)),
                (genotype_library, pool_waste))
    graph.discard(pool_waste, "left in the bead supernatant")
    graph.sequence(genotype_library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers),
        read_lengths={"Read 1": 151, "Read 2": 151}))
    graph.declare_reviewed_scope((
        "single_cell", "lysis_buffer", "ercc", "protease_inactivation", "lyse_and_inactivate",
        "oligo_dt_ispcr", "mrna_target_primers", "reverse_transcriptase", "reverse_transcription",
        "tso_lna", "template_switching", "ispcr_primer", "gdna_target_primers",
        "cdna_target_primers", "seqamp", "multiplex_preamp", "split_preamp",
        "ampure_06x", "purify_transcriptome", "loaded_nextera_adapters", "tn5",
        "tagment_transcriptome", "nextera_pcr_mix", "fill_tagmentation_gaps",
        "nextera_forward_primer", "nextera_reverse_primer",
        "transcriptome_pcr", "amplified_transcriptome", "ampure_06x_library",
        "purify_transcriptome_library", "transcriptome_library", "pcr1_forward_pool",
        "pcr1_reverse_pool", "kapa_2g_robust", "nested_pcr1", "access_array_pe1", "access_array_pe2",
        "faststart_hf", "nested_pcr2", "ampure_08x", "pool_and_purify", "genotype_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'target-seq'
    return target_seq_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
