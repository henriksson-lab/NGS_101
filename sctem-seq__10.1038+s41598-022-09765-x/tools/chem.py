"""Protocol-owned source, molecular workflow and reaction graph."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, complement_segments, feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

DEFINITION = SimpleNamespace(
    title='scTEM-seq',
    doi='10.1038/s41598-022-09765-x',
    citation='Hunt et al., Scientific Reports (2022)',
    summary='Bisulfite-converted single-cell DNA is amplified with indexed SINE Alu or LINE-1 primers while RNA from the same cell is processed in parallel.',
    endpoint='Paired reads span targeted transposable-element CpGs; internal and adapter indexes identify the cell library.',
    notes='01_sctem_seq.html',
)
# Published Table S1 oligos; Table S2 A1 selects i501_F + i701_R.
# This is an explicit second-generation Alu example, not the LINE-1 or first-generation route.
ALU_LEFT_ARM = "CTACACGACGCTCTTCCGATCT"
ALU_RIGHT_ARM = "CAGACGTGTGCTCTTCCGATCT"
ALU_FORWARD_SITE = "AAAAATATAAAAAATTAGTTGGG"
ALU_REVERSE_SITE = "CCCAAACTAAAATACAATAAC"
ALU_I501 = "TAGATCGC"
ALU_I701 = "TCGCCTTA"
ALU_A1_FORWARD = ALU_LEFT_ARM + ALU_I501 + ALU_FORWARD_SITE
ALU_A1_REVERSE = ALU_RIGHT_ARM + ALU_I701 + ALU_REVERSE_SITE

def alu_primer_parts(reverse=False):
    # Stagger is zero for these two printed A1 members. It is not a UMI.
    return [Segment("partial adapter", ALU_RIGHT_ARM if reverse else ALU_LEFT_ARM),
            Segment("internal well index", ALU_I701 if reverse else ALU_I501, "cbc", length_bp=8,
                    feature=feature("sctem_internal_right" if reverse else "sctem_internal_left",
                                    "cell_barcode", "unknown")),
            Segment("Alu binding region", ALU_REVERSE_SITE if reverse else ALU_FORWARD_SITE)]

def alu_amplicon_parts():
    # Polymerase installs the printed primer sequences at the successful candidate's ends.
    left = alu_primer_parts()
    right = complement_segments(alu_primer_parts(reverse=True))
    return [*left, role("bisulfite-converted Alu interval; variable genomic target"), *right]

PAYLOAD = tuple(alu_amplicon_parts())

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("The DNA diagram explicitly selects the published second-generation SINE Alu A1 primer pair "
          "(Table S1 i501_F/i701_R, Table S2 well A1). Its zero-stagger members and two eight-base "
          "internal well indexes are distinct from the later NEBNext indexes. Variable genomic Alu "
          "targets, selected NEBNext oligos and the separate matched Smart-seq2 route remain model limits; "
          "the first-generation and LINE-1 ordering primers below are a separate source inventory.")
FINAL_LIBRARY = Construct([role("NEBNext outer left adapter/index region; unassigned", molecular_feature=feature("sctem_outer_left_index", "sample_index", "unknown")),
                           *PAYLOAD, role("NEBNext outer right adapter/index region; unassigned", molecular_feature=feature("sctem_outer_right_index", "sample_index", "unknown"))],
                          name="scTEM second-generation Alu A1 candidate endpoint")
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "DNA uses paired MiSeq reads and RNA uses NextSeq, but their exact primer assignments and acquired cycle counts are not declared in this model."

def oligos():
    return ("Second-generation SINE Alu A1 example forward, 5′-" + ALU_A1_FORWARD + "-3′",
            "Second-generation SINE Alu A1 example reverse, 5′-" + ALU_A1_REVERSE + "-3′",'Sine.Alu_F — 5′-CTACACGACGCTCTTCCGATCTAAAAATATAAAAAATTAGTTGGG-3′',
            'Sine.Alu_R — 5′-CAGACGTGTGCTCTTCCGATCTNNNNNNNNCCCAAACTAAAATACAATAAC-3′',
            'LINE.L1_F — 5′-CTACACGACGCTCTTCCGATCTTTTTGAGTTAGGTGTGGGATATA-3′',
            'LINE.L1_R — 5′-TGCTGAACCGCTCTTCCGATCTNNNNNNNNAAAATCAAAAAATTCCCTTTC-3′')


def sorted_cell_population():
    """Independent DNA and RNA populations share a cell, not a covalent strand."""
    sc=Scene()
    sc.strand("genomic DNA population",[role("genomic DNA")],label="one sorted cell: genomic DNA")
    sc.strand("mRNA population",[role("mRNA")],label="one sorted cell: mRNA")
    return sc


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one sorted cell containing genomic DNA and mRNA",
        MaterialKind.BIOLOGICAL,
        state=sorted_cell_population(),
        primary=True)
    dna = graph.product("dna_branch", "single-cell genomic DNA branch",
                        MaterialKind.MOLECULAR, state=Construct([role("single-cell genomic DNA")], name="DNA branch"))
    rna = graph.product("rna_branch", "matched single-cell mRNA branch",
                        MaterialKind.MOLECULAR, state=molecule("matched mRNA branch", (role("single-cell mRNA"),)))
    graph.react("separate_modalities", "Separate genomic DNA and mRNA",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    bisulfite = graph.input("bisulfite", "EZ-96 DNA Methylation-Direct MagPrep Kit D5054, half manufacturer volumes",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_dna", "bisulfite-converted transposable-element templates",
        MaterialKind.MOLECULAR,
        state=Construct([role("bisulfite-converted genomic Alu target population")], name="converted gDNA PCR substrate"))
    graph.react("bisulfite_conversion", "Bisulfite-convert genomic DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    te_forward = graph.input("te_forward", "Sine.Alu_i501_F, selected Table S2 A1 example",
                             MaterialKind.OLIGO, state=Construct(alu_primer_parts()))
    te_reverse = graph.input("te_reverse", "Sine.Alu_i701_R, selected Table S2 A1 example",
                             MaterialKind.OLIGO, state=Construct(alu_primer_parts(reverse=True)))
    te_polymerase = graph.input("te_polymerase", "KAPA HiFi HotStart Uracil+ polymerase",
                                MaterialKind.CATALYST)
    amplicons = graph.product("te_amplicons", "scTEM partial-adapter TE amplicons",
                              MaterialKind.MOLECULAR, state=Construct(PAYLOAD, name="Alu A1 partial-adapter PCR candidate"))
    graph.react("te_pcr", "Amplify converted Alu with the selected A1 example pair (35 cycles)",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (te_forward, InputRole.STRUCTURAL),
                 (te_reverse, InputRole.STRUCTURAL), (te_polymerase, InputRole.CATALYST)),
                (amplicons,), geometry=PCRGeometry.unavailable(
                    converted.state,
                    reason="The printed A1 primers are source-positive, but exact sites across the genomic bisulfite-converted Alu population are not declared. Primer-derived candidate ends do not establish genomic target coordinates or amplification specificity."))
    cleanup = graph.input("ampure_cleanup", "1.2x AMPure purification and normalized pooling", MaterialKind.SELECTION)
    pooled = graph.product("pooled_amplicons", "purified pooled TE amplicons", MaterialKind.MOLECULAR, state=amplicons.state)
    wash = graph.product("ampure_wash", "AMPure supernatant and wash fraction", MaterialKind.WASTE)
    graph.react("purify_pool", "Purify TE amplicons (1.2x AMPure), normalize and pool", StepKind.CLEANUP,
                ((amplicons, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)), (pooled, wash), changes_molecule=False)
    graph.discard(wash, "physical AMPure purification fraction")
    index_primers = graph.input("nebnext_indexes", "NEBNext dual-index oligos E7780S", MaterialKind.OLIGO)
    index_reverse = graph.input("nebnext_reverse", "NEBNext opposite-end dual-index oligo E7780S", MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "KAPA HiFi HotStart polymerase", MaterialKind.CATALYST)
    dna_library = graph.product("dna_library", "dual-indexed scTEM TE library", MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Add adapters and dual indexes to pooled TE amplicons (5 cycles)", StepKind.PCR,
                ((pooled, InputRole.TEMPLATE), (index_primers, InputRole.STRUCTURAL), (index_reverse, InputRole.STRUCTURAL), (index_polymerase, InputRole.CATALYST)),
                (dna_library,), geometry=PCRGeometry.unavailable(pooled.state, reason="Selected NEBNext E7780S primer species and their binding footprints have not been mapped on this actual partial-adapter Alu PCR substrate."))
    postbeads = graph.input("post_index_beads", "0.9X AMPure XP", MaterialKind.SELECTION)
    purified_library = graph.product("dna_purified_library", "purified dual-indexed Alu A1 candidate library", MaterialKind.LIBRARY, state=final)
    postwaste = graph.product("post_index_waste", "post-index-PCR supernatant and washes", MaterialKind.WASTE)
    graph.react("post_index_cleanup", "Purify indexed TE libraries (0.9X AMPure), normalize and combine", StepKind.CLEANUP,
                ((dna_library, InputRole.TEMPLATE), (postbeads, InputRole.SELECTION)), (purified_library, postwaste), changes_molecule=False)
    graph.discard(postwaste, "removed during source-stated 0.9X AMPure cleanup")
    graph.sequence(purified_library, SequencingPlan(library=final, acquired_roles=("Read 1", "Read 2"),
        limitation="DNA branch: paired MiSeq reads are source-stated, but acquired cycle counts and selected NEBNext primer/site assignments are not declared here. The two internal eight-base well indexes are not the outer kit indexes."))
    # Published Methods points to refs 9/17 for the independent matched RNA route.
    # No shared RT/PCR endpoint is used to invent a selected RNA implementation.
    rna_final = Construct([role("Smart-seq2 transcriptome left library end; not yet declared"),
                           role("full-length transcript-derived insert"),
                           role("Smart-seq2 transcriptome right library end; not yet declared")],
                          name="scTEM matched transcriptome endpoint roles")
    rna_library = graph.product("rna_library", "matched scTEM Smart-seq2 transcriptome library", MaterialKind.LIBRARY, state=rna_final)
    graph.react("rna_preparation", "Prepare matched scRNA-seq libraries by the cited Smart-seq2 route", StepKind.MAGIC,
                ((rna, InputRole.TEMPLATE),), (rna_library,),
                reason="Hunt et al. published Methods, Library preparation, delegates matched RNA libraries to refs 9/17. "
                       "The current declaration has not mapped the cited implementation's RT, template switch, cDNA amplification, "
                       "tagmentation, library PCR or sample-index choices; it adds no generic terminal PCR or synthetic cell barcode.")
    graph.sequence(rna_library, SequencingPlan(library=rna_final, acquisition="unspecified",
        limitation="Matched RNA libraries use NextSeq according to published Methods; the acquired read roles/cycles and implementation-specific primers are not declared here."))
    reviewed = ["cell", "separate_modalities", "bisulfite", "bisulfite_conversion",
                "te_forward", "te_reverse", "te_polymerase", "te_pcr", "dna_library",
                "ampure_cleanup", "purify_pool", "nebnext_indexes", "nebnext_reverse", "index_polymerase", "index_pcr",
                "post_index_beads", "post_index_cleanup", "dna_purified_library", "rna_preparation", "rna_library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
