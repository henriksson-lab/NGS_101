"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('scTEM-seq', '10.1038/s41598-022-09765-x', 'Hunt et al., Scientific Reports (2022)', 'Bisulfite-converted single-cell DNA is amplified with indexed SINE Alu or LINE-1 primers while RNA from the same cell is processed in parallel.', 'Paired reads span targeted transposable-element CpGs; internal and adapter indexes identify the cell library.', '01_sctem_seq.html')
PAYLOAD = (CELL, role("internal dual index"), role("TE amplicon"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ('Published Supplementary Table S1 prints SINE Alu/LINE-1 primers and '
          'indexed stagger families. This role-level diagram does not yet select '
          'one well-specific member or map the inherited Smart-seq2 and NEBNext '
          'oligos onto the final libraries. The published article and supplements '
          'define this page; superseded preprint primer designs are not used.')
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))), (
        ("Separate genomic DNA and mRNA for parallel processing",
         molecule("matched DNA and RNA branches", (CELL, role("genomic DNA"), role("full-length cDNA"))),
         "The split permits Smart-seq-type transcriptomics alongside targeted methylation."),
        ("Bisulfite-convert genomic DNA",
         molecule("converted DNA", (role("converted SINE Alu / LINE-1 templates"),)),
         "Unmethylated cytosines convert while methylated CpGs remain cytosine."),
        ("Amplify transposable elements with staggered, internally indexed primer pairs",
         molecule("scTEM amplicons", (CELL, role("internal dual index"), role("SINE Alu insert"), SAMPLE)),
         "A 0–5 N stagger and reversed primer orientation diversify the low-complexity amplicon pool."),
    ))


def oligos():
    return ('Sine.Alu_F — 5′-CTACACGACGCTCTTCCGATCTAAAAATATAAAAAATTAGTTGGG-3′',
            'Sine.Alu_R — 5′-CAGACGTGTGCTCTTCCGATCTNNNNNNNNCCCAAACTAAAATACAATAAC-3′',
            'LINE.L1_F — 5′-CTACACGACGCTCTTCCGATCTTTTTGAGTTAGGTGTGGGATATA-3′',
            'LINE.L1_R — 5′-TGCTGAACCGCTCTTCCGATCTNNNNNNNNAAAATCAAAAAATTCCCTTTC-3′')


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one sorted cell containing genomic DNA and mRNA",
        MaterialKind.BIOLOGICAL,
        state=molecule("one sorted cell", (role("genomic DNA"), role("mRNA"))),
        primary=True)
    dna = graph.product("dna_branch", "single-cell genomic DNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    rna = graph.product("rna_branch", "matched single-cell mRNA branch",
                        MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_modalities", "Separate genomic DNA and mRNA",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE),), (dna, rna),
                changes_molecule=False)
    bisulfite = graph.input("bisulfite", "bisulfite conversion condition",
                            MaterialKind.CONDITION)
    converted = graph.product(
        "converted_dna", "bisulfite-converted transposable-element templates",
        MaterialKind.MOLECULAR,
        state=molecule("converted DNA", (
            role("converted SINE Alu / LINE-1 templates"),)))
    graph.react("bisulfite_conversion", "Bisulfite-convert genomic DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((dna, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    te_forward = graph.input("te_forward", "staggered internally indexed TE forward-primer pool",
                             MaterialKind.OLIGO)
    te_reverse = graph.input("te_reverse", "staggered internally indexed TE reverse-primer pool",
                             MaterialKind.OLIGO)
    te_polymerase = graph.input("te_polymerase", "KAPA HiFi HotStart Uracil+ polymerase",
                                MaterialKind.CATALYST)
    amplicons = graph.product("te_amplicons", "scTEM partial-adapter TE amplicons",
                              MaterialKind.MOLECULAR, state=molecule("TE amplicons", PAYLOAD))
    graph.react("te_pcr", "Amplify converted transposable elements (35 cycles)",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (te_forward, InputRole.STRUCTURAL),
                 (te_reverse, InputRole.STRUCTURAL), (te_polymerase, InputRole.CATALYST)),
                (amplicons,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The panel uses staggered internally indexed primer families rather than one universal amplicon."))
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
                (dna_library,), geometry=PCRGeometry.unavailable(final, reason="Selected NEBNext dual-index primer geometry has not been mapped in this role-level model."))
    graph.sequence(dna_library, SequencingPlan(limitation="DNA branch: " + d.endpoint))
    rna_final = role_library(
        "scTEM matched transcriptome library roles",
        (CELL, role("full-length transcript cDNA")), index=SAMPLE)
    reviewed = [
        "cell", "separate_modalities", "bisulfite", "bisulfite_conversion",
        "te_forward", "te_reverse", "te_polymerase", "te_pcr", "dna_library",
        "ampure_cleanup", "purify_pool", "nebnext_indexes", "nebnext_reverse", "index_polymerase", "index_pcr",
    ]
    _rt_then_pcr(graph, "rna", rna, rna_final,
                 "matched RNA branch: Smart-seq-type full-length transcript reads",
                 reviewed,
                 rt_primer_label="SMART-seq2 oligo(dT) RT primer (cited protocol; not printed)",
                 rt_enzyme_label="SMART-seq2 reverse transcriptase (cited protocol)",
                 label="scTEM transcriptome",
                 between_rt_and_pcr=(
                     "Methods: gDNA and RNA are separated by G&T-seq (ref. 45), and matched "
                     "scRNA-seq libraries are prepared 'using the SMART-seq2 protocol' (refs 17 "
                     "and 19). The paper prints no RNA-branch oligos or conditions. The cited "
                     "SMART-seq2 route (template switching, cDNA pre-amplification and "
                     "tagmentation) is not redrawn here."),
                 step_label="Final library PCR of the SMART-seq2 transcriptome library (cited protocol)",
                 forward_label="SMART-seq2 library forward PCR primer (cited protocol; not printed)",
                 reverse_label="SMART-seq2 library reverse PCR primer (cited protocol; not printed)")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
