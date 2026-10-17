from chemdraw import Segment, Workflow
import illumina

# Supplementary Table 5 (PMC6782071 sup_tab_5.xlsx rows 5-6), 5'->3'.
P5_GENERIC = illumina.P5
RPI_X_3PRIME = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
from types import SimpleNamespace
from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, SourceCoverage, SourceItem, StepKind, custom_library, single, token, unknown_region
from rna_cell_assays import feature

from rna_cell_assays import branches, chain

def got_workflow():
    initial = single("10x first-strand cDNA", (
        token("target transcript cDNA"),
        token("10x cell barcode", tag="cbc", length=16,
              molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                        whitelist="10x Chromium chemistry-specific whitelist")),
        unknown_region("chemistry-specific UMI", tag="umi",
              molecular_feature=feature("umi", "umi", "random")),))
    return chain("barcoded 10x cDNA", initial, (
        ("Optionally spike a target primer into RT or cDNA amplification", initial,
         "The paper specifies 5′ RT/additive spike-ins and a 3′ cDNA-amplification spike-in when needed."),
        ("Split amplified cDNA into expression and genotyping branches",
         branches("same cell barcode and UMI", (
             ("expression branch", (token("whole-transcriptome cDNA"),)),
             ("genotype branch", (
                 token("target cDNA"),
                 token("cell barcode", tag="cbc", length=16),
                 unknown_region("chemistry-specific UMI", tag="umi",
                                molecular_feature=feature("umi", "umi", "random")),
             )),
         )),
         "Only a small cDNA aliquot enters targeted genotyping; the remainder follows standard 10x preparation."),
        ("Target PCR: generic forward SI-PCR primer with a staggered locus-specific reverse primer",
         single("GoT amplicon", (token("variant-containing target"),
                                  token("cell barcode", tag="cbc", length=16),
                                  unknown_region("chemistry-specific UMI", tag="umi",
                                      molecular_feature=feature("umi", "umi", "random")))),
         "SI-PCR anneals the partial Read 1 end and preserves cell barcode and UMI; the reverse primer adds a partial Read 2 handle, a 0-3 nt stagger and the gene-specific region. Unincorporated primers are then removed by SPRI."),
        ("Second PCR with P5_generic and an indexed RPI-x primer",
         single("GoT library", (token("sequencing arm"), token("variant-containing target"),
                                 token("cell barcode", tag="cbc", length=16),
                                 unknown_region("chemistry-specific UMI", tag="umi",
                                     molecular_feature=feature("umi", "umi", "random")),
                                 token("indexed sequencing arm"))),
         "P5_generic retains the cell barcode and UMI; RPI-x completes P7 and adds the 8-nt sample index. Read 1 reports cellular identifiers and Read 2 reaches the selected variant."),
    ))


DEFINITION = Definition("GoT — Genotyping of Transcriptomes", "10.1038/s41586-019-1367-0",
        "Nam et al., <i>Nature</i> (2019)", "01_got.html",
        "A targeted amplicon branch is made from 10x cDNA while retaining its cell barcode and UMI, linking a selected transcribed variant to the cell's whole transcriptome.",
        "Read 1 is 26 cycles (3′ v2 and 5′) or 28 cycles (3′ v3) for cell barcode and UMI, the 8-cycle index read selects the library, and the 98- or 130-cycle Read 2 reaches the targeted variant (HiSeq 2500 spike-in or separate MiSeq run).",
        "The defining paper supports 10x 3′ v2/v3 and 5′ inputs. The schematic shows the common linear identifier-preserving branch; exact UMI length follows the chosen 10x chemistry. For 3′ v2 only, the GoT aliquot (3 µL; 10 µL for 3′ v3 and 5′) is first re-amplified 3–4 cycles with KAPA HiFi HotStart ReadyMix and 10x primer mix; this optional step is not drawn. The 10x SI-PCR primer sequence is not printed in the paper or Supplementary Table 5. Circularization GoT is not drawn; its published handle assignment differs between sources: Methods give Fw-primer#2 AGGUCAGTCU and Rv-primer#2 AGACUGACCU, whereas Supplementary Table 5 lists Handle PCR2_fw AGACUGACCU and PCR2_rv AGGUCAGTCU; Methods print Rv-primer#1 without the 3′ phosphorothioate shown in Table 5 (…CCGATC*T). No choice is modelled.")

d = DEFINITION
final = custom_library('GoT targeted library', (token('cell barcode', tag='cbc', length=16, molecular_feature=feature('cell_barcode', 'cell_barcode', 'whitelist', whitelist='10x Chromium chemistry-specific whitelist')), unknown_region('chemistry-specific UMI', tag='umi', molecular_feature=feature('umi', 'umi', 'random')), token('variant-containing cDNA', length=44)))
primers = ()
reads = None
OLIGOS = (
    "<b>P5generic</b> (Supplementary Table 5) — 5′-AATGATACGGCGACCACCGAGATCTACAC-3′",
    "<b>RPI-X</b> (Supplementary Table 5) — 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′; XXXXXXXX is the sample index",
    "<b>Staggered GoT handles</b> (Supplementary Table 5) — 5′-CACCCGAGAATTCCA-, CACCCGAGAATTCCAA-, CACCCGAGAATTCCATT- or CACCCGAGAATTCCACAT-[gene-specific primer]-3′",
    "<b>CALR L367fs*46 gene-specific portion</b> — 5′-[staggered handle]-GCAGCAGAGAAACAAATGAAGGAC-3′",
    "<b>SI-PCR forward primer</b> — 10x Genomics reagent annealing partial Read 1; sequence not printed in the defining paper",
)
workflow_fn = got_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def got_graph(definition, final: Construct) -> ProtocolGraph:
    """GoT graph retaining both expression and targeted-genotype outputs."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    initial = got_workflow().initial
    cdna = graph.input("barcoded_10x_cdna", "amplified barcoded 10x cDNA",
                       MaterialKind.MOLECULAR, state=initial, primary=True)
    expression_cdna = graph.product("expression_cdna", "whole-transcriptome cDNA aliquot",
                                    MaterialKind.MOLECULAR,
                                    state=single("expression branch", (
                                        token("whole-transcriptome cDNA"),)))
    genotype_cdna = graph.product("genotype_cdna", "targeted-genotyping cDNA aliquot",
                                  MaterialKind.MOLECULAR,
                                  state=single("genotype branch", (
                                      token("target cDNA"),
                                      token("cell barcode", tag="cbc", length=16),
                                      unknown_region("chemistry-specific UMI", tag="umi"))))
    graph.react("split_cdna", "Split cDNA into expression and genotyping branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),),
                (expression_cdna, genotype_cdna))

    expression_forward = graph.input("expression_forward_primer",
                                     "10x expression-library forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "10x expression-library indexed reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    polymerase = graph.input("library_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    expression_final = custom_library("GoT matched expression library", (
        token("whole-transcriptome cDNA", length=44),))
    expression_library = graph.product(
        "expression_library", "matched 10x expression library",
        MaterialKind.LIBRARY, state=expression_final, inferred=True)
    graph.react("expression_library_pcr", "Complete the standard expression library",
                StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The exact primers depend on the selected 10x 3-prime or 5-prime chemistry."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The matched expression endpoint follows the chosen 10x chemistry."))

    target_forward = graph.input("target_forward_primer",
                                 "10x SI-PCR generic forward primer (partial Read 1)",
                                 MaterialKind.OLIGO)
    target_reverse = graph.input("target_reverse_primer",
                                 "staggered locus-specific reverse primer (partial Read 2 handle)",
                                 MaterialKind.OLIGO)
    amplicon = graph.product(
        "target_amplicon", "cell-barcode- and UMI-retaining GoT amplicon",
        MaterialKind.MOLECULAR,
        state=single("GoT amplicon", (
            token("variant-containing target"),
            token("cell barcode", tag="cbc", length=16),
            unknown_region("chemistry-specific UMI", tag="umi"))))
    graph.react("target_pcr", "Amplify from the variant toward the 10x end",
                StepKind.PCR,
                ((genotype_cdna, InputRole.TEMPLATE),
                 (target_forward, InputRole.STRUCTURAL),
                 (target_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplicon,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Target-specific primer sequence varies with the selected locus."))
    spri = graph.input("target_spri", "SPRI purification", MaterialKind.CONDITION)
    purified = graph.product(
        "purified_target_amplicon", "SPRI-purified GoT amplicon", MaterialKind.MOLECULAR,
        state=single("GoT amplicon", (
            token("variant-containing target"),
            token("cell barcode", tag="cbc", length=16),
            unknown_region("chemistry-specific UMI", tag="umi"))))
    primer_waste = graph.product("unincorporated_target_primers",
                                 "unincorporated target-PCR primers", MaterialKind.WASTE)
    graph.react("target_spri_cleanup", "SPRI purification to remove unincorporated primers",
                StepKind.CLEANUP, ((amplicon, InputRole.TEMPLATE), (spri, InputRole.CONDITION)),
                (purified, primer_waste), changes_molecule=False)
    graph.discard(primer_waste, "removed by SPRI purification before the second PCR")
    handle_forward = graph.input(
        "handle_forward_primer", "P5_generic forward primer", MaterialKind.OLIGO,
        state=single("P5generic", (Segment("P5generic", P5_GENERIC),)))
    handle_reverse = graph.input(
        "handle_reverse_primer", "RPI-x indexed reverse primer (Illumina)", MaterialKind.OLIGO,
        state=single("RPI-X", (
            Segment("P7", illumina.P7),
            Segment("sample index", "XXXXXXXX", "idx", placeholder=True,
                    feature=feature("sample_index", "sample_index", "unknown")),
            Segment("RPI-X 3′ handle", RPI_X_3PRIME))))
    genotype_library = graph.product("genotype_library", "GoT targeted genotype library",
                                     MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("handle_pcr", "P5_generic/RPI-x PCR completes P7 and adds the sample index", StepKind.PCR,
                ((purified, InputRole.TEMPLATE),
                 (handle_forward, InputRole.STRUCTURAL),
                 (handle_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (genotype_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The endpoint varies with the selected 10x chemistry and target."))
    graph.sequence(genotype_library, SequencingPlan(limitation=definition.endpoint))

    graph.declare_reviewed_scope((
        "barcoded_10x_cdna", "split_cdna", "expression_forward_primer",
        "expression_reverse_primer", "library_polymerase", "expression_library_pcr",
        "expression_library", "target_forward_primer", "target_reverse_primer",
        "target_pcr", "target_spri", "target_spri_cleanup", "purified_target_amplicon",
        "unincorporated_target_primers", "handle_forward_primer", "handle_reverse_primer", "handle_pcr",
        "genotype_library",
    ), source)
    optional = SourceItem("optional_rt_spikein", "optional target-primer spike-in", source)
    graph.add_source_item(optional)
    graph.cover_source_item(SourceCoverage(
        optional.id, "excluded",
        reason="Optional target-primer spike-ins are chemistry- and locus-dependent; this graph begins at amplified barcoded 10x cDNA."))
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'got'
    return got_graph(d, final)
SPEC.protocol_graph = protocol_graph
