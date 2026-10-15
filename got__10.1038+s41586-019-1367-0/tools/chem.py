from chemdraw import Workflow
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
        ("Target-PCR from upstream of the variant toward the 10x end",
         single("GoT amplicon", (token("variant-containing target"),
                                  token("cell barcode", tag="cbc", length=16),
                                  unknown_region("chemistry-specific UMI", tag="umi",
                                      molecular_feature=feature("umi", "umi", "random")))),
         "The amplicon deliberately retains the cell barcode and UMI."),
        ("Add sequencing handles and sequence with the expression library",
         single("GoT library", (token("sequencing arm"), token("variant-containing target"),
                                 token("cell barcode", tag="cbc", length=16),
                                 unknown_region("chemistry-specific UMI", tag="umi",
                                     molecular_feature=feature("umi", "umi", "random")),
                                 token("indexed sequencing arm"))),
         "Read 1 reports cellular identifiers and Read 2 reaches the selected variant."),
    ))


DEFINITION = Definition("GoT — Genotyping of Transcriptomes", "10.1038/s41586-019-1367-0",
        "Nam et al., <i>Nature</i> (2019)", "01_got.html",
        "A targeted amplicon branch is made from 10x cDNA while retaining its cell barcode and UMI, linking a selected transcribed variant to the cell's whole transcriptome.",
        "For 3′ v3, Read 1 uses 28 cycles for cell barcode and UMI, the 8-cycle index selects the library, and the 98- or 130-cycle Read 2 reaches the targeted variant.",
        "The defining paper supports 10x 3′ v2/v3 and 5′ inputs. The schematic shows the common identifier-preserving branch; exact UMI length follows the chosen 10x chemistry.")

d = DEFINITION
final = custom_library('GoT targeted library', (token('cell barcode', tag='cbc', length=16, molecular_feature=feature('cell_barcode', 'cell_barcode', 'whitelist', whitelist='10x Chromium chemistry-specific whitelist')), unknown_region('chemistry-specific UMI', tag='umi', molecular_feature=feature('umi', 'umi', 'random')), token('variant-containing cDNA', length=44)))
primers = ()
reads = None
OLIGOS = ()
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

    target_forward = graph.input("target_forward_primer", "variant-upstream target primer",
                                 MaterialKind.OLIGO)
    target_reverse = graph.input("target_reverse_primer", "10x-end target primer",
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
    handle_forward = graph.input("handle_forward_primer", "GoT handle PCR forward primer",
                                 MaterialKind.OLIGO, inferred=True)
    handle_reverse = graph.input("handle_reverse_primer", "GoT indexed handle PCR reverse primer",
                                 MaterialKind.OLIGO, inferred=True)
    genotype_library = graph.product("genotype_library", "GoT targeted genotype library",
                                     MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("handle_pcr", "Add sequencing handles and sample index", StepKind.PCR,
                ((amplicon, InputRole.TEMPLATE),
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
        "target_pcr", "handle_forward_primer", "handle_reverse_primer", "handle_pcr",
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
