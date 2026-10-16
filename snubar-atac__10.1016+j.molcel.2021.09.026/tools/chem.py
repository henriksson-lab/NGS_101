"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from chemdraw import Segment, feature, oligo
from batch_ngs import truseq_library
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_specs import role
from rna_library_operations import finish_pcr
from rna_library_operations import molecule

DEFINITION = Definition("SNuBar-ATAC — oligo sample multiplexing for scATAC", "10.1016/j.molcel.2021.09.026", "Wang et al., Molecular Cell (2021)", "An Illumina TDE1 transposome tags each sample, then a sample-barcode oligo anneals to its universal tail before nuclei are pooled for 10x scATAC.", "nextera", ("cell barcode", "accessible DNA"), "After 10x GEM cleanup, separate ATAC and SBO barcode templates are PCR-amplified and the two libraries are mixed for sequencing. ATAC reads report genomic accessibility and a cell barcode; the SBO branch reports sample barcode and UMI.", "DNA", "The displayed final molecule is the ATAC branch. Original STAR Methods print the generic SBO scaffold and barcode-library PCR primers; that separate SBO branch appears in the reaction graph. Table S1 with individual 18-base sample barcodes remains inaccessible.")


def snubar_oligos():
    """Published ATAC-only SBO family and barcode-library primers (STAR Methods)."""
    return (
        oligo("SNuBar-ATAC SBO family", [
            Segment("PCR handle", "CCTTGGCACCCGAGAATTCCA"),
            Segment("UMI", "N" * 8, "umi", placeholder=True, length_bp=8,
                    feature=feature("snubar_umi", "umi", "random")),
            Segment("sample barcode", "N" * 18, "cbc", placeholder=True, length_bp=18,
                    feature=feature("snubar_sample", "sample_index", "unknown")),
            Segment("CST transposome-tail complement", "CTACTCTGCGTTGATACCAGACGCTGCCGACGA"),
        ]),
        oligo("TruSeq RPIX barcode-library primer", [
            Segment("P7 and index-primer arm", "CAAGCAGAAGACGGCATACGAGAT"),
            Segment("sample index", "N" * 6, "cbc", placeholder=True, length_bp=6,
                    feature=feature("snubar_pcr_index", "sample_index", "unknown")),
            Segment("barcode-library arm", "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"),
        ]),
        oligo("P5 barcode-library primer", [
            Segment("P5 primer", "AATGATACGGCGACCACCGAGATCTACAC"),
        ]),
    )


from protocol_specs import workflow, two_branch


def _source_workflow():
    return workflow(molecule("sample-specific nuclei", ("chromatin",)), (
        ("Tagment each sample with Illumina TDE1",
         molecule("tagmented chromatin", ("genomic fragment", "transposome universal tail")),
         "The sample barcode is not carried by Tn5 itself."),
        ("Anneal a SNuBar oligo to the transposome tail",
         molecule("sample-tagged fragment", ("genomic fragment", "sample barcode", "UMI", "PCR handle"), joined=True),
         "The SBO transfers sample identity before samples are pooled."),
        ("Pool samples, run 10x GEM barcoding, then separate ATAC and SBO barcode templates",
         two_branch("same cell", "ATAC library", "SBO barcode library"),
         "The two products are separately amplified, then mixed for sequencing."),
    ))

SPEC = get('snubar-atac', definition=DEFINITION,
           workflow_fn=_source_workflow, oligos_fn=snubar_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "sample-separated nuclei containing accessible chromatin",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("sample-specific nuclei", ("accessible chromatin",)),
                         primary=True)
    transposome = graph.input("transposome", "Illumina TDE1 loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented chromatin carrying the universal tail",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagmented chromatin", (
                               "genomic fragment", "transposome universal tail")))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment each sample with Illumina TDE1",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an accessible-DNA insert with the universal tail")
    snubar = graph.input("snubar", "sample-barcode/UMI SNuBar oligo",
                         MaterialKind.OLIGO)
    sample_tagged = graph.product("sample_tagged", "SNuBar-tagged ATAC fragments",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("sample-tagged fragment", (
                                      "genomic fragment", "sample barcode", "UMI", "PCR handle")))
    graph.react("anneal_snubar", "Anneal the SNuBar oligo to the transposome tail",
                StepKind.HYBRIDIZATION,
                ((tagged, InputRole.TEMPLATE), (snubar, InputRole.STRUCTURAL)),
                (sample_tagged,))
    droplets = graph.input("droplets", "10x single-cell ATAC droplet condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "cell- and sample-barcoded ATAC molecules",
                             MaterialKind.MOLECULAR,
                             state=molecule("SNuBar-ATAC products", (
                                 "10x cell barcode", "sample barcode", "UMI", "ATAC fragment")))
    collisions = graph.product("collisions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("partition", "Pool samples and partition nuclei into 10x droplets",
                StepKind.PARTITION,
                ((sample_tagged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, collisions))
    graph.discard(collisions, "does not retain an interpretable cell barcode")
    atac_template = graph.product("atac_template", "cell-barcoded ATAC template",
                                  MaterialKind.MOLECULAR,
                                  state=molecule("ATAC template", ("10x cell barcode", "ATAC fragment")))
    barcode_template = graph.product("barcode_template", "cell-barcoded SBO sample-tag template",
                                     MaterialKind.MOLECULAR,
                                     state=molecule("SBO barcode template", (
                                         "10x cell barcode", "sample barcode", "UMI")))
    graph.react("separate_templates", "Separate ATAC and SBO barcode templates after GEM cleanup",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),),
                (atac_template, barcode_template), changes_molecule=False)
    reviewed = ["nuclei", "transposome", "tn5", "tagment", "snubar", "anneal_snubar",
                "droplets", "partition", "separate_templates"]
    finish_pcr(graph, "atac", atac_template, final, primers, d.sequencing, reviewed,
               label="SNuBar-ATAC accessible-DNA")
    barcode_final, barcode_primers = truseq_library([
        role("cell barcode", molecular_feature=feature("snubar_cell", "cell_barcode", "unknown")),
        role("sample barcode", length_bp=18,
             molecular_feature=feature("snubar_sample", "sample_index", "unknown")),
        role("UMI", length_bp=8,
             molecular_feature=feature("snubar_umi", "umi", "random"))],
        "SNuBar-ATAC SBO barcode library")
    finish_pcr(graph, "sbo", barcode_template, barcode_final, barcode_primers,
               "separate SBO sample-barcode read", reviewed, label="SNuBar-ATAC SBO barcode")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
