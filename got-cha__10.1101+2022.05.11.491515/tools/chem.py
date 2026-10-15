from types import SimpleNamespace
from protocol_model import chain, role_library
from multimodal_protocols import CELL, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SAMPLE, StepKind, _finish_pcr, molecule, role, role_library
DEFINITION = Definition("GoT-ChA — genotyping of targeted loci with chromatin accessibility", "10.1101/2022.05.11.491515",
        "Namikawa et al., bioRxiv (2022)", "A locus-specific primer participates in single-cell ATAC droplet barcoding, after which hemi-nested biotin capture enriches the genomic variant without losing the ATAC cell barcode.",
        "The genotype read spans the chosen variant; the cell barcode links it to the parent accessibility profile.", "01_got_cha.html")


def _gotcha_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "single-cell ATAC reaction containing a target genomic locus",
                        MaterialKind.BIOLOGICAL, state=_gotcha().initial, primary=True)
    bead_primer = graph.input("bead_primer", "droplet cell-barcode extension primer",
                              MaterialKind.OLIGO)
    locus_primer = graph.input("locus_primer", "locus-specific spike-in primer",
                               MaterialKind.OLIGO)
    polymerase = graph.input("extension_polymerase", "droplet extension polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and target-locus products",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded ATAC and target molecules", (
                                 CELL, role("accessible-chromatin fragments"),
                                 role("target-locus amplicon"))))
    graph.react("barcode_extension", "Copy the droplet barcode onto ATAC and target molecules",
                StepKind.EXTENSION,
                ((cells, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (locus_primer, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (barcoded,))
    atac = graph.product("atac_branch", "parent scATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    target = graph.product("target_branch", "target-enrichment aliquot", MaterialKind.MOLECULAR,
                           state=barcoded.state)
    graph.react("split_branches", "Split parent ATAC and target-genotype branches",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, target),
                changes_molecule=False)
    reviewed = ["cells", "bead_primer", "locus_primer", "extension_polymerase",
                "barcode_extension", "split_branches"]
    atac_final = role_library("GoT-ChA parent ATAC roles",
                              (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final,
                "Parent scATAC branch linked by the same cell barcode.", reviewed,
                label="GoT-ChA parent ATAC")
    nested_forward = graph.input("nested_forward", "biotinylated hemi-nested target primer",
                                 MaterialKind.OLIGO)
    nested_reverse = graph.input("nested_reverse", "target PCR reverse primer",
                                 MaterialKind.OLIGO)
    nested_polymerase = graph.input("nested_polymerase", "target PCR polymerase",
                                    MaterialKind.CATALYST)
    enriched = graph.product("enriched_target", "biotinylated variant-spanning genotype amplicon",
                             MaterialKind.MOLECULAR,
                             state=molecule("enriched genotype amplicon", (
                                 CELL, role("variant-spanning genomic insert"),
                                 role("biotinylated nested end"))))
    graph.react("nested_pcr", "Perform biotinylated hemi-nested target PCR", StepKind.PCR,
                ((target, InputRole.TEMPLATE), (nested_forward, InputRole.STRUCTURAL),
                 (nested_reverse, InputRole.STRUCTURAL),
                 (nested_polymerase, InputRole.CATALYST)), (enriched,),
                geometry=PCRGeometry.unavailable(
                    enriched.state, reason="The locus-specific panel determines the target-primer geometry."))
    beads = graph.input("streptavidin", "streptavidin capture beads", MaterialKind.SELECTION)
    captured = graph.product("captured_target", "streptavidin-captured genotype amplicon",
                             MaterialKind.MOLECULAR, state=enriched.state)
    uncaptured = graph.product("uncaptured_target", "uncaptured target-PCR material",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated target amplicons on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((enriched, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not carry the selected biotinylated target amplicon")
    reviewed.extend(("nested_forward", "nested_reverse", "nested_polymerase",
                     "nested_pcr", "streptavidin", "capture"))
    _finish_pcr(graph, "genotype", captured, final, d.endpoint, reviewed,
                label="GoT-ChA genotype")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'got-cha'
    return _gotcha_graph(d, final)

def _gotcha():
    return chain(molecule("single-cell ATAC reaction", (role("tagmented accessible DNA"),
                                                         role("target genomic locus"))), (
        ("Spike a locus-specific primer into droplet barcoding extension",
         molecule("barcoded ATAC and target molecules", (CELL,
             role("ordinary accessible-chromatin fragments"), role("target-locus amplicon"))),
         "The target primer copies the same droplet cell barcode used by the scATAC library onto the genomic locus."),
        ("Split off an aliquot and perform biotinylated hemi-nested target PCR",
         molecule("enriched genotype amplicon", (CELL, role("variant-spanning genomic insert"),
                                                  role("biotinylated nested end"))),
         "Nested amplification enriches the selected locus while ordinary ATAC fragments are only linearly carried through."),
        ("Streptavidin-capture and index on-bead",
         molecule("GoT-ChA library", (role("P5 / Read 1 end"), CELL,
                                       role("variant-spanning insert"), SAMPLE,
                                       role("Read 2 / P7 end"))),
         "The genotype library can be pooled with the parent scATAC library and joined computationally by cell barcode."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (CELL, role("variant-spanning genomic insert")), index=SAMPLE)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT="",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Target panels, pool barcodes and kit revisions vary; the defining source establishes the read roles, while this page does not invent a universal primer sequence.",
    workflow=_gotcha,
    oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
