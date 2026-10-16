"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from batch_ngs import nextera_library
from batch_ngs import truseq_library
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_graph import SourceClaim, SourceDiscrepancy
from protocol_graph import SequencingPlan
from protocol_specs import role
from rna_library_operations import finish_pcr
from rna_library_operations import molecule
from pcr import PCRGeometry

DEFINITION = Definition("ReDeeM — deep mitochondrial mutation profiling", "10.1038/s41586-024-07066-z", "Weng et al., Nature (2024)", "A whole-cell 10x Multiome modification retains mtDNA, then tiled hybrid capture produces a deeply sequenced mtDNA branch matched to RNA and ATAC.", "truseq", ("cell barcode", "mtDNA fragment"), "The mtDNA, ATAC and RNA libraries carry matched cell barcodes; mtDNA fragment endpoints support consensus error correction.", "multi", "The published supplements disagree on mitochondrial capture-probe lengths and cite both Rev A and Rev B of 10x CG000338. The model uses Supplementary Data 1 for probe bases; both available 10x guide revisions publish the same inspected fixed RNA/ATAC bead and library shells. Branch-specific ReDeeM junctions remain to be reconstructed.")


from protocol_specs import workflow, two_branch


def _source_workflow():
    return workflow(two_branch("fixed permeabilized whole cell", "nuclear chromatin", "retained mtDNA"), (
        ("Run 10x Multiome and amplify its ATAC branch",
         two_branch("matched cell barcode", "RNA library", "ATAC library carrying mtDNA"),
         "The mtDNA capture input is an aliquot of the indexed ATAC library."),
        ("Hybrid-capture mtDNA from the ATAC library in four staggered probe pools",
         molecule("captured mtDNA branch", ("cell barcode", "mtDNA fragment", "capture probe")),
         "Biotinylated probes and streptavidin beads select mtDNA; fragment endpoints act as endogenous molecular identifiers."),
        ("Sequence mtDNA, ATAC and RNA as three matched libraries",
         two_branch("same cell barcode", "mtDNA library", "ATAC + RNA libraries"),
         "The mtDNA branch receives deep error-corrected sequencing."),
    ))

SPEC = get('redeem', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "whole cells retaining nuclear chromatin, RNA and mtDNA",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("whole cell", ("nuclear chromatin", "RNA", "mtDNA")),
                        primary=True)
    droplets = graph.input("multiome_droplets", "10x Multiome droplet/barcoding condition",
                           MaterialKind.CONDITION)
    barcoded = graph.product("barcoded", "matched cell-barcoded RNA and ATAC pre-amplification products",
                             MaterialKind.MOLECULAR,
                             state=molecule("matched cell products", (
                                 "cell barcode", "RNA products", "ATAC fragments including mtDNA")))
    failed = graph.product("failed_partitions", "empty or ambiguous droplets", MaterialKind.WASTE)
    graph.react("multiome", "Run 10x Multiome while retaining mitochondrial DNA",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable shared cell barcode")
    atac = graph.product("atac_branch", "matched ATAC branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    rna = graph.product("rna_branch", "matched RNA branch", MaterialKind.MOLECULAR,
                        state=barcoded.state)
    graph.react("split_branches", "Separate ATAC and RNA branches",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, rna),
                changes_molecule=False)
    atac_final, atac_primers = nextera_library(
        [role("cell barcode"), role("accessible-DNA insert")], "ReDeeM ATAC library")
    rna_final, rna_primers = truseq_library(
        [role("cell barcode"), role("UMI"), role("transcript cDNA")], "ReDeeM RNA library")
    reviewed = ["cells", "multiome_droplets", "multiome", "split_branches",
                "capture_probes", "capture_mt"]
    atac_forward = graph.input("atac_forward", "ATAC indexing forward primer", MaterialKind.OLIGO)
    atac_reverse = graph.input("atac_reverse", "ATAC indexing reverse primer", MaterialKind.OLIGO)
    atac_polymerase = graph.input("atac_polymerase", "ATAC indexing PCR polymerase", MaterialKind.CATALYST)
    indexed_atac = graph.product("atac_indexed", "indexed ATAC library including mtDNA fragments",
                                 MaterialKind.LIBRARY, state=atac_final)
    graph.react("atac_pcr", "Amplify and index the ATAC library", StepKind.PCR,
                ((atac, InputRole.TEMPLATE), (atac_forward, InputRole.STRUCTURAL),
                 (atac_reverse, InputRole.STRUCTURAL), (atac_polymerase, InputRole.CATALYST)),
                (indexed_atac,), geometry=PCRGeometry.unavailable(
                    atac_final, reason="Exact ReDeeM ATAC indexing primer-to-product geometry has not been reconstructed."))
    atac_for_sequencing = graph.product("atac_library", "ATAC library aliquot for sequencing",
                                        MaterialKind.LIBRARY, state=atac_final)
    atac_for_capture = graph.product("atac_capture_input", "ATAC library aliquot for mtDNA capture",
                                     MaterialKind.LIBRARY, state=atac_final)
    graph.react("split_atac_library", "Reserve indexed ATAC aliquots for sequencing and mtDNA capture",
                StepKind.SPLIT, ((indexed_atac, InputRole.TEMPLATE),),
                (atac_for_sequencing, atac_for_capture), changes_molecule=False)
    graph.sequence(atac_for_sequencing, SequencingPlan(library=atac_final, primers=tuple(atac_primers)))
    reviewed.extend(("atac_forward", "atac_reverse", "atac_polymerase", "atac_pcr",
                     "atac_library", "split_atac_library"))
    finish_pcr(graph, "rna", rna, rna_final, rna_primers,
               "matched RNA branch: " + d.sequencing, reviewed, label="ReDeeM RNA")
    probes = graph.input("capture_probes", "four staggered 5′-biotinylated mtDNA capture-probe pools (bases in Supplementary Data 1)",
                         MaterialKind.SELECTION)
    captured = graph.product("captured_mt", "hybrid-captured cell-barcoded mtDNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured mtDNA", ("cell barcode", "mtDNA fragment")))
    rejected = graph.product("rejected_mt", "non-mitochondrial capture material",
                             MaterialKind.WASTE)
    graph.react("capture_mt", "Hybrid-capture mtDNA from indexed ATAC library in four staggered probe pools",
                StepKind.AFFINITY_CAPTURE,
                ((atac_for_capture, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not belong to the selected mitochondrial branch")
    finish_pcr(graph, "mt", captured, final, primers, "mtDNA branch: " + d.sequencing,
                reviewed, label="ReDeeM mtDNA")
    graph.declare_reviewed_scope(reviewed, source)
    graph.annotate(SourceDiscrepancy(
        "capture_probes", "ReDeeM capture-probe specification changed across published supplements",
        (SourceClaim("Weng et al., Nature (2024), supplementary ReDeeM protocol, extracted lines 64–65",
                     "Four staggered sets of 137 probes are described as 100 nt, with sequences forthcoming."),
         SourceClaim("Same article, Supplementary Methods, extracted lines 456–465",
                     "Each of four sets is described as 138 biotinylated probes of 120 nt."),
         SourceClaim("Same article, Supplementary Data 1, eight MitoProbe sheets",
                     "Individual human and mouse probe bases are supplied, usually 120 nt; sheet counts vary and a few human v1 rows are shorter.")),
        "The model treats Supplementary Data 1 as the sequence inventory for the published version.",
        "The protocol text predates the workbook and disagrees with the final supplementary Methods and sequence inventory.",
        "The workbook lists bases but does not by itself establish each oligo's synthesis modification or a sample-specific pool assignment."))
    graph.annotate(SourceDiscrepancy(
        "multiome", "ReDeeM supplements cite different revisions of the 10x Multiome guide",
        (SourceClaim("Weng et al., Nature (2024), Supplementary Methods, extracted lines 433–440",
                     "GEM generation and RNA library preparation cite CG000338 Rev A."),
         SourceClaim("Same article, supplementary ReDeeM protocol, extracted lines 144–156 and 416",
                     "The corresponding steps cite CG000338 Rev B.")),
        "The model uses the fixed 10x oligo shells printed identically in the available Rev A and Rev B Sequences appendices.",
        "The supplements do not say which revision was used in each experiment.",
        "The inspected sequence appendices agree, but the ReDeeM-specific junctions and any procedural differences between guide revisions remain unresolved."))
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
