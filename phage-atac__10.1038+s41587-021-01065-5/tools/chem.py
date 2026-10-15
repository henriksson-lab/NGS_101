"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, molecule, role, role_library
from extended_protocols import (Definition, CELL, UMI, SAMPLE, TARGET, GUIDE, RT_WELL, LIG_WELL, PCR_WELL, _finish_pcr, _rt_then_pcr)

DEFINITION = Definition('PHAGE-ATAC', '10.1038/s41587-021-01065-5', 'Fiskin et al., Nature Biotechnology (2022)', 'Nanobody-displaying phagemids carry 10x-compatible PAC tags so one droplet barcode joins protein-derived tags to ATAC and mitochondrial DNA.', 'Separate protein-tag and ATAC libraries share the 10x cell barcode.', '01_phage_atac.html')
PAYLOAD = (CELL, TARGET, role("ATAC or phage insert"))

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ""
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=SAMPLE)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The source defines the read roles, but the primer sequence is panel-, kit- or implementation-dependent; no universal binding oligo is asserted."

def workflow():
    return chain(molecule("cell stained with nanobody-displaying M13 phage", (role("surface or intracellular epitope"), TARGET)), (
        ("Fix, lyse and tagment accessible chromatin in bulk",
         molecule("phage-bound tagmented cell", (TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))),
         "The engineered PAC/RD1 sequence makes the phagemid compatible with 10x ATAC capture."),
        ("Co-encapsulate cells and 10x ATAC gel beads",
         molecule("droplet-barcoded products", (CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))),
         "Linear amplification applies one droplet barcode to chromatin and phage-derived tags."),
        ("Split protein-derived-tag and ATAC libraries",
         molecule("PHAGE-ATAC readouts", (CELL, TARGET, role("ATAC or phage insert"), SAMPLE)),
         "The common cell barcode links protein abundance, accessibility and captured mitochondrial variants."),
    ))


def oligos():
    return ()


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input(
        "stained_cells", "cells stained with nanobody-displaying M13 phage",
        MaterialKind.BIOLOGICAL,
        state=molecule("phage-stained cell", (
            role("surface or intracellular epitope"), TARGET, role("PAC-tag phagemid"))),
        primary=True)
    transposome = graph.input("atac_transposome", "loaded ATAC transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "ATAC Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product(
        "tagmented_cells", "phage-bound cells with tagmented accessible chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("phage-bound tagmented cell", (
            TARGET, role("PAC-tag phagemid"), role("ATAC fragments"))))
    losses = graph.product("tn5_losses", "non-library ATAC products", MaterialKind.WASTE)
    graph.react("tagment", "Fix, lyse and tagment accessible chromatin",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, losses))
    graph.discard(losses, "does not retain an amplifiable accessible-chromatin fragment")
    droplets = graph.input("droplet_condition", "10x ATAC gel-bead co-encapsulation condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "cells co-encapsulated with ATAC gel beads",
                             MaterialKind.MOLECULAR, state=tagmented.state)
    empty = graph.product("empty_droplets", "empty or collision droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and ATAC gel beads",
                StepKind.PARTITION,
                ((tagmented, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty), changes_molecule=False)
    graph.discard(empty, "does not contain one interpretable cell/bead combination")
    bead_primer = graph.input("bead_primer", "10x ATAC cell-barcode bead primer",
                              MaterialKind.OLIGO)
    polymerase = graph.input("linear_polymerase", "droplet linear-amplification polymerase",
                             MaterialKind.CATALYST)
    barcoded = graph.product(
        "barcoded_products", "cell-barcoded ATAC and phagemid products",
        MaterialKind.MOLECULAR,
        state=molecule("droplet-barcoded products", (
            CELL, TARGET, role("phagemid tag"), role("ATAC fragment"))))
    graph.react("linear_amplification", "Apply the droplet barcode to chromatin and phage tags",
                StepKind.EXTENSION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (barcoded,))
    atac = graph.product("atac_branch", "cell-barcoded ATAC products",
                         MaterialKind.MOLECULAR, state=barcoded.state)
    protein = graph.product("protein_branch", "cell-barcoded protein-derived phage tags",
                            MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("split_readouts", "Split ATAC and protein-derived-tag libraries",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, protein),
                changes_molecule=False)
    reviewed = [
        "stained_cells", "atac_transposome", "tn5", "tagment",
        "droplet_condition", "partition", "bead_primer", "linear_polymerase",
        "linear_amplification", "split_readouts",
    ]
    _finish_pcr(graph, "atac", atac, final, "ATAC branch: " + d.endpoint,
                reviewed, label="PHAGE-ATAC chromatin")
    _finish_pcr(graph, "protein", protein, final, "protein-tag branch: " + d.endpoint,
                reviewed, label="PHAGE-ATAC protein tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
