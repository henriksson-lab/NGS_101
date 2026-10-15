from types import SimpleNamespace
from multimodal_protocols import SAMPLE
from protocol_model import chain, role_library
from multimodal_protocols import ANTIBODY, CELL, Definition, InputRole, MaterialKind, POOL, ProtocolGraph, StepKind, UMI, _finish_pcr, molecule, role
DEFINITION = Definition("SCITO-seq — single-cell combinatorial indexed cytometry sequencing", "10.1038/s41592-021-01222-3",
        "Hwang et al., Nature Methods (2021)", "Hybridized splints add antibody and staining-pool barcodes before droplet feature capture adds a second cell barcode and UMI.",
        "Protein identity and staining-pool identity are read with the droplet barcode and UMI.", "01_scito_seq.html")


def _scito_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cell_pools", "cells divided into staining pools",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("cell pools", (role("surface proteins"),)), primary=True)
    antibodies = graph.input("antibodies", "universally oligo-conjugated antibodies",
                             MaterialKind.MOLECULAR)
    splints = graph.input("splints", "antibody- and pool-identifying splint oligos",
                          MaterialKind.OLIGO)
    tags = graph.product("antibody_tags", "SCITO antibody tags carrying target and pool identity",
                         MaterialKind.SELECTION,
                         state=molecule("SCITO antibody tag", (ANTIBODY, POOL,
                             role("10x feature-capture sequence"))))
    graph.react("hybridize_splints", "Hybridize pool-specific splints to universal antibody handles",
                StepKind.HYBRIDIZATION,
                ((antibodies, InputRole.TEMPLATE), (splints, InputRole.STRUCTURAL)), (tags,))
    stained = graph.product("stained_cells", "SCITO-tagged cells", MaterialKind.MOLECULAR,
                            state=molecule("tagged cell", (ANTIBODY, POOL,
                                role("surface-protein counts"))))
    unbound = graph.product("unbound_tags", "unbound antibody tags", MaterialKind.WASTE)
    graph.react("stain", "Stain cells with target- and pool-barcoded antibodies",
                StepKind.BINDING,
                ((cells, InputRole.TEMPLATE), (tags, InputRole.SELECTION)),
                (stained, unbound))
    graph.discard(unbound, "does not report a cell-associated protein measurement")
    droplets = graph.input("droplets", "10x feature-barcoding droplet condition",
                           MaterialKind.CONDITION)
    occupied = graph.product("occupied_droplets", "overloaded droplets containing pool-resolved cells",
                             MaterialKind.MOLECULAR, state=stained.state)
    empty = graph.product("empty_droplets", "empty or unresolved droplets", MaterialKind.WASTE)
    graph.react("partition", "Combine stained pools and partition cells into droplets",
                StepKind.PARTITION,
                ((stained, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (occupied, empty), changes_molecule=False)
    graph.discard(empty, "does not contain a resolvable pool/cell combination")
    bead_primer = graph.input("bead_primer", "10x cell-barcode/UMI feature-capture primer",
                              MaterialKind.OLIGO)
    polymerase = graph.input("extension_polymerase", "feature-tag extension polymerase",
                             MaterialKind.CATALYST)
    captured = graph.product("captured_tags", "droplet-barcoded SCITO feature tags",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured SCITO tags", (CELL, UMI,
                                 ANTIBODY, POOL)))
    graph.react("capture_tags", "Copy droplet barcode and UMI onto antibody tags",
                StepKind.EXTENSION,
                ((occupied, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (captured,))
    reviewed = ["cell_pools", "antibodies", "splints", "hybridize_splints", "stain",
                "droplets", "partition", "bead_primer", "extension_polymerase",
                "capture_tags"]
    _finish_pcr(graph, "feature", captured, final, d.endpoint, reviewed,
                label="SCITO feature-tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'scito-seq'
    return _scito_graph(d, final)

def _scito():
    return chain(molecule("universally oligo-conjugated antibody", (
        role("antibody"), role("20-nt antibody handle"))), (
        ("Hybridize a pool-specific splint carrying antibody and pool identity",
         molecule("SCITO antibody tag", (role("antibody handle complement"),
                                          role("TruSeq Read 2 handle"), ANTIBODY,
                                          POOL, role("10x feature-capture sequence"))),
         "One antibody conjugation can be reused; the splint supplies both target and staining-pool identity."),
        ("Stain split cell pools, combine them and overload droplets",
         molecule("tagged cell", (ANTIBODY, POOL, role("surface-protein counts"))),
         "Pool identity makes many otherwise colliding cells distinguishable."),
        ("Capture tags on gel beads and copy droplet barcode plus UMI",
         molecule("SCITO-seq library", (role("P5 end"), CELL, UMI,
                                         ANTIBODY, POOL, role("Read 2 / P7 end"))),
         "Antibody-plus-pool barcode and droplet barcode form the combinatorial cell identity."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (CELL, UMI, ANTIBODY, POOL), index=SAMPLE)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT="The defining paper fixes the splint architecture and lengths, "
                       "including a 20-nt antibody handle, TruSeq Read 2 region, 10-nt "
                       "compound antibody/pool barcode and 22-nt 10x feature-capture "
                       "region. It does not print the handle or barcode bases, and the "
                       "completed library follows the commercial 10x 3′ v3 Feature "
                       "Barcoding workflow.",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="The original article, author manuscript, "
                                       "Supplementary Tables 1–4 and authors' public "
                                       "analysis repository do not disclose a universal "
                                       "SCITO handle/barcode set. Read geometry also "
                                       "inherits the commercial 10x 3′ v3 endpoint.",
    workflow=_scito,
    oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
