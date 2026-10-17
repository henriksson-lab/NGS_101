from types import SimpleNamespace
from multimodal_protocols import SAMPLE
from protocol_model import chain, role_library
from multimodal_protocols import ANTIBODY, CELL, Definition, InputRole, MaterialKind, POOL, ProtocolGraph, StepKind, UMI, molecule, role
from pcr import PCRGeometry
from protocol_graph import SequencingPlan
import seqprimers as sp
DEFINITION = Definition("SCITO-seq — single-cell combinatorial indexed cytometry sequencing", "10.1038/s41592-021-01222-3",
        "Hwang et al., Nature Methods (2021)", "Hybridized splints add antibody and staining-pool barcodes before droplet feature capture adds a second cell barcode and UMI.",
        "Protein identity and staining-pool identity are read with the droplet barcode and UMI.", "01_scito_seq.html")


def _scito_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cell_pools", "cells divided into staining pools",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("cell pools", (role("surface proteins"),)), primary=True)
    # Methods "Primary antibody oligonucleotide conjugation": each clone is
    # ThunderLink-conjugated (1 antibody : 3 oligos) to its own 20-bp 3'-amine,
    # HPLC-purified handle oligo; handle bases are not printed.
    raw_antibodies = graph.input("antibodies", "unconjugated antibody clones",
                                 MaterialKind.MOLECULAR)
    handles = graph.input("handle_oligos", "antibody-specific 20-bp 3'-amine handle oligos",
                          MaterialKind.OLIGO)
    conjugation = graph.input("thunderlink", "ThunderLink antibody-oligo conjugation kit",
                              MaterialKind.CATALYST)
    antibodies = graph.product("handle_conjugates", "antibodies carrying one clone-specific 20-bp handle",
                               MaterialKind.MOLECULAR,
                               state=molecule("handle-conjugated antibody",
                                              (role("antibody"), role("20-nt antibody handle"))))
    graph.react("conjugate_handles", "Conjugate each antibody clone to its 3'-amine handle oligo",
                StepKind.CHEMICAL_CONVERSION,
                ((raw_antibodies, InputRole.TEMPLATE), (handles, InputRole.STRUCTURAL),
                 (conjugation, InputRole.CATALYST)), (antibodies,))
    splints = graph.input("splints", "antibody- and pool-identifying splint oligos",
                          MaterialKind.OLIGO)
    tag_state = molecule("SCITO antibody tag", (ANTIBODY, POOL,
                                                role("10x feature-capture sequence")))
    hybrids = graph.product("hybridized_tags", "hybridized antibody tags with excess free splints",
                            MaterialKind.MOLECULAR, state=tag_state)
    # Methods "Antibody-specific hybridization design": 1 ug of each antibody
    # pooled, 1 ul of 1 uM secondary oligo, room temperature 15 min, then an
    # Amicon 50K MWCO column removes excess free oligonucleotides.
    graph.react("hybridize_splints", "Hybridize pool-specific splints to antibody handles (RT, 15 min)",
                StepKind.HYBRIDIZATION,
                ((antibodies, InputRole.TEMPLATE), (splints, InputRole.STRUCTURAL)), (hybrids,))
    tags = graph.product("antibody_tags", "purified SCITO antibody tags carrying target and pool identity",
                         MaterialKind.SELECTION, state=tag_state)
    free_splints = graph.product("free_splints", "excess free splint oligos", MaterialKind.WASTE)
    amicon = graph.input("amicon_column", "Amicon 50K MWCO centrifugal filter",
                         MaterialKind.CONDITION)
    graph.react("purify_tags", "Remove excess free splints on an Amicon 50K MWCO column",
                StepKind.CLEANUP,
                ((hybrids, InputRole.TEMPLATE), (amicon, InputRole.CONDITION)),
                (tags, free_splints), changes_molecule=False)
    graph.discard(free_splints, "filtrate removed by the 50K MWCO column before staining")
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
    # Methods "10x Genomics run for SCITO-seq": 3' v3 Feature Barcoding for Cell
    # Surface Proteins per manufacturer; the index-PCR primers and cycling live in
    # that guide (CG000185 Rev D cited for the TS-C variant), not in the paper.
    forward = graph.input("feature_forward", "kit P5-side feature PCR primer", MaterialKind.OLIGO)
    reverse = graph.input("feature_reverse", "kit i7 sample-index PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("feature_polymerase", "kit index-PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("feature_library", "SCITO feature-tag sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("feature_pcr", "Index PCR per the 10x 3' v3 Feature Barcoding guide",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The paper delegates index PCR to the 10x 3' v3 Feature "
                    "Barcoding guide; the applicable guide revision was not inspected."))
    graph.sequence(library, SequencingPlan(limitation=d.endpoint))
    reviewed = ["cell_pools", "antibodies", "handle_oligos", "thunderlink",
                "conjugate_handles", "splints", "hybridize_splints", "amicon_column",
                "purify_tags", "stain", "droplets", "partition", "bead_primer",
                "extension_polymerase", "capture_tags", "feature_forward",
                "feature_reverse", "feature_polymerase", "feature_pcr", "feature_library"]
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
         "One conjugation per clone is reused; the splint supplies both target and staining-pool identity. "
         "Excess free splints are removed on an Amicon 50K MWCO column before staining."),
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
    READ_SPECS=(
        sp.read_spec("Read 1", 28, "droplet barcode and UMI from the 10x 3′ v3 gel bead",
                     "extract"),
        sp.read_spec("Index 1 (i7)", 8, "10x sample index", "demultiplex"),
        sp.read_spec("Read 2", 98, "antibody-plus-pool barcode after the TruSeq Read 2 site; "
                     "the paper notes this read can be shortened", "extract"),
    ),
    SEQUENCING_ENDING=(d.endpoint + " Methods (10x Genomics run for SCITO-seq): NovaSeq 6000, "
                       "Read 1 28, i7 index 8 and Read 2 98 cycles."),
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
