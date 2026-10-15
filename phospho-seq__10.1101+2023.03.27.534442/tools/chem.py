from types import SimpleNamespace
from protocol_model import chain, role_library
from multimodal_protocols import ANTIBODY, CELL, Definition, InputRole, MaterialKind, ProtocolGraph, SAMPLE, StepKind, UMI, _finish_pcr, ident, molecule, role, role_library
DEFINITION = Definition("Phospho-seq — intracellular protein dynamics with chromatin accessibility", "10.1101/2023.03.27.534442",
        "Chung et al., bioRxiv (2023; published 2025)", "DNA-tagged intracellular antibodies are bridged into a 10x scATAC capture reaction; tag and ATAC branches are separated after GEM barcoding while retaining cell identity.",
        "The ATAC library and antibody/hash tag libraries share a droplet cell barcode; antibody barcodes identify total or phosphorylated targets.", "01_phospho_seq.html")


def _phospho_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "fixed permeabilized hashed cells with accessible chromatin",
                        MaterialKind.BIOLOGICAL, state=_phospho().initial, primary=True)
    antibodies = graph.input("antibodies", "DNA-tagged intracellular antibodies and hash tags",
                             MaterialKind.SELECTION)
    stained = graph.product("stained_cells", "antibody-tagged hashed cells",
                            MaterialKind.MOLECULAR,
                            state=molecule("tagged cell", (ANTIBODY,
                                ident("hash tag", "sample_hash", "sample_index"),
                                role("accessible chromatin"))))
    unbound = graph.product("unbound_antibodies", "unbound antibody tags", MaterialKind.WASTE)
    graph.react("stain", "Stain intracellular phosphoprotein targets",
                StepKind.BINDING,
                ((cells, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION)),
                (stained, unbound))
    graph.discard(unbound, "does not report a cell-associated intracellular target")
    transposome = graph.input("transposome", "10x scATAC loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product("tagmented_cells", "tagged cells with tagmented accessible chromatin",
                              MaterialKind.MOLECULAR, state=stained.state)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment accessible chromatin", StepKind.TAGMENTATION,
                ((stained, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, tn5_waste))
    graph.discard(tn5_waste, "does not retain an amplifiable accessibility fragment")
    bridge = graph.input("bridge_oligos", "antibody-tag-to-10x-capture bridge oligos",
                         MaterialKind.OLIGO)
    bridged = graph.product("bridged_cells", "cells with capturable ATAC and antibody-tag branches",
                            MaterialKind.MOLECULAR, state=_phospho().reactions[1].after)
    graph.react("bridge_tags", "Hybridize antibody tags to 10x-compatible bridge oligos",
                StepKind.HYBRIDIZATION,
                ((tagmented, InputRole.TEMPLATE), (bridge, InputRole.STRUCTURAL)), (bridged,))
    droplets = graph.input("droplets", "10x ATAC GEM barcoding condition", MaterialKind.CONDITION)
    barcoded = graph.product("barcoded_products", "cell-barcoded ATAC and antibody-tag products",
                             MaterialKind.MOLECULAR,
                             state=molecule("barcoded Phospho-seq products", (
                                 CELL, role("ATAC insert"), UMI, ANTIBODY,
                                 role("phosphoprotein tag"))))
    failed = graph.product("failed_partitions", "empty or ambiguous GEM products", MaterialKind.WASTE)
    graph.react("partition", "Partition and barcode both molecular branches in GEMs",
                StepKind.PARTITION,
                ((bridged, InputRole.TEMPLATE), (droplets, InputRole.CONDITION)),
                (barcoded, failed), changes_molecule=False)
    graph.discard(failed, "does not retain an interpretable cell barcode")
    atac = graph.product("atac_branch", "ATAC bead-bound branch", MaterialKind.MOLECULAR,
                         state=barcoded.state)
    tags = graph.product("tag_branch", "antibody/hash-tag supernatant branch",
                         MaterialKind.MOLECULAR, state=barcoded.state)
    graph.react("split_branches", "Separate ATAC DNA from antibody/hash-tag supernatant",
                StepKind.SPLIT, ((barcoded, InputRole.TEMPLATE),), (atac, tags),
                changes_molecule=False)
    reviewed = ["cells", "antibodies", "stain", "transposome", "tn5", "tagment",
                "bridge_oligos", "bridge_tags", "droplets", "partition", "split_branches"]
    atac_final = role_library("Phospho-seq ATAC-library roles",
                              (CELL, role("accessible-chromatin insert")), index=SAMPLE)
    tag_final = role_library("Phospho-seq tag-library roles",
                             (CELL, UMI, ANTIBODY, role("hash or phosphoprotein tag")), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final, "ATAC branch: " + d.endpoint,
                reviewed, label="Phospho-seq ATAC")
    _finish_pcr(graph, "tags", tags, tag_final, "Antibody/hash branch: " + d.endpoint,
                reviewed, label="Phospho-seq tag")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'phospho-seq'
    return _phospho_graph(d, final)

PHOSPHO_BOA = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
PHOSPHO_BOB = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC"
PHOSPHO_P5 = "AATGATACGGCGACCACCGA"
PHOSPHO_P7_TSA = "CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
PHOSPHO_P7_TSB = "CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGC"
PHOSPHO_HTO_ADDITIVE = "CCTTGGCACCCGAGAATTCC"
PHOSPHO_ADT_ADDITIVE = "GTGACTGGAGTTCAGACGTGTGCTC"
PHOSPHO_MULTIOME_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"

def _phospho_oligos():
    """Published Phospho-seq tag, bridge and amplification oligos.

    The bridge sequences are inherited explicitly from ASAP-seq, the method cited for
    this molecular conversion. Antibody-specific 15-base barcodes and sample-specific
    i7 choices remain variable windows rather than invented representative bases.
    """
    return (
        "<b>TotalSeq-A-style hashing tag</b> — "
        "5′-/5AmMC12/CCTTGGCACCCGAGAATTCCA[15-nt antibody barcode]B"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*A*A-3′",
        "<b>TotalSeq-B-style antibody tag</b> — "
        "5′-/5AmMC12/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[10-nt UMI]"
        "[15-nt antibody barcode][9-nt UMI]GCTTTAAGGCCGGTCCTAGC*A*A-3′",
        "<b>30-nt blocking oligo</b> — 5′-NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/3ddC/-3′",
        "<b>Bridge oligo A (BOA)</b> — 5′-" + PHOSPHO_BOA + "/3InvdT/-3′",
        "<b>Bridge oligo B (BOB)</b> — 5′-" + PHOSPHO_BOB + "/3InvdT/-3′",
        "<b>P5 tag-library primer</b> — 5′-" + PHOSPHO_P5 + "-3′",
        "<b>TotalSeq-A P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSA + "-3′",
        "<b>TotalSeq-B P7 index-primer family</b> — 5′-" + PHOSPHO_P7_TSB + "-3′",
        "<b>Multiome HTO additive</b> — 5′-" + PHOSPHO_HTO_ADDITIVE + "*C*C-3′",
        "<b>Multiome ADT additive</b> — 5′-GTGACTGGAGTTCAGACGTGTGC*T*C-3′",
        "<b>Multiome SI-PCR P5</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C-3′",
    )


def _phospho():
    return chain(molecule("fixed permeabilized hashed cell", (
        role("accessible chromatin"), role("intracellular phosphoproteins"))), (
        ("Stain intracellular targets with DNA-tagged antibodies",
         molecule("tagged cell", (ANTIBODY, role("phosphoprotein epitope"),
                                   ident("hash tag", "sample_hash", "sample_index"))),
         "Click-conjugated antibodies, including TotalSeq-A and TotalSeq-B formats, carry distinct feature barcodes."),
        ("Tagment chromatin and add bridge oligos before 10x ATAC GEM barcoding",
         molecule("capturable branches", (role("tagmented accessible DNA"),
                                           role("bridge-bound antibody tag"))),
         "Bridge oligos make the antibody tags compatible with the gel-bead capture sequence."),
        ("Separate ATAC DNA from the tag-containing supernatant after SPRI cleanup",
         molecule("matched Phospho-seq libraries", (CELL, role("ATAC insert"), UMI,
                                                     ANTIBODY, role("phosphoprotein tag"))),
         "ATAC fragments remain on beads; hashing and antibody tags are amplified from the saved supernatant."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (CELL, UMI, ANTIBODY, role("ATAC or protein insert")), index=SAMPLE)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT="The Phospho-seq paper prints the complete home-conjugated tag "
                       "architectures, tag-library PCR primers and run lengths. It names "
                       "BOA and BOB but inherits their exact sequences from the cited "
                       "ASAP-seq bridge method; those published sequences are listed "
                       "above. Antibody barcodes and i7 bases are intentionally variable.",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="The antibody-tag branch has published PCR "
                                       "oligos and a 34-base Read 1, 8-base i7, 16-base "
                                       "i5 and 34-base Read 2 run. Its cell-barcode end "
                                       "is supplied by the versioned commercial 10x "
                                       "scATAC gel-bead chemistry, so this page "
                                       "does not fabricate that proprietary oligo.",
    workflow=_phospho,
    oligos=_phospho_oligos,
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
