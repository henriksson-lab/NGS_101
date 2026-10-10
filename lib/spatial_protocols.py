"""Shared molecular models for spatial capture and sample hashing protocols."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import feature
from protocol_model import chain, molecule, role, role_library


def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))


SPATIAL_A = ident("spatial barcode A", "spatial_a", "spatial_barcode")
SPATIAL_B = ident("spatial barcode B", "spatial_b", "spatial_barcode")
SPATIAL_BEAD = ident("decoded bead barcode", "spatial_bead", "spatial_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
HASH = ident("nuclear hash", "sample_hash", "sample_index")
CELL_A = ident("first-round cell barcode", "cell_round1", "cell_barcode")
CELL_B = ident("second-round cell barcode", "cell_round2", "cell_barcode")


def _dbit(*, cuttag: bool):
    target = (role("antibody-bound histone mark"), role("poly(A) RNA")) if cuttag else (
        role("accessible chromatin"), role("poly(A) RNA"))
    enzyme = "antibody-tethered pA–Tn5" if cuttag else "Tn5"
    dna_kind = "histone-mark DNA" if cuttag else "accessible DNA"
    return chain(molecule("fixed permeabilized tissue section", target), (
        (f"Tag {dna_kind} in situ with {enzyme}",
         molecule("tagged DNA branch", (role("universal ligation linker"), role(dna_kind),
                                         role("transposase end"))),
         "The DNA linker is the landing site for the first spatial barcode."),
        ("Reverse-transcribe poly(A) RNA with a biotinylated poly(dT) primer",
         molecule("DNA and cDNA branches in tissue", (role("tagged DNA branch"),
             role("biotinylated RT handle"), UMI, role("transcript cDNA"))),
         "DNA and RNA are still physically registered in the same tissue section."),
        ("Flow the first barcode set through parallel microfluidic channels",
         molecule("A-barcoded branches", (SPATIAL_A, role("DNA or cDNA insert"))),
         "Barcode A records one tissue coordinate on both molecular branches."),
        ("Rotate the device and ligate the orthogonal barcode set",
         molecule("two-dimensional pixel address", (SPATIAL_A, SPATIAL_B,
                                                       role("DNA or cDNA insert")),
                  junctions=(("spatial barcode A", "spatial barcode B", "ligation"),)),
         "The ordered A–B pair identifies the channel intersection, not an inferred cell."),
        ("Reverse crosslinks, recover molecules and build separate DNA and RNA libraries",
         molecule("matched spatial libraries", (SPATIAL_A, SPATIAL_B, UMI,
                                                   role(f"{dna_kind} or transcript insert"))),
         "The physical libraries are separated while retaining the same two-dimensional address."),
    ))


def _slide_seq():
    return chain(molecule("10-µm bead puck", (SPATIAL_BEAD, UMI, role("poly(dT) capture oligo"))), (
        ("Decode each bead position by SOLiD sequencing-by-ligation",
         molecule("mapped bead puck", (SPATIAL_BEAD, role("known x–y coordinate"))),
         "The decoded bead sequence is a spatial address before tissue is applied."),
        ("Transfer a fresh-frozen tissue section and release RNA",
         molecule("RNA captured on bead", (SPATIAL_BEAD, UMI, role("poly(dT)"),
                                             role("poly(A) RNA"))),
         "RNA inherits position by hybridizing to an oligo on a mapped bead."),
        ("Reverse-transcribe and prepare a 3′ barcoded RNA-seq library",
         molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI, role("3-prime transcript cDNA"))),
         "The original method follows the Drop-seq-style 3′ capture logic."),
        ("Amplify, tagment and add the sample index",
         molecule("Slide-seq sequencing library", (SPATIAL_BEAD, UMI,
             role("3-prime transcript insert"), SAMPLE)),
         "Sequencing joins each transcript molecule to the bead position decoded on the puck."),
    ))


def _sciplex_atac():
    return chain(molecule("permeabilized nuclei separated by sample", (role("accessible chromatin"),)), (
        ("Let nuclei absorb a sample-specific unmodified ssDNA hash oligo",
         molecule("hashed nucleus", (HASH, role("accessible chromatin"))),
         "The hash records specimen or treatment; it is not a cell barcode."),
        ("Extend the polyadenylated hash and tagment accessible chromatin in indexed wells",
         molecule("first-round indexed nucleus", (HASH, CELL_A, role("tagmented accessible DNA"))),
         "Hash extension makes the label recoverable alongside the ATAC assay."),
        ("Pool, redistribute and add the second combinatorial index",
         molecule("two-level cell identity", (HASH, CELL_A, CELL_B,
                                                role("accessible-DNA fragment"))),
         "The pair of well indices identifies a nucleus; the independent hash identifies its sample."),
        ("Amplify ATAC fragments and hash reads",
         molecule("sciPlex-ATAC readouts", (HASH, CELL_A, CELL_B,
                                             role("ATAC insert"), SAMPLE)),
         "Two-level experiments sequence hash and ATAC products separately; later three-level variants can intermix them."),
    ))


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str


DEFS = {
    "spatial-atac-rna": Definition("spatial-ATAC–RNA-seq", "10.1038/s41586-023-05795-1",
        "Zhang et al., Nature (2023)", "Tn5 tags accessible DNA and in-situ poly(dT) RT captures RNA before orthogonal microfluidic barcode ligations assign both branches the same tissue pixel.",
        "Separate ATAC and RNA libraries retain spatial barcode A plus B; the RNA branch additionally carries a UMI.", "01_spatial_atac_rna_seq.html"),
    "spatial-cuttag-rna": Definition("spatial-CUT&Tag–RNA-seq", "10.1038/s41586-023-05795-1",
        "Zhang et al., Nature (2023)", "Antibody-tethered pA–Tn5 tags histone-marked DNA and in-situ poly(dT) RT captures RNA before both branches receive orthogonal spatial barcodes.",
        "Separate CUT&Tag and RNA libraries retain spatial barcode A plus B; the RNA branch additionally carries a UMI.", "01_spatial_cuttag_rna_seq.html"),
    "slide-seq": Definition("Slide-seq (v1)", "10.1126/science.aaw1219",
        "Rodriques et al., Science (2019)", "A pre-decoded monolayer of 10-µm DNA-barcoded beads captures RNA released from a tissue section for spatially registered 3′ RNA sequencing.",
        "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript.", "01_slide_seq.html"),
    "sciplex-atac": Definition("sciPlex-ATAC-seq", "10.1101/2023.03.05.531201",
        "Booth et al., bioRxiv (2023; Genome Research 2023)", "Permeabilized nuclei absorb sample-specific ssDNA hashes, then two-level combinatorial ATAC indexing assigns a separate cell identity.",
        "Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.", "01_sciplex_atac_seq.html"),
}


def get(key: str):
    d = DEFS[key]
    workflows = {
        "spatial-atac-rna": lambda: _dbit(cuttag=False),
        "spatial-cuttag-rna": lambda: _dbit(cuttag=True),
        "slide-seq": _slide_seq,
        "sciplex-atac": _sciplex_atac,
    }
    payloads = {
        "spatial-atac-rna": (SPATIAL_A, SPATIAL_B, UMI, role("ATAC or RNA insert")),
        "spatial-cuttag-rna": (SPATIAL_A, SPATIAL_B, UMI, role("CUT&Tag or RNA insert")),
        "slide-seq": (SPATIAL_BEAD, UMI, role("3-prime transcript insert")),
        "sciplex-atac": (HASH, CELL_A, CELL_B, role("ATAC insert")),
    }
    final = role_library(d.title + " final library roles", payloads[key], index=SAMPLE)
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT="", FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE="The defining source establishes the read roles; kit revisions and barcode panels prevent one universal sequencing-primer sequence from being asserted here.",
        workflow=workflows[key], oligos=lambda: (),
    )
