"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (ident, SPATIAL_BEAD, UMI, CELL_A, CELL_B, HASH, SAMPLE,
    SECTION_HASH, SPACE_HASH, SPATIAL_A, SPATIAL_B, SPATIAL_A_ROUND, SPATIAL_B_ROUND,
    _finish_pcr, _spatial_dna_libraries, _spatial_rows, hdst_library, slide_seq_library,
    sci_space_library)
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState, Scene, Workflow, duplex_rows
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'sciPlex-ATAC-seq'
NOTES = '01_sciplex_atac_seq.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1101/2023.03.05.531201">Booth et al., bioRxiv (2023; Genome Research 2023)</a>.'
SUMMARY = 'Permeabilized nuclei absorb sample-specific ssDNA hashes, then two-level combinatorial ATAC indexing assigns a separate cell identity.'
CAVEAT = ''
FINAL_CAPTION = 'Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.'
SEQUENCING_INTRO = 'Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.'
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = 'Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.'
SEQUENCING_UNAVAILABLE = 'The defining source establishes the read roles; kit revisions and barcode panels prevent one universal sequencing-primer sequence from being asserted here.'
oligos = lambda: ()
FINAL_LIBRARY = role_library(TITLE + " final library roles", (HASH, CELL_A, CELL_B, role("ATAC insert")), index=SAMPLE)
SEQ_PRIMERS = ()
_definition = SimpleNamespace(title=TITLE, doi='10.1101/2023.03.05.531201', citation='Booth et al., bioRxiv (2023; Genome Research 2023)', endpoint='Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.')


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

workflow = _sciplex_atac

def _sciplex_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "permeabilized nuclei separated by sample",
                         MaterialKind.BIOLOGICAL, state=_sciplex_atac().initial, primary=True)
    hash_oligo = graph.input("hash_oligo", "sample-specific polyadenylated ssDNA hash",
                             MaterialKind.OLIGO)
    hashed = graph.product("hashed", "nuclei carrying absorbed sample hashes",
                           MaterialKind.MOLECULAR,
                           state=molecule("hashed nucleus", (HASH, role("accessible chromatin"))))
    graph.react("hash", "Let nuclei absorb the sample-specific hash oligo",
                StepKind.HYBRIDIZATION,
                ((nuclei, InputRole.TEMPLATE), (hash_oligo, InputRole.STRUCTURAL)), (hashed,))
    hash_primer = graph.input("hash_primer", "hash extension primer", MaterialKind.OLIGO)
    hash_polymerase = graph.input("hash_polymerase", "hash extension polymerase",
                                  MaterialKind.CATALYST)
    extended = graph.product("extended", "extended recoverable hash in each nucleus",
                             MaterialKind.MOLECULAR, state=hashed.state)
    graph.react("extend_hash", "Extend the absorbed sample hash", StepKind.EXTENSION,
                ((hashed, InputRole.TEMPLATE), (hash_primer, InputRole.STRUCTURAL),
                 (hash_polymerase, InputRole.CATALYST)), (extended,))
    transposome = graph.input("transposome", "first-well-indexed ATAC transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "sample-hashed first-round-indexed nuclei",
                           MaterialKind.MOLECULAR,
                           state=molecule("first-round indexed nucleus", (HASH, CELL_A,
                               role("tagmented accessible DNA"))))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment chromatin with a first-round well index",
                StepKind.TAGMENTATION,
                ((extended, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain a recoverable hash or ATAC insert")
    second = graph.input("second_round", "second-round combinatorial-index wells",
                         MaterialKind.CONDITION)
    indexed = graph.product("indexed", "two-level-indexed hashed nuclei",
                            MaterialKind.MOLECULAR,
                            state=molecule("two-level cell identity", (HASH, CELL_A, CELL_B,
                                role("accessible-DNA fragment"))))
    collisions = graph.product("collisions", "ambiguous combinatorial-index collisions",
                               MaterialKind.WASTE)
    graph.react("partition", "Pool and redistribute into second-index wells",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (second, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not have an unambiguous two-level nucleus identity")
    atac = graph.product("atac_branch", "ATAC fragment branch", MaterialKind.MOLECULAR,
                         state=indexed.state)
    hashes = graph.product("hash_branch", "sample-hash readout branch", MaterialKind.MOLECULAR,
                           state=indexed.state)
    graph.react("split_readouts", "Separate ATAC-fragment and hash readouts",
                StepKind.SPLIT, ((indexed, InputRole.TEMPLATE),), (atac, hashes),
                changes_molecule=False)
    reviewed = ["nuclei", "hash_oligo", "hash", "hash_primer", "hash_polymerase",
                "extend_hash", "transposome", "tn5", "tagment", "second_round",
                "partition", "split_readouts"]
    atac_final = role_library("sciPlex-ATAC fragment roles",
                              (HASH, CELL_A, CELL_B, role("ATAC insert")), index=SAMPLE)
    hash_final = role_library("sciPlex-ATAC hash-read roles",
                              (HASH, CELL_A, CELL_B), index=SAMPLE)
    _finish_pcr(graph, "atac", atac, atac_final, (), "ATAC branch: " + d.endpoint,
                reviewed, label="sciPlex-ATAC")
    _finish_pcr(graph, "hashes", hashes, hash_final, (), "hash branch: " + d.endpoint,
                reviewed, label="sciPlex hash")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _sciplex_graph(_definition, FINAL_LIBRARY)
