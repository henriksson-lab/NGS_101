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

TITLE = 'sci-Space'
NOTES = '01_sci_space.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.abb9536">Srivatsan et al., Science (2021)</a>.'
SUMMARY = 'Polyadenylated DNA hashes diffuse from a registered printed grid into permeabilized tissue nuclei; sci-RNA-seq then associates each grid identity with a single-nucleus transcriptome.'
CAVEAT = 'The printed hash, section hash and transcript are separate molecules that are associated through shared nucleus indices; the final construct shows the spatial-hash branch, not a fictitious covalent fusion to a transcript.'
FINAL_CAPTION = 'Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.'
SEQUENCING_INTRO = 'Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.'
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: ()
FINAL_LIBRARY, SEQ_PRIMERS = sci_space_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1126/science.abb9536', citation='Srivatsan et al., Science (2021)', endpoint='Separate transcript, printed-grid-hash and section-hash molecules share the same combinatorial nucleus indices; the representative final construct below is a spatial-hash molecule.')


def _sci_space():
    return chain(molecule("fresh-frozen tissue section", (
        role("poly(A) RNA in intact nuclei"), role("anatomical position"))), (
        ("Permeabilize the section and add one polyadenylated section hash",
         molecule("section-labelled nuclei in tissue", (SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "The section hash identifies the slide after sections are pooled; it is distinct from the printed coordinate hashes."),
        ("Juxtapose tissue with the agarose slide carrying a gridded hash pattern",
         molecule("spatially hashed nuclei", (SPACE_HASH, SECTION_HASH,
             role("poly(A) RNA in intact nuclei"))),
         "Single-stranded poly(A)-tailed hash oligos diffuse from each printed spot into nearby permeabilized nuclei. Fluorescent fiducials register the grid to the tissue image."),
        ("Recover nuclei, fix them and perform barcoded sci-RNA-seq reverse transcription",
         molecule("first-round indexed nuclear molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, role("transcript cDNA or copied hash"))),
         "The same oligo(dT) reaction copies endogenous mRNA, the grid hash and the section hash, associating all three with a nucleus-specific RT-well barcode."),
        ("Pool, redistribute, second-strand synthesize, tagment and PCR-index",
         molecule("sci-Space sequencing molecules", (
             ident("RT-well cell barcode", "sci_space.rt", "cell_barcode"),
             UMI, SPACE_HASH, role("transcript or hash payload"), SAMPLE)),
         "The downstream combinatorial sci-RNA-seq indices identify one nucleus. Hash molecules and transcript molecules are separate physical reads joined computationally by that shared nucleus identity."),
    ))

workflow = _sci_space

def _sci_space_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    tissue = graph.input("tissue", "fresh-frozen tissue section containing intact nuclei",
                         MaterialKind.BIOLOGICAL, state=_sci_space().initial, primary=True)
    section_hash = graph.input("section_hash", "polyadenylated section-hash oligo",
                               MaterialKind.OLIGO)
    sectioned = graph.product("sectioned", "section-labelled tissue nuclei",
                              MaterialKind.MOLECULAR,
                              state=molecule("section-labelled nuclei", (SECTION_HASH,
                                  role("poly(A) RNA"))))
    graph.react("label_section", "Hybridize the section hash in permeabilized tissue",
                StepKind.HYBRIDIZATION,
                ((tissue, InputRole.TEMPLATE), (section_hash, InputRole.STRUCTURAL)),
                (sectioned,))
    grid_hash = graph.input("grid_hash", "registered printed spatial-hash grid",
                            MaterialKind.OLIGO)
    spatial = graph.product("spatial_nuclei", "nuclei carrying section and grid hashes",
                            MaterialKind.MOLECULAR,
                            state=molecule("spatially hashed nuclei", (SPACE_HASH,
                                SECTION_HASH, role("poly(A) RNA"))))
    graph.react("transfer_grid", "Transfer printed spatial hashes into nearby nuclei",
                StepKind.HYBRIDIZATION,
                ((sectioned, InputRole.TEMPLATE), (grid_hash, InputRole.STRUCTURAL)),
                (spatial,))
    rt_primer = graph.input("rt_primer", "barcoded sci-RNA-seq oligo(dT) RT primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("first_strand", "RT-well-indexed transcript and hash cDNA",
                          MaterialKind.MOLECULAR,
                          state=molecule("first-round indexed nuclear molecules", (
                              CELL_A, UMI, role("transcript cDNA or copied hash"))))
    graph.react("reverse_transcribe", "Copy endogenous RNA and both hash classes",
                StepKind.REVERSE_TRANSCRIPTION,
                ((spatial, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_primer = graph.input("second_strand_primer", "second-strand synthesis primer",
                                MaterialKind.OLIGO)
    polymerase = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    duplex = graph.product("duplex", "double-stranded indexed transcript/hash products",
                           MaterialKind.MOLECULAR, state=first.state)
    graph.react("second_strand", "Synthesize the second strand after pool redistribution",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (duplex,))
    transposome = graph.input("transposome", "indexed library transposome", MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "tagmented sci-Space products", MaterialKind.MOLECULAR,
                           state=duplex.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment transcript and hash products", StepKind.TAGMENTATION,
                ((duplex, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable transcript or hash molecule")
    reviewed = ["tissue", "section_hash", "label_section", "grid_hash", "transfer_grid",
                "rt_primer", "reverse_transcriptase", "reverse_transcribe",
                "second_strand_primer", "second_strand_polymerase", "second_strand",
                "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, final, primers, d.endpoint, reviewed,
                label="sci-Space")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _sci_space_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
