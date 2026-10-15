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

TITLE = 'Slide-seq (v1)'
NOTES = '01_slide_seq.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.aaw1219">Rodriques et al., Science (2019)</a>.'
SUMMARY = 'A pre-decoded monolayer of 10-µm DNA-barcoded beads captures RNA released from a tissue section for spatially registered 3′ RNA sequencing.'
CAVEAT = ''
FINAL_CAPTION = "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript."
SEQUENCING_INTRO = "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript."
READ_LENGTHS = {'Read 1': 42, 'Index 1 (i7)': 8, 'Read 2': 50}
READ_USES = {'Read 1': 'extract', 'Index 1 (i7)': 'demultiplex', 'Read 2': 'include'}
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: ("<p>Bead oligo sequence 2 (representative): 5'- TTTTTTCTACACGACGCTCTTCCGATCT<cbc>JJJJJJJJ</cbc>TCTTCAGCGTTCCCGAGA<cbc>JJJJJJJ</cbc><umi>NNNNNNNN</umi>TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -3'</p>", "<p>SMART PCR primer: 5'- AAGCAGTGGTATCAACGCAGAGT -3'</p>", "<p>TruSeq PCR handle primer: 5'- CTACACGACGCTCTTCCGATCT -3'</p>", "<p>TruSeq5 tagmentation PCR primer: 5'- AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'</p>")
FINAL_LIBRARY, SEQ_PRIMERS = slide_seq_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1126/science.aaw1219', citation='Rodriques et al., Science (2019)', endpoint="The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript.")


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

workflow = _slide_seq

def _slide_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    puck = graph.input("puck", "DNA-barcoded 10-micrometre bead puck",
                       MaterialKind.BIOLOGICAL, state=_slide_seq().initial, primary=True)
    decoder = graph.input("decoder", "SOLiD sequencing-by-ligation decoder probes",
                          MaterialKind.OLIGO)
    mapped = graph.product("mapped_puck", "position-decoded bead puck",
                           MaterialKind.MOLECULAR,
                           state=molecule("mapped bead puck", (SPATIAL_BEAD,
                               role("known x-y coordinate"))))
    graph.react("decode", "Decode each bead position by sequencing-by-ligation",
                StepKind.HYBRIDIZATION,
                ((puck, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (mapped,))
    rna = graph.input("tissue_rna", "poly(A) RNA released from the tissue section",
                      MaterialKind.MOLECULAR)
    captured = graph.product("captured_rna", "RNA captured on mapped beads",
                             MaterialKind.MOLECULAR,
                             state=molecule("RNA captured on bead", (SPATIAL_BEAD, UMI,
                                 role("poly(A) RNA"))))
    graph.react("capture", "Capture tissue RNA on mapped bead oligo(dT)",
                StepKind.HYBRIDIZATION,
                ((mapped, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL)), (captured,))
    rt_primer = graph.input("rt_primer", "bead oligo(dT) RT primer", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("cdna", "bead-barcoded 3-prime cDNA", MaterialKind.MOLECULAR,
                         state=molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI,
                             role("3-prime transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    transposome = graph.input("transposome", "library-adapter-loaded transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented Slide-seq cDNA", MaterialKind.MOLECULAR,
                           state=cdna.state)
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified bead-barcoded cDNA", StepKind.TAGMENTATION,
                ((cdna, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain an amplifiable spatial transcript insert")
    reviewed = ["puck", "decoder", "decode", "tissue_rna", "capture", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, final, primers, d.endpoint, reviewed,
                label="Slide-seq")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _slide_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
