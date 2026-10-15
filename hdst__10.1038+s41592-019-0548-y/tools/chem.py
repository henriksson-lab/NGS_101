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

TITLE = 'HDST — High-Definition Spatial Transcriptomics'
NOTES = '01_hdst.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0548-y">Vickovic et al., Nature Methods (2019)</a>.'
SUMMARY = "A decoded array of 2-micrometre wells holds randomly placed beads bearing three-part spatial barcodes, UMIs and poly(dT); captured transcripts inherit the bead's imaged coordinate."
CAVEAT = '<inf>Individual barcode pools and decoder sequences were proprietary. Their published lengths and Read 1 coordinates are preserved, while unavailable bridge bases and the compatible completed Illumina shell are marked inferred.</inf>'
FINAL_CAPTION = 'Read 1 encodes the three barcode blocks and six-base UMI; Read 2 is transcript sequence. The representative construct marks unavailable proprietary bridge and barcode-panel bases.'
SEQUENCING_INTRO = 'Read 1 encodes the three barcode blocks and six-base UMI; Read 2 is transcript sequence. The representative construct marks unavailable proprietary bridge and barcode-panel bases.'
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: ('Bead primer precursor, as written — /AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1',)
FINAL_LIBRARY, SEQ_PRIMERS = hdst_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1038/s41592-019-0548-y', citation='Vickovic et al., Nature Methods (2019)', endpoint='Read 1 encodes the three barcode blocks and six-base UMI; Read 2 is transcript sequence. The representative construct marks unavailable proprietary bridge and barcode-panel bases.')


def _hdst():
    return chain(molecule("2-micrometre amine-functionalized silica bead", (
        role("surface amines"),)), (
        ("Couple the uracil-bearing T7 / Illumina precursor with barcode pool 1",
         molecule("first-round bead oligo", (role("bead linkage and U5 release site"),
             role("T7 promoter and Illumina handle"),
             ident("15-nt barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("bridge 1"))),
         "The published precursor skeleton is exact; the individual Illumina-supplied barcode and decoder sequences are proprietary."),
        ("Split, ligate barcode pools 2 and 3, then append UMI and poly(dT) capture",
         molecule("completed bead capture oligo", (
             ident("barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("bridge 1"), ident("barcode pool 2", "hdst.bc2", "spatial_barcode", length=15),
             role("bridge 2"), ident("barcode pool 3", "hdst.bc3", "spatial_barcode", length=15),
             ident("six-base UMI", "hdst.umi", "umi", length=6, tag="umi"),
             role("poly(dT) capture"))),
         "The three split-and-pool rounds provide 65 × 211 × 211 combinations. Every oligo on one bead shares its three-part spatial identity but has an independent UMI."),
        ("Load beads randomly into 2-micrometre wells and decode each position",
         molecule("decoded HDST array", (role("known x-y well"),
             role("one three-part spatial barcode"), role("many capture oligos"))),
         "Fourteen cycles of FAM, Cy3 or dark decoder hybridization map each randomly placed bead to an image coordinate before tissue is applied."),
        ("Place, image and permeabilize tissue; reverse-transcribe captured poly(A) RNA",
         molecule("bead-bound spatial cDNA", (
             role("bead"), role("three-part spatial barcode"),
             role("six-base UMI"), role("poly(dT)"), role("antisense transcript cDNA"))),
         "Superscript III extends the bead oligo overnight. Barcode, UMI and transcript information are now covalently linked in one first strand."),
        ("Digest tissue, release at the five-uracil site with USER and prepare RNA-seq library",
         molecule("paired-end HDST library", (
             role("Read 1 arm"), role("three barcode blocks and bridges"),
             role("six-base UMI"), role("transcript insert"), role("Read 2 and indexed arm"))),
         "The released cDNA follows the published Spatial Transcriptomics CEL-seq-like amplification and library workflow. Read 1 reports position and UMI; Read 2 reports transcript."),
    ))

workflow = _hdst

def _hdst_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    bead = graph.input("bead", "amine-functionalized silica bead bearing the T7/Illumina precursor",
                       MaterialKind.BIOLOGICAL, state=_hdst().initial, primary=True)
    pools = graph.input("barcode_pools", "three split-and-pool barcode/bridge oligo pools",
                        MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "barcode-assembly DNA ligase", MaterialKind.CATALYST)
    capture_bead = graph.product("capture_bead", "bead bearing three-part spatial barcode, UMI and poly(dT)",
                                 MaterialKind.MOLECULAR,
                                 state=molecule("completed bead capture oligo", (
                                     role("three barcode blocks and bridges"), UMI,
                                     role("poly(dT) capture"))))
    graph.react("assemble_bead", "Ligate the split-and-pool barcode rounds and capture tail",
                StepKind.LIGATION,
                ((bead, InputRole.TEMPLATE), (pools, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (capture_bead,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The public source gives the precursor skeleton and component lengths but withholds the proprietary barcode and decoder sequences."))
    decoder = graph.input("decoder_oligos", "fluorescent decoder oligo series",
                          MaterialKind.OLIGO)
    decoded = graph.product("decoded_array", "position-decoded HDST bead array",
                            MaterialKind.MOLECULAR,
                            state=molecule("decoded HDST array", (
                                role("known x-y well"), role("three-part spatial barcode"))))
    graph.react("decode", "Decode bead positions by iterative hybridization",
                StepKind.HYBRIDIZATION,
                ((capture_bead, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (decoded,))
    tissue_rna = graph.input("tissue_rna", "poly(A) RNA released from the overlaid tissue section",
                             MaterialKind.MOLECULAR)
    captured = graph.product("captured_rna", "spatially barcoded bead-bound RNA",
                             MaterialKind.MOLECULAR,
                             state=molecule("RNA captured on bead", (SPATIAL_BEAD, UMI,
                                 role("poly(dT)"), role("poly(A) RNA"))))
    graph.react("capture_rna", "Hybridize released tissue RNA to bead poly(dT)",
                StepKind.HYBRIDIZATION,
                ((decoded, InputRole.TEMPLATE), (tissue_rna, InputRole.STRUCTURAL)), (captured,))
    rt_primer = graph.input("rt_primer", "bead-bound poly(dT) RT primer", MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "Superscript III reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product("cdna", "bead-bound spatial cDNA", MaterialKind.MOLECULAR,
                         state=molecule("bead-bound spatial cDNA", (SPATIAL_BEAD, UMI,
                             role("antisense transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe captured poly(A) RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    user = graph.input("user_enzyme", "USER release enzyme", MaterialKind.CATALYST)
    released = graph.product("released_cdna", "cDNA released at the five-uracil site",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("release", "Release cDNA at the five-uracil site", StepKind.DIGESTION,
                ((cdna, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (released,))
    reviewed = ["bead", "barcode_pools", "ligase", "assemble_bead", "decoder_oligos",
                "decode", "tissue_rna", "capture_rna", "rt_primer", "reverse_transcriptase",
                "reverse_transcribe", "user_enzyme", "release"]
    _finish_pcr(graph, "library", released, final, primers, d.endpoint, reviewed,
                label="HDST")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _hdst_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
