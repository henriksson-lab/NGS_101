"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (ident, SPATIAL_BEAD, UMI, _finish_pcr,
    hdst_library, HDST_BRIDGE12, HDST_BRIDGE34)
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import MolecularState, Scene, Workflow, duplex_rows
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'HDST — High-Definition Spatial Transcriptomics'
NOTES = '01_hdst.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0548-y">Vickovic et al., Nature Methods (2019)</a>.'
SUMMARY = "A decoded array of 2-micrometre wells holds randomly placed beads bearing three-part spatial barcodes, UMIs and poly(dT); captured transcripts inherit the bead's imaged coordinate."
CAVEAT = ('<inf>Individual barcode pools and decoder sequences were proprietary; the completed Illumina shell is inferred.</inf> '
          'Published Methods disagree: Bead production (JATS line 121) specifies 15/14/15-base barcode pools and a 5-base UMI, '
          'while Sequence preprocessing (JATS line 145) extracts Read 1 bases 1–15, 31–45, 61–76 and 77–82 as 15/15/16-base barcodes and a 6-base UMI. '
          'The authors’ bioRxiv v2 Bead production (PDF text lines 220–235) instead says a 7-base UMI; its preprocessing (lines 302–319) also selects bases 77–82. '
          'The molecular drawing follows the final published synthesis recipe; both read masks are retained as written. '
          'The offset between them remains unresolved.')
FINAL_CAPTION = 'The physical assembly follows the Bead production Methods (15/14/15-base barcode pools, published bridges and 5-base UMI); the reported Read 1 mask uses different positions and a six-base UMI. Read 2 is transcript sequence.'
SEQUENCING_INTRO = 'The paper reports Read 1 barcode and UMI positions 1–15, 31–45, 61–76 and 77–82, respectively. These coordinates conflict with the published synthesis lengths; the physical construct below follows the synthesis recipe.'
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: (
    'Bead primer precursor, as written — /AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1',
    'Published ligated bridge 1+2 — ' + HDST_BRIDGE12,
    'Published ligated bridge 3+4 — ' + HDST_BRIDGE34,
)
FINAL_LIBRARY, SEQ_PRIMERS = hdst_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1038/s41592-019-0548-y', citation='Vickovic et al., Nature Methods (2019)', endpoint=FINAL_CAPTION)


def _hdst():
    return chain(molecule("2-micrometre amine-functionalized silica bead", (
        role("surface amines"),)), (
        ("Couple the uracil-bearing T7 / Illumina precursor with barcode pool 1",
         molecule("first-round bead oligo", (role("bead linkage and U5 release site"),
             role("T7 promoter and Illumina handle"),
             ident("15-nt barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("bridge 1"))),
         "The published precursor skeleton is exact; the individual Illumina-supplied barcode and decoder sequences are proprietary."),
        ("Split and pool; splint-ligate 5'-phosphorylated barcode pools 2 and 3 (pool 3 is followed by the UMI and dT20VN)",
         molecule("completed bead capture oligo", (
             ident("barcode pool 1", "hdst.bc1", "spatial_barcode", length=15),
             role("ligated bridge 1+2"), ident("barcode pool 2", "hdst.bc2", "spatial_barcode", length=14),
             role("ligated bridge 3+4"), ident("barcode pool 3", "hdst.bc3", "spatial_barcode", length=15),
             ident("five-base UMI", "hdst.umi", "umi", length=5, tag="umi"),
             role("poly(dT) capture"))),
         "T4 DNA ligase joins /Pho/Bridge2-Pool2-Bridge3 (2:1 over precursor) through a 14-nt Bridge1'Bridge2' helper, then /Pho/Bridge4-Pool3 through a Bridge3'Bridge4' helper. The three split-and-pool rounds provide 65 × 211 × 211 combinations. Every oligo on one bead shares its three-part spatial identity but has an independent UMI."),
        ("Load beads randomly into 2-micrometre wells and decode each position",
         molecule("decoded HDST array", (role("known x-y well"),
             role("one three-part spatial barcode"), role("many capture oligos"))),
         "Fourteen cycles of FAM, Cy3 or dark decoder hybridization map each randomly placed bead to an image coordinate before tissue is applied."),
        ("Place, image and permeabilize tissue; reverse-transcribe captured poly(A) RNA",
         molecule("bead-bound spatial cDNA", (
             role("bead"), role("three-part spatial barcode"),
             role("five-base UMI"), role("poly(dT)"), role("antisense transcript cDNA"))),
         "Superscript III extends the bead oligo overnight. Barcode, UMI and transcript information are now covalently linked in one first strand."),
        ("Proteinase K-digest tissue, release at the five-uracil site with USER and prepare RNA-seq library",
         molecule("paired-end HDST library", (
             role("Read 1 arm"), role("three barcode blocks and bridges"),
             role("five-base UMI"), role("transcript insert"), role("Read 2 and indexed arm"))),
         "The cited Jemt et al. workflow performs second-strand synthesis, end repair, T7 in vitro transcription, RNA-adapter ligation, a second cDNA synthesis and indexing PCR. Its exact oligos are delegated to Ståhl et al. (2016); the HDST-specific final junction remains unresolved. Read 1 reports position and UMI; Read 2 reports transcript."),
    ))

workflow = _hdst


def _round2_geometry():
    """Round-2 splint ligation end as described in Online Methods "Bead production".

    The precursor's Bridge1 is stated as 7 bp and the ligated Bridge1Bridge2 is printed as
    14 bases, so the 5'-phosphorylated Bridge2 donor end is the last seven printed bases
    (a derived split, not a printed oligo).  The Bridge3/Bridge4 split is not stated.
    """
    return reported_adapter_end(
        "pool-2 oligo /Pho/Bridge2-Spatial_barcode_Pool2-Bridge3", kind="splinted",
        polarity="5-prime", terminal_sequence=HDST_BRIDGE12[7:],
        statement=("Online Methods 'Bead production': /Pho/Bridge2-Pool2-Bridge3 is ligated to the bead's "
                   "Bridge1 through a complementary 14-bp Bridge1'Bridge2' helper with T4 DNA ligase, and "
                   "/Pho/Bridge4-Pool3 through Bridge3'Bridge4'. Bridge2 = last 7 of the printed "
                   "GACTTGTCTAGAGC (derived from the stated 7-bp Bridge1); helper bases are not printed and "
                   "the Bridge3/Bridge4 split and barcode sequences are not disclosed."))

def _hdst_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    bead = graph.input("bead", "amine-functionalized silica bead bearing the T7/Illumina precursor",
                       MaterialKind.BIOLOGICAL, state=_hdst().initial, primary=True)
    pools = graph.input("barcode_pools", "three split-and-pool barcode/bridge oligo pools",
                        MaterialKind.ADAPTER)
    helpers = graph.input("bridge_helpers", "Bridge1'Bridge2' (14 nt) and Bridge3'Bridge4' helper oligos",
                          MaterialKind.OLIGO)
    ligase = graph.input("ligase", "T4 DNA ligase (NEB)", MaterialKind.CATALYST)
    capture_bead = graph.product("capture_bead", "bead bearing three-part spatial barcode, UMI and poly(dT)",
                                 MaterialKind.MOLECULAR,
                                 state=molecule("completed bead capture oligo", (
                                     role("three barcode blocks and bridges"), UMI,
                                     role("poly(dT) capture"))))
    graph.react("assemble_bead", "Ligate the split-and-pool barcode rounds and capture tail",
                StepKind.LIGATION,
                ((bead, InputRole.TEMPLATE), (pools, InputRole.STRUCTURAL),
                 (helpers, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (capture_bead,),
                geometry=_round2_geometry())
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
    proteinase = graph.input("proteinase_k", "proteinase K (Qiagen), 1 h at 56 °C", MaterialKind.CATALYST)
    user = graph.input("user_enzyme", "USER enzyme (NEB), 2 h at 37 °C", MaterialKind.CATALYST)
    released = graph.product("released_cdna", "cDNA released at the five-uracil site",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("release", "Digest tissue and release cDNA at the five-uracil site", StepKind.DIGESTION,
                ((cdna, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST),
                 (user, InputRole.CATALYST)), (released,))
    processed = graph.product("adapter_cdna", "second-round cDNA from adapter-ligated aRNA (Jemt et al. workflow)",
                              MaterialKind.MOLECULAR,
                              state=molecule("second-round HDST cDNA", (
                                  role("barcode-side handle"), SPATIAL_BEAD, UMI,
                                  role("transcript sequence"), role("ligated RNA-adapter sequence"))))
    graph.react("jemt_processing",
                "Second-strand synthesis, end repair, T7 IVT, RNA-adapter ligation and second cDNA synthesis",
                StepKind.MAGIC, ((released, InputRole.TEMPLATE),), (processed,),
                reason=("HDST Online Methods delegate library processing to Jemt et al. (2016), whose Methods "
                        "(line 89) in turn delegate step details and oligo sequences to the Ståhl et al. (2016) "
                        "supplement. Those adapter, RT and indexing oligos were not inspected, so the molecular "
                        "transitions and end geometry are not reconstructed."))
    reviewed = ["bead", "barcode_pools", "bridge_helpers", "ligase", "assemble_bead", "decoder_oligos",
                "decode", "tissue_rna", "capture_rna", "rt_primer", "reverse_transcriptase",
                "reverse_transcribe", "proteinase_k", "user_enzyme", "release",
                "adapter_cdna", "jemt_processing"]
    _finish_pcr(graph, "library", processed, final, primers, d.endpoint, reviewed,
                label="HDST")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _hdst_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
