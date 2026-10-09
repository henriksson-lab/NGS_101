"""Visium Spatial Gene Expression (fresh-frozen workflow) architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, spatial_rt_scene, truseq_library
from chemdraw import Row, feature

TITLE = "10x Visium Spatial Gene Expression — array capture"
NOTES = "01_visium.html"
SOURCE = 'Commercial protocol: 10x Genomics <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000239_Visium_Spatial_Gene_Expression_User_Guide_Rev_F.pdf">CG000239 Rev F</a>.'
SUMMARY = "A tissue section is permeabilized over an array of spatially indexed oligo-dT spots. Released polyadenylated RNA is copied on the slide so each cDNA receives the spot barcode and a UMI before standard short-read library construction."
CAVEAT = "The spot barcode whitelist and complete surface oligo are commercial. Their published roles and lengths are drawn as placeholders; no undisclosed bases are invented."
SPATIAL=feature("spatial_barcode","spatial_barcode","whitelist",whitelist="10x-visium-spot-barcode-set")
UMI_FEATURE=feature("umi","umi","random")

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([
    seg("spatial barcode", "B" * 16, "cbc", placeholder=True, feature=SPATIAL),
    seg("UMI", "U" * 12, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT) junction", "T" * 12, placeholder=True),
    seg("captured cDNA", "X" * 38, placeholder=True)],
    "Visium gene-expression library", dual_index=False, inferred_adapters=True)
FINAL_CAPTION = "Read 1 contains the 16-nt spatial barcode followed by the 12-nt UMI; Read 2 enters the captured transcript. Adapter bases are canonical where the commercial guide specifies standard Illumina functions."
SEQUENCING_INTRO = "The 28-cycle Read 1 reports spot barcode plus UMI. Read 2 reports transcript sequence, and the i7 read identifies the sample library."

def sections():
    return [
        ("Place, image and permeabilize tissue over the capture array",
         [Row(chunks=[("tissue section:       [cell][cell][cell]", None, False)]),
          Row(chunks=[("                         ↓ released poly(A) RNA", None, False)]),
          Row(chunks=[("capture area:    ● ● ● ● ● spatially barcoded spots", "cbc", False)])],
         "Each printed spot carries many copies of one spatial-barcode species. Imaging preserves the mapping from spot coordinates to tissue morphology."),
        ("Capture RNA and reverse-transcribe on the slide",
         spatial_rt_scene(barcode_parts=(("16-nt spatial barcode",16),), umi=12,
                          surface="slide attachment").rows(),
         "The surface oligo contributes spatial barcode, UMI and oligo-dT. Reverse transcription copies the locally captured RNA while it remains associated with its spot."),
        ("Release cDNA and construct the sequencing library",
         [Row(chunks=[("barcoded cDNA → second strand → amplification → fragmentation", None, False)]),
          Row(chunks=[("→ end repair / dA-tail → adapter ligation → sample-index PCR", "me", False)])],
         "The spatial identity is already inside the cDNA before bulk recovery; downstream fragmentation must retain the barcode-bearing end used for Read 1."),
    ]
