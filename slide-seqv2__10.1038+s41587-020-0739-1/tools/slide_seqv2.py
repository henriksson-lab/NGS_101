"""Slide-seqV2 bead-puck capture and library architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import surface_capture_rows

TITLE="Slide-seqV2"
NOTES="01_slide-seqv2.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41587-020-0739-1">Stickels et al., <i>Nature Biotechnology</i> (2021)</a>.'
SUMMARY="A lawn of randomly positioned 10-µm beads captures tissue mRNA. Each bead carries a 15-base split-pool barcode that is decoded in situ, an 8-base UMI and oligo-dT; template switching plus second-strand synthesis increases recovery before Nextera library preparation."
CAVEAT="The paper prints alternative bead versions. This schematic follows the published Chemgenes bead architecture and shows its split-pool barcode as a role-preserving placeholder rather than selecting one random barcode sequence."
SPATIAL=feature("spatial_barcode","spatial_barcode","combinatorial")
UMI_FEATURE=feature("umi","umi","random")

LIB, PRIMERS=truseq_library([
    seg("bead barcode", "J"*15, "cbc", placeholder=True, feature=SPATIAL),
    seg("bead constant", "TCTTCAGCGTTCCCGAGA"),
    seg("UMI", "N"*8, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("poly(dT) junction", "T"*12, placeholder=True),
    seg("captured cDNA", "X"*36, placeholder=True)],
    "Slide-seqV2 library", dual_index=False)
FINAL_LIBRARIES=(("Final Slide-seqV2 library",LIB,PRIMERS,
                  "The barcode-bearing read contains the bead barcode and UMI; the opposite read enters transcript-derived sequence.",
                  "The final duplex shows the Illumina primer sites used to read bead identity and transcript sequence."),)

def sections():
    return [
        ("Decode the randomly packed bead puck",
         [Row(chunks=[("10-µm bead lawn:  ● ● ● ● ●", "cbc", False)]),
          Row(chunks=[("in-situ ligation sequencing → each 15-base bead barcode gets an (x,y) coordinate", None, False)])],
         "Fourteen of the fifteen split-pool barcode bases are decoded on the puck and mapped to bead coordinates."),
        ("Capture tissue RNA on spatially decoded beads",
         surface_capture_rows("bead", "15-base bead barcode", "tissue mRNA"),
         "Oligo-dT captures polyadenylated RNA at the bead beneath the tissue; the same oligo contributes bead barcode and UMI."),
        ("Copy, add the distal handle and amplify",
         [Row(chunks=[("first-strand RT + template switch → tissue digestion → bead recovery", None, False)]),
          Row(chunks=[("→ second-strand synthesis → PCR → Nextera tagmentation", "me", False)])],
         "Slide-seqV2 adds bead-bound second-strand synthesis after RT, then amplifies and converts cDNA to an Illumina library."),
    ]
