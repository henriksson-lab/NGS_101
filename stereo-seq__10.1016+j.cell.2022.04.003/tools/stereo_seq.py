"""Stereo-seq DNA-nanoball patterned-array spatial transcriptomics."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, circle_rows, feature
from circular import circularize_ssdna
from multimodal_spatial import surface_capture_rows

TITLE="Stereo-seq"
NOTES="01_stereo-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2022.04.003">Chen et al., <i>Cell</i> (2022)</a>.'
SUMMARY="Spatially indexed DNA nanoballs are arrayed and decoded before receiving a UMI and oligo-dT capture sequence. Tissue mRNA is reverse-transcribed on this surface; released cDNA retains the coordinate identifier and is converted to a DNBSEQ library."
CAVEAT="The paper establishes coordinate-identifier and molecular-identifier roles, but complete STOmics chip and sequencing-primer bases are commercial. They are shown as inferred placeholders and no proprietary sequence is invented."
SPATIAL=feature("coordinate_identifier","spatial_barcode","unknown")
UMI_FEATURE=feature("molecular_identifier","umi","random")
SEQ_PRIMERS=()
FINAL_LIBRARY=None
SEQUENCING_ENDING="The completed library is circularized and amplified into sequencing DNA nanoballs. A DNBSEQ/cPAS primer anneals to the repeated platform-adapter site immediately before the coordinate-identifier/UMI-bearing insert. The current primer and adapter bases are proprietary, so the final panel shows the binding geometry and extension direction rather than a fabricated sequence."

LINEAR=Construct([seg("DNBSEQ primer site","P"*18,"r1",placeholder=True,inferred=True),seg("coordinate identifier","C"*25,"cbc",placeholder=True,feature=SPATIAL),seg("molecular identifier","M"*10,"umi",placeholder=True,feature=UMI_FEATURE),seg("captured cDNA","X"*38,placeholder=True),seg("distal DNBSEQ adapter","Q"*18,placeholder=True,inferred=True)],name="Stereo-seq ssDNA library")
CIRCLE=circularize_ssdna(LINEAR,five_prime_phosphate=True)

def sections():
    return [
        ("Pattern and decode coordinate-identifier DNA nanoballs",
         [Row(chunks=[("DNB array:  ●CID₁  ●CID₂  ●CID₃  …", "cbc", False)]),
          Row(chunks=[("sequencing-by-synthesis before tissue placement → CID-to-(x,y) map", None, False)])],
         "Each surface DNB carries many copies of one coordinate identifier; pre-decoding maps molecular barcode to a physical location."),
        ("Add UMI/oligo-dT and capture tissue RNA",
         surface_capture_rows("decoded DNB", "coordinate identifier", "tissue mRNA"),
         "A UMI-bearing poly(dT) capture sequence is ligated onto the decoded chip oligo before tissue RNA capture."),
        ("Reverse-transcribe, release and amplify spatial cDNA",
         [Row(chunks=[("surface CID—MID—poly(dT) || RNA → reverse transcription + template switch", "umi", False)]),
          Row(chunks=[("→ release cDNA → PCR → fragmentation / platform-library conversion", None, False)])],
         "The coordinate identifier (CID) and molecular identifier (MID) remain attached to the cDNA through bulk library preparation."),
        ("Circularize the sequencing strand and form library DNBs",
         circle_rows(CIRCLE.linear,closure_label="splint ligation"),
         "INFERRED — A platform library strand is splint-ligated into a circle; undisclosed platform-adapter bases are visibly inferred."),
        ("Bind the sequencing primer at every repeat",
         [Row(chunks=[("…[primer site]—CID—MID—cDNA—[adapter][primer site]—CID—MID—cDNA…", "r1", False)]),
          Row(chunks=[("   sequencing primer →                    sequencing primer →", "r2", False)])],
         "INFERRED — The undisclosed cPAS primer binds the repeated adapter role and extends toward the coordinate/UMI-bearing insert."),]
