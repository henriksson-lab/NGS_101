"""Spatial-CUT&Tag molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows,spatial_barcoded_library
from multimodal_spatial import barcode_rounds
FINAL_LIBRARY,SEQ_PRIMERS=spatial_barcoded_library("CUT&Tag chromatin fragment")
TITLE="Spatial-CUT&Tag — orthogonal in-tissue barcoding of chromatin marks"
NOTES="01_spatial-cut-tag.html"
SOURCE='<a href="https://doi.org/10.1126/science.abg7216">Deng et al., <i>Science</i> (2022)</a>.'
SUMMARY="Antibody-tethered pA–Tn5 tags a histone mark in a fixed tissue section. Two perpendicular microfluidic flows ligate barcode A and then barcode B, so their ordered pair defines a tissue pixel."
CAVEAT="INFERRED — the defining article establishes two in-tissue ligation rounds, but the accessible paper text does not print complete sequencing-adapter oligos. Outer TruSeq geometry is shown as an inferred shell."
FINAL_CAPTION="Representative doubly barcoded CUT&Tag fragment. Barcode A and barcode B form one combinatorial spatial-pixel identity."
SEQUENCING_INTRO="Paired reads map the CUT&Tag fragment; the two ligated barcode blocks identify its row-by-column tissue pixel."
def sections(): return [
 ("Bind antibody and pA–Tn5 in the section",antibody_enzyme_rows("Tn5","targeted tagmentation"),"CUT&Tag installs adapter-bearing ends at the chosen histone modification."),
 ("Flow barcode A through parallel channels",barcode_rounds(1,label="spatial barcode A"),"Each channel ligates one A identity in situ."),
 ("Rotate the chip and flow barcode B",barcode_rounds(1,label="spatial barcode B"),"The orthogonal B channel crosses A channels; the A×B pair defines a pixel."),
 ("Reverse crosslinks and PCR",[Row(chunks=[("P5 — A — B — targeted genomic insert — P7",None,False)])],"Released fragments are PCR-completed for paired-end sequencing."),]
