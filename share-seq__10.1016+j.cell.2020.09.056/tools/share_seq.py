"""SHARE-seq joint RNA/accessibility combinatorial indexing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Row, feature
from multimodal_spatial import barcode_rounds, modality_split

TITLE="SHARE-seq"
NOTES="01_share-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2020.09.056">Ma et al., <i>Cell</i> (2020)</a>.'
SUMMARY="Fixed permeabilized cells are first tagmented and reverse-transcribed. Three split-pool ligation rounds then append the same barcode combination to ATAC fragments and biotinylated cDNA before the two modalities are separated and built as independent libraries."
CAVEAT="The paper supplies plate oligos in Supplementary Table S1. This page models their verified three-round ligation roles; individual well-barcode bases remain placeholders."
BC1=feature("cell_barcode_round1","cell_barcode","combinatorial",group="cell_barcode",part="round 1")
BC2=feature("cell_barcode_round2","cell_barcode","combinatorial",group="cell_barcode",part="round 2")
BC3=feature("cell_barcode_round3","cell_barcode","combinatorial",group="cell_barcode",part="round 3")
UMI_FEATURE=feature("umi","umi","random")

RNA,RNA_P=truseq_library([seg("barcode round 3","C"*8,"cbc",placeholder=True,feature=BC3),seg("barcode round 2","B"*8,"cbc",placeholder=True,feature=BC2),seg("barcode round 1","A"*8,"cbc",placeholder=True,feature=BC1),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI_FEATURE),seg("cDNA","X"*34,placeholder=True)],"SHARE-seq RNA library")
ATAC,ATAC_P=nextera_library([seg("barcode round 3","C"*8,"cbc",placeholder=True,feature=BC3),seg("barcode round 2","B"*8,"cbc",placeholder=True,feature=BC2),seg("barcode round 1","A"*8,"cbc",placeholder=True,feature=BC1),seg("accessible DNA","X"*34,placeholder=True)],"SHARE-seq ATAC library")
FINAL_LIBRARIES=(
    ("Final RNA library",RNA,RNA_P,"Three ligated well barcodes and the RT UMI identify each transcript and cell.","The RNA library's declared primer sites are computed on the completed duplex."),
    ("Final ATAC library",ATAC,ATAC_P,"The same three well barcodes identify the cell of origin of each accessible-DNA fragment.","The ATAC library's declared primer sites are computed on the completed duplex."),)

def sections():
    return [
        ("Tag chromatin and reverse-transcribe RNA in the same fixed cell",
         [Row(chunks=[("open chromatin → Tn5-tagged DNA", "me", False)]),Row(chunks=[("poly(A) RNA → biotin-UMI oligo-dT → cDNA", "umi", False)])],
         "The two molecular classes acquire modality-specific handles before cellular indexing."),
        ("Install three shared split-pool barcodes by ligation",barcode_rounds(3),
         "Cells are redistributed through three 96-well plates. Hybridization aligns each well barcode, and ligase covalently joins it to both RNA-derived and ATAC-derived molecules."),
        ("Release cells and separate modalities",modality_split("streptavidin-captured RNA-derived cDNA","ATAC-derived genomic DNA"),
         "Reverse crosslinking releases molecules; biotin enables physical separation before modality-specific amplification."),]
