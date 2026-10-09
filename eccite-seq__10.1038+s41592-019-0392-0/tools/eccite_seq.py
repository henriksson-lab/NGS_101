"""ECCITE-seq 5-prime multimodal tag and direct-guide capture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row, feature
from multimodal_spatial import modality_split

TITLE="ECCITE-seq"
NOTES="01_eccite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-019-0392-0">Mimitou et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY="ECCITE-seq adapts 10x 5′ chemistry so transcript, antibody tag, sample hashtag and directly captured sgRNA products receive a common cell barcode and UMI. Modality-specific handles then support separate libraries, alongside standard immune-receptor enrichment when desired."
CAVEAT="The schematic shows the three construct classes that change library chemistry (transcript, antibody/hashtag tag and guide-derived oligo). Immune-receptor enrichment is the unmodified 10x V(D)J branch and is linked in the notes rather than redrawn as a new chemistry."
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-5prime")
UMI=feature("umi","umi","random")
TAG_FEATURE=feature("antibody_or_hashtag_barcode","feature_barcode","unknown")
GUIDE_FEATURE=feature("guide","guide_barcode","unknown")

def library(name,feature):
    return truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),*feature],name,dual_index=False)
GEX,GEX_P=library("ECCITE-seq transcript library",[seg("5-prime cDNA","X"*36,placeholder=True)])
ADT,ADT_P=library("ECCITE-seq ADT/hashtag library",[seg("tag-specific handle","H"*18,placeholder=True),seg("antibody or hashtag barcode","A"*15,"cbc",placeholder=True,feature=TAG_FEATURE)])
GDO,GDO_P=library("ECCITE-seq guide-derived oligo library",[seg("sgRNA protospacer","G"*20,"cbc",placeholder=True,feature=GUIDE_FEATURE),seg("guide scaffold capture region","S"*18,placeholder=True)])
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"5′ transcript sequence follows cell barcode and UMI.","The transcript library's primer geometry is computed on the finished duplex."),
    ("Final antibody / hashtag library",ADT,ADT_P,"Tag-specific PCR separates protein and sample-hashtag products.","The tag library's primer geometry is computed on the finished duplex."),
    ("Final guide-derived oligo library",GDO,GDO_P,"Direct guide capture reports sgRNA identity with the same cell barcode and UMI.","The guide library's primer geometry is computed on the finished duplex."),)

def sections():
    return [
        ("Prime tags and sgRNAs during 5′ reverse transcription",
         [Row(chunks=[("antibody/hashtag oligo ⇄ bead-TSO-complementary end → cell barcode + UMI", "cbc", False)]),
          Row(chunks=[("sgRNA protospacer—scaffold ⇄ scaffold RT primer → template switch → cell barcode + UMI", "umi", False)])],
         "An antibody tag anneals to the bead-associated template-switch system; a separate RT primer binds the sgRNA scaffold, copies the guide and then template-switches."),
        ("Amplify all barcoded products, then enrich by handle",
         modality_split("5′ gene-expression / optional V(D)J","antibody-derived tags","sample hashtags","guide-derived oligos"),
         "A shared cell identity is installed before modality-specific enrichment. Different handles keep the short tag libraries distinct."),]
