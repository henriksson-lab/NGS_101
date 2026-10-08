"""Original CITE-seq with 10x 3-prime transcriptome and ADT readout."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from multimodal_spatial import modality_split

TITLE="CITE-seq"
NOTES="01_cite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.4380">Stoeckius et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY="Antibodies carry poly(A)-tailed DNA tags containing an antibody identity and an amplification handle. The same oligo-dT bead primer that captures mRNA also copies these antibody-derived tags, transferring the droplet cell barcode and UMI before ADT and transcript libraries are separated by size."
CAVEAT="The original paper demonstrates both Drop-seq and 10x implementations. This page follows its 10x 3′ implementation; antibody barcode identities are placeholders, not a claim that all panels use one fixed sequence."

GEX,GEX_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("poly(dT) junction","T"*10,placeholder=True),seg("mRNA cDNA","X"*34,placeholder=True)],"CITE-seq transcript library",dual_index=False)
ADT,ADT_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("ADT PCR handle","H"*18,placeholder=True),seg("antibody barcode","A"*15,"cbc",placeholder=True)],"CITE-seq ADT library",dual_index=False)
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"The 10x cell barcode and UMI identify each transcript molecule.","The final transcript duplex shows the declared Illumina primer geometry."),
    ("Final antibody-derived-tag library",ADT,ADT_P,"The same cell barcode and UMI are followed by the antibody-specific tag.","The final ADT duplex shows the declared Illumina primer geometry."),)

def sections():
    return [
        ("Label surface proteins with poly(A)-tailed DNA tags",
         [Row(chunks=[("antibody—[ADT PCR handle][antibody barcode][poly(A)]", "cbc", False)])],
         "Each antibody species is associated with a distinct DNA barcode; the poly(A) tail makes the tag compatible with an oligo-dT single-cell RNA workflow."),
        ("Copy mRNA and antibody tags in the same droplet",
         [Row(chunks=[("bead [cell barcode][UMI]—oligo(dT) || mRNA poly(A)", "umi", False)]),
          Row(chunks=[("bead [cell barcode][UMI]—oligo(dT) || ADT poly(A)", "cbc", False)])],
         "Reverse transcriptase copies RNA and extends on the DNA ADT, transferring the same droplet barcode system to both modalities."),
        ("Separate by size and amplify independently",modality_split("large transcript-derived cDNA","short antibody-derived tags"),
         "SPRI size separation sends mRNA cDNA through the standard gene-expression branch and ADTs through tag-specific PCR."),]
