"""Qin et al. 2016 TGIRT-seq, with the detailed 2021 protocol oligos."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]
import illumina as il
import seqprimers as sp
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow, feature,
                      revcomp)
from rna_prep import group_ii_starter_scene, rna_template, single_strand

TITLE="TGIRT-seq — group-II-intron RT template switching"
NOTES="01_tgirt-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1261/rna.054809.115">Qin et al., <i>RNA</i> (2016)</a>; exact current oligos: <a href="https://doi.org/10.21769/BioProtoc.4239">Xu et al., <i>Bio-protocol</i> (2021)</a>.'
SUMMARY="TGIRT-III switches from an R2 RNA/R2R DNA starter duplex onto a target RNA, copies through structured RNA, then a 5′-adenylated R1R oligo is ligated to cDNA before indexed PCR."

R2_RNA_WRITTEN="rArArGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrUrCrArC/3SpC3/"
R1R_WRITTEN="/5Phos/GATCGTCGGACTGTAGAACTCTGAACGTGTAG/3SpC3/"
R2_RNA="AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R2R_DNA="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTN"
R1R_DNA="GATCGTCGGACTGTAGAACTCTGAACGTGTAG"
R1_SITE=revcomp(R1R_DNA)
R2_SITE=revcomp(il.TRUSEQ_READ2)

INITIAL_NAME="target RNA pool"
INITIAL_ROWS=rna_template("target RNA",38)

def _library():
    lib=Construct([
      Segment("P5 partial",il.P5[:23],"p5"),
      Segment("i5","J"*6,"cbc",placeholder=True,feature=feature("sample_index_i5","sample_index","unknown")),
      Segment("Read 1 site",R1_SITE,"r1"),
      Segment("UMI","U"*6,"umi",placeholder=True,feature=feature("umi","umi","random")),
      Segment("RNA-derived insert","X"*38,placeholder=True),
      Segment("Index 1 / Read 2 arm",R2_SITE,"r2"),
      Segment("i7 reverse complement","I"*6,"cbc",placeholder=True,feature=feature("sample_index_i7","sample_index","unknown")),
      Segment("P7 reverse complement",il.P7_RC,"p7")],name="TGIRT-seq UMI library")
    primers=(
      sp.custom("Read 1","TGIRT/NEB small-RNA Read 1",R1_SITE,"Xu et al. Bio-protocol 2021"),
      sp.TRUSEQ["I1"],
      sp.custom("Index 2 (i5)","TGIRT i5 sequencing primer",revcomp(R1_SITE),"derived from the published P5 PCR arm"),
      sp.TRUSEQ["R2"],)
    problems=sp.verify(lib,primers)
    if problems: raise ValueError("TGIRT library: "+"; ".join(problems))
    return lib,primers

FINAL_LIBRARY,SEQ_PRIMERS=_library()
FINAL_CAPTION="The detailed protocol's optional-6N-UMI, six-base dual-index example after PCR; all adapter and primer sequences are printed by that protocol."
SEQUENCING_INTRO="This page selects the published 6N-UMI R1R option and the protocol's six-base TruSeq-barcode example on both PCR primers. Read 1 reports the UMI and then RNA-derived insert."
READ_LENGTHS={"Read 1":75,"Index 1 (i7)":6,"Index 2 (i5)":6,"Read 2":75}

def _starter():
    sc=group_ii_starter_scene(R2_RNA,R2R_DNA,rna_3_block="3SpC3")
    sc.strand("target RNA",[Segment("target RNA","X"*38,placeholder=True)],label="acceptor RNA")
    sc.labels("target RNA")
    return sc.rows()

def workflow():
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    w.react("Add blocked TGIRT starter duplex",_starter(),name="starter plus target RNA",note="The 3′-blocked R2 RNA pairs to R2R DNA, leaving a mixed one-nucleotide DNA overhang beside the target RNA.")
    w.react("Template switch and reverse transcribe",single_strand("cDNA",[Segment("R2R adapter",R2R_DNA[:-1],"r2"),Segment("RNA-derived cDNA","X"*38,placeholder=True)],"completed cDNA"),name="adapter-linked cDNA",note="TGIRT-III pairs the 1-nt overhang to the RNA 3′ nucleotide and copies the target at 60 °C.")
    w.react("Ligate 5′-App R1R to cDNA",single_strand("ligated cDNA",[Segment("R1R + optional UMI",R1R_DNA,"r1"),Segment("RNA-derived cDNA","X"*38,placeholder=True),Segment("R2R adapter",R2R_DNA[:-1],"r2")],"adapter-flanked cDNA"),name="two-adapter cDNA",note="Thermostable 5′ App DNA/RNA ligase joins 5′-adenylated, 3′-blocked R1R to the completed cDNA 3′ end.")
    w.react("Minimal indexed PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),name="final library",note="No more than 12 PCR cycles add P5, P7 and sample indices.")
    return w
