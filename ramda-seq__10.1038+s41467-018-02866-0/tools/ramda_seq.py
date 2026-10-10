"""Hayashi et al. 2018 RamDA-seq chemistry."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]
from batch_ngs import nextera_library, seg
from chemdraw import MolecularState, Scene, Workflow
from rna_prep import rna_template, rna_dna_hybrid, single_strand

TITLE="RamDA-seq — single-cell full-length total RNA"
NOTES="01_ramda-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-018-02866-0">Hayashi et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY="Not-so-random and oligo-dT primers initiate RT; DNase I nicks cDNA while RNase-H-minus RT displaces and amplifies it, then second-strand synthesis and Nextera XT make the sequencing library."
INITIAL_NAME="single-cell total RNA"
INITIAL_ROWS=rna_template("total RNA",38)
FINAL_LIBRARY, SEQ_PRIMERS=nextera_library([seg("RamDA cDNA insert","X"*38,placeholder=True)],"RamDA-seq library")
FINAL_CAPTION="Nextera XT library made from double-stranded RT-RamDA product."
SEQUENCING_INTRO="The defining study used single-read 76-cycle NextSeq sequencing; the page also locates the complete standard Nextera primer set present on the finished library."
READ_LENGTHS={"Read 1":76}

def workflow():
    amplified=single_strand("displaced cDNA",[seg("gp32-protected displaced cDNA","X"*46,placeholder=True)],"displaced cDNA products")
    ds=Scene.duplex([seg("amplified cDNA","X"*42,placeholder=True)],label="double-stranded cDNA").rows()
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    for action,rows,note in [
      ("NSR + oligo-dT reverse transcription",rna_dna_hybrid(length=38),"A pool of 408 first-strand NSR hexamers avoids exact rRNA matches while oligo(dT) captures poly(A) RNA."),
      ("DNase-I nicking and RT displacement",amplified,"DNase I nicks cDNA; RNase-H-minus PrimeScript extends from nicked 3' ends, while T4 gp32 promotes displacement and protects products."),
      ("Second-strand NSR synthesis",ds,"Complementary second-strand NSRs prime Klenow Fragment (3'→5' exo−)."),
      ("Nextera XT tagmentation and PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),"The study used scaled Nextera XT and 13 PCR cycles (14 in its advanced method)."),
    ]: w.react(action,rows,note=note)
    return w
