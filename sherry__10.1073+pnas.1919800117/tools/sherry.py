"""Di et al. 2020 SHERRY chemistry."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]

from batch_ngs import nextera_library, seg
from chemdraw import MolecularState, Scene, Workflow
from rna_prep import rna_template, rna_dna_hybrid
import nextera as nx

TITLE="SHERRY — direct tagmentation of RNA/DNA hybrids"
NOTES="01_sherry.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1919800117">Di et al., <i>PNAS</i> (2020)</a>.'
SUMMARY="Oligo-dT reverse transcription makes an RNA/cDNA hybrid; Tn5 directly fragments and tags that hybrid before gap extension and indexed PCR."

INITIAL_NAME="poly(A) RNA"
INITIAL_ROWS=rna_template("poly(A) RNA", 36)
FINAL_LIBRARY, SEQ_PRIMERS=nextera_library([seg("transcript-derived insert","X"*38,placeholder=True)],"SHERRY library")
FINAL_CAPTION="The amplifiable S5–insert–S7 product after gap extension and indexed PCR."
SEQUENCING_INTRO="The paper used Illumina NextSeq 500 or HiSeq 4000; standard Nextera primer sites are installed by the published transposomes and indexed PCR primers."
def _tagged():
    sc=Scene()
    sc.strand("tagged RNA",[seg("transferred S5–ME",nx.ADAPTOR_S5,"s5+me"),seg("RNA fragment","X"*30,placeholder=True)],label="tagged RNA product")
    sc.strand("tagged cDNA",[seg("transferred S7–ME",nx.ADAPTOR_S7,"s7+me"),seg("cDNA fragment","X"*30,placeholder=True)],label="tagged cDNA product")
    sc.labels("tagged RNA"); sc.labels("tagged cDNA")
    return sc.rows()

def workflow():
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    for action,rows,note in [
      ("Oligo-dT reverse transcription",rna_dna_hybrid(length=36),"SuperScript II copies poly(A) RNA; no second-strand synthesis is performed."),
      ("Tn5 tagment RNA/DNA hybrid",_tagged(),"Strand tests demonstrated adapter transfer to both fragmented RNA and cDNA; these are transferred-strand products, not a gap-filled duplex."),
      ("Gap extension and indexed PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),"Polymerase extension makes the tagmented cDNA amplifiable; indexed common primers complete P5 and P7."),
    ]: w.react(action,rows,note=note)
    return w
