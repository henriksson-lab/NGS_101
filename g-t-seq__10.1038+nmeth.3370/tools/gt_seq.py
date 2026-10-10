"""Macaulay et al. G&T-seq separation and branch model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg
from chemdraw import Row, Scene, Segment
from multimodal_spatial import modality_split
from rna_special import template_switch_scene

TITLE="G&T-seq — separate genome and transcriptome from one cell"
NOTES="01_g-t-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nmeth.3370">Macaulay et al., <i>Nature Methods</i> (2015)</a>.'
SUMMARY="A biotinylated oligo(dT) captures a cell's poly(A) RNA on streptavidin beads, physically separating it from genomic DNA; the two fractions then undergo Smart-seq2-like WTA and MDA/PicoPLEX WGA, respectively."
CAVEAT="G&T-seq is the separation protocol plus two downstream amplification branches. The page does not imply that one molecular construct contains both RNA and DNA."
RNA,RNA_P=nextera_library([seg("full-length cDNA","X"*44,placeholder=True)],"G&T-seq transcriptome library")
DNA,DNA_P=nextera_library([seg("WGA genomic fragment","X"*44,placeholder=True)],"G&T-seq genome library")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RNA_P,"Bead-captured mRNA undergoes Smart-seq2 reverse transcription/WTA followed by Nextera library construction.","Standard Nextera primers read the amplified cDNA insert from both ends."),("Genome library",DNA,DNA_P,"The mRNA-depleted supernatant undergoes MDA or PicoPLEX WGA followed by Nextera library construction.","Standard Nextera primers read the WGA-derived genomic insert from both ends."))

def capture_scene():
    rna=[Segment("mRNA body","X"*28,placeholder=True),Segment("poly(A)","A"*16)]
    oligo=[Segment("biotin capture handle","X"*10,"w1",placeholder=True),Segment("poly(dT)","T"*16)]
    sc=Scene(); sc.strand("mRNA",rna,label="cellular mRNA")
    sc.anneal("capture oligo",oligo,to="mRNA",pair=("poly(dT)","poly(A)"),label="biotin-oligo(dT)")
    sc.mark("capture oligo","biotin capture handle","streptavidin bead capture")
    return sc

def sections(): return [
 ("Lyse one cell and hybridize biotinylated oligo(dT)",capture_scene().rows(),"The capture oligo binds poly(A) RNA while genomic DNA remains unbound."),
 ("Immobilize RNA and separate the fractions",modality_split("bead-bound mRNA","gDNA-containing supernatant"),"Magnetic separation creates physically distinct transcriptome and genome inputs from the same cell."),
 ("Amplify the transcriptome",template_switch_scene(inferred=True).rows(),"The RNA fraction follows a Smart-seq2-like oligo(dT) RT, template-switch and whole-transcriptome amplification."),
 ("Amplify the genome",[Row(chunks=[("single-cell gDNA -> MDA or PicoPLEX whole-genome amplification",None,False)])],"The paper validated both MDA and PicoPLEX; these are alternative branches, not consecutive reactions."),
 ("Construct separate Nextera libraries",modality_split("tagmented transcriptome library","tagmented genome library"),"The two amplified products remain separate through indexing and sequencing."),
]
