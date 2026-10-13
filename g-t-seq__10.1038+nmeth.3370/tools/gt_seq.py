"""Macaulay et al. G&T-seq separation and branch model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg
from chemdraw import Row, Scene, Segment
from multimodal_spatial import finish_pcr_branch, modality_split
from pcr import PCRGeometry
from rna_special import template_switch_scene
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

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

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 lysate=g.input("lysate","single-cell lysate containing genomic DNA and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("single-cell lysate",(role("genomic DNA"),role("poly(A) RNA"))),primary=True)
 capture=g.input("capture_oligo","biotinylated oligo(dT) capture oligo",MaterialKind.OLIGO)
 hybrid=g.product("hybrid","lysate with poly(A) RNA hybridized to biotin-oligo(dT)",MaterialKind.MOLECULAR,state=capture_scene())
 g.react("hybridize","Hybridize poly(A) RNA to biotinylated oligo(dT)",StepKind.HYBRIDIZATION,((lysate,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL)),(hybrid,))
 beads=g.input("streptavidin","streptavidin magnetic beads",MaterialKind.SELECTION)
 rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=molecule("RNA fraction",(role("poly(A) RNA"),role("biotin-oligo(dT)"))))
 dna=g.product("dna_fraction","mRNA-depleted genomic-DNA supernatant",MaterialKind.MOLECULAR,state=molecule("DNA fraction",(role("genomic DNA"),)))
 g.react("separate","Physically separate RNA and genomic-DNA fractions",StepKind.AFFINITY_CAPTURE,((hybrid,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna))
 rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
 first=g.product("first_strand","oligo(dT)-primed transcript cDNA",MaterialKind.MOLECULAR,state=rna.state)
 g.react("reverse_transcribe","Reverse-transcribe bead-captured mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 tso=g.input("tso","Smart-seq2-like template-switching oligo",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (G&T-seq permits the Smart-seq2 branch)")
 switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=template_switch_scene(inferred=True))
 g.react("template_switch","Template-switch at the transcript 5-prime end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 rf=g.input("rna_wta_forward","RNA-WTA forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (downstream Smart-seq2 configuration)")
 rr=g.input("rna_wta_reverse","RNA-WTA reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (downstream Smart-seq2 configuration)")
 rp=g.input("rna_wta_polymerase","RNA-WTA polymerase",MaterialKind.CATALYST)
 rwta=g.product("rna_wta","whole-transcriptome-amplified cDNA",MaterialKind.MOLECULAR,state=switched.state)
 g.react("amplify_rna","Amplify the transcriptome",StepKind.PCR,((switched,InputRole.TEMPLATE),(rf,InputRole.STRUCTURAL),(rr,InputRole.STRUCTURAL),(rp,InputRole.CATALYST)),(rwta,),geometry=PCRGeometry.unavailable(rwta.state,reason="The paper permits a Smart-seq2-like downstream configuration."))
 random=g.input("wga_primer","MDA or PicoPLEX WGA priming reagent",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (alternative whole-genome-amplification branches)")
 wga_pol=g.input("wga_polymerase","MDA or PicoPLEX amplification polymerase",MaterialKind.CATALYST)
 wga=g.product("wga","whole-genome-amplified DNA",MaterialKind.MOLECULAR,state=dna.state)
 g.react("amplify_dna","Amplify the genome by MDA or PicoPLEX",StepKind.EXTENSION,((dna,InputRole.TEMPLATE),(random,InputRole.STRUCTURAL),(wga_pol,InputRole.CATALYST)),(wga,))
 reviewed=["lysate","capture_oligo","hybridize","streptavidin","separate","rt","reverse_transcribe","tso","template_switch","rna_wta_forward","rna_wta_reverse","rna_wta_polymerase","amplify_rna","wga_primer","wga_polymerase","amplify_dna"]
 branches=(("rna",rwta,RNA,RNA_P,"transcriptome"),("dna",wga,DNA,DNA_P,"genome"))
 for prefix,template,final,primers,label in branches:
  transposome=g.input(prefix+"_transposome",label+" loaded Nextera transposome",MaterialKind.ADAPTER)
  tn5=g.input(prefix+"_tn5",label+" Tn5 transposase",MaterialKind.CATALYST)
  tagged=g.product(prefix+"_tagged",label+" tagmented product",MaterialKind.MOLECULAR,state=final)
  waste=g.product(prefix+"_tn5_waste",label+" non-library tagmentation products",MaterialKind.WASTE)
  g.react(prefix+"_tagment","Tagment the amplified "+label,StepKind.TAGMENTATION,((template,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste)); g.discard(waste,"does not retain an amplifiable Nextera fragment")
  reviewed.extend((transposome.id,tn5.id,prefix+"_tagment"))
  finish_pcr_branch(g,prefix,tagged,final,primers,reviewed,label="G&T-seq "+label)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
