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
SOURCE=('Defining source: <a href="https://doi.org/10.1038/nmeth.3370">Macaulay et al., <i>Nature Methods</i> (2015)</a>; '
        '<a href="https://www.researchgate.net/profile/Mubeen-Goolam/publication/275525502_GT-seq_Parallel_sequencing_of_single-cell_genomes_and_transcriptomes/links/5757e80f08ae5c65490724bc/G-T-seq-Parallel-sequencing-of-single-cell-genomes-and-transcriptomes.pdf">full original Online Methods</a>; '
        '<a href="https://www.pure.ed.ac.uk/ws/files/181359881/Separation_and_parallel_sequencing_of_the_genomes_and_transcriptomes_of_single_cells_using_G_T_seq.pdf">authors\' detailed protocol, <i>Nature Protocols</i> (2016)</a>.')
SUMMARY="A biotinylated oligo(dT) captures a cell's poly(A) RNA on streptavidin beads, physically separating it from genomic DNA; the two fractions then undergo Smart-seq2-like WTA and MDA/PicoPLEX WGA, respectively."
CAVEAT="G&T-seq separates one cell into genome and transcriptome branches. The original 2015 Online Methods print the capture, template-switch and PCR oligos. Their TSO is AAGCAGTGGTATCAACGCAGAGTACrGrG+G (Online Methods lines 597–599); the authors' 2016 detailed protocol prints AAGCAGTGGTATCAACGCAGAGTACATrGrG+G (oligo list lines 652–653), adding AT immediately before rGrG+G. The model uses the 2015 value for this defining version; whether the 2016 value is an intentional revision or a source error is unresolved. Downstream selected Nextera XT index members and WGA kit-specific oligos still need source checks."
CAPTURE_HANDLE="AAGCAGTGGTATCAACGCAGAGTAC"
CAPTURE_POLYT="T"*30
TSO_DNA="AAGCAGTGGTATCAACGCAGAGTAC"
TSO_WRITTEN=TSO_DNA+"rGrG+G"
TSO_2016="AAGCAGTGGTATCAACGCAGAGTACATrGrG+G"
ISPCR="AAGCAGTGGTATCAACGCAGAGT"
SOURCE_OLIGO_EXCEPTIONS={
    "tso": "The original 2015 Online Methods (author-uploaded PDF, lines 597–599) print the shorter 2015 TSO. The locally downloadable 2016 author protocol instead prints an extra AT before rGrG+G; the 2015 PDF is web-readable but direct download is blocked in this environment."
}
RNA,RNA_P=nextera_library([seg("full-length cDNA","X"*44,placeholder=True)],"G&T-seq transcriptome library")
DNA,DNA_P=nextera_library([seg("WGA genomic fragment","X"*44,placeholder=True)],"G&T-seq genome library")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RNA_P,"Bead-captured mRNA undergoes Smart-seq2 reverse transcription/WTA followed by Nextera library construction.","Standard Nextera primers read the amplified cDNA insert from both ends."),("Genome library",DNA,DNA_P,"The mRNA-depleted supernatant undergoes MDA or PicoPLEX WGA followed by Nextera library construction.","Standard Nextera primers read the WGA-derived genomic insert from both ends."))

def capture_scene():
    rna=[Segment("mRNA body","X"*28,placeholder=True),Segment("poly(A)","A"*30)]
    oligo=[Segment("biotin capture handle",CAPTURE_HANDLE,"w1"),
           Segment("poly(dT)",CAPTURE_POLYT),Segment("3′ VN anchor","VN",placeholder=True)]
    sc=Scene(); sc.strand("mRNA",rna,label="cellular mRNA")
    sc.anneal("capture oligo",oligo,to="mRNA",pair=("poly(dT)","poly(A)"),
              label="biotin-oligo(dT)",unpaired=("3′ VN anchor",))
    sc.mark("capture oligo","biotin capture handle","streptavidin bead capture")
    return sc

def oligos():
    return [
        "Original 2015 Online Methods, lines 573–605:",
        f"Biotinylated oligo-dT30VN: 5′-Biotin-TEG-{CAPTURE_HANDLE}{CAPTURE_POLYT}VN-3′",
        f"Template switching oligo: 5′-{TSO_WRITTEN}-3′",
        f"PCR primer: 5′-{ISPCR}-3′",
        f"2016 author protocol TSO: 5′-{TSO_2016}-3′ (extra AT before rGrG+G; model uses 2015 TSO; reason for change unresolved)",
    ]

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
 capture=g.input("capture_oligo","biotinylated oligo(dT) capture oligo",MaterialKind.OLIGO,
                 state=(Segment("5′ biotin-TEG capture handle",CAPTURE_HANDLE),
                        Segment("dT30",CAPTURE_POLYT),Segment("VN anchor","VN",placeholder=True)),
                 source=SOURCE+" (2015 Online Methods lines 573–574)")
 hybrid=g.product("hybrid","lysate with poly(A) RNA hybridized to biotin-oligo(dT)",MaterialKind.MOLECULAR,state=capture_scene())
 g.react("hybridize","Hybridize poly(A) RNA to biotinylated oligo(dT)",StepKind.HYBRIDIZATION,((lysate,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL)),(hybrid,))
 beads=g.input("streptavidin","streptavidin magnetic beads",MaterialKind.SELECTION)
 rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=molecule("RNA fraction",(role("poly(A) RNA"),role("biotin-oligo(dT)"))))
 dna=g.product("dna_fraction","mRNA-depleted genomic-DNA supernatant",MaterialKind.MOLECULAR,state=molecule("DNA fraction",(role("genomic DNA"),)))
 g.react("separate","Physically separate RNA and genomic-DNA fractions",StepKind.AFFINITY_CAPTURE,((hybrid,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna))
 rt=g.input("rt","template-switching reverse transcriptase",MaterialKind.CATALYST)
 first=g.product("first_strand","oligo(dT)-primed transcript cDNA",MaterialKind.MOLECULAR,state=rna.state)
 g.react("reverse_transcribe","Reverse-transcribe bead-captured mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 tso=g.input("tso","2015 G&T-seq template-switching oligo",MaterialKind.OLIGO,
             state=(Segment("TSO DNA handle",TSO_DNA),Segment("3′ rGrG+G","GGG")),
             source=SOURCE+" (2015 Online Methods lines 597–599)")
 switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=template_switch_scene(inferred=True))
 g.react("template_switch","Template-switch at the transcript 5-prime end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 rf=g.input("rna_wta_forward","2015 G&T-seq PCR primer, forward landing",MaterialKind.OLIGO,
            state=(Segment("G&T universal PCR primer",ISPCR),),
            source=SOURCE+" (2015 Online Methods lines 604–605)")
 rr=g.input("rna_wta_reverse","2015 G&T-seq PCR primer, opposite landing",MaterialKind.OLIGO,
            state=(Segment("G&T universal PCR primer",ISPCR),),
            source=SOURCE+" (2015 Online Methods lines 604–605)")
 rp=g.input("rna_wta_polymerase","RNA-WTA polymerase",MaterialKind.CATALYST)
 rwta=g.product("rna_wta","whole-transcriptome-amplified cDNA",MaterialKind.MOLECULAR,state=switched.state)
 g.react("amplify_rna","Amplify the transcriptome",StepKind.PCR,((switched,InputRole.TEMPLATE),(rf,InputRole.STRUCTURAL),(rr,InputRole.STRUCTURAL),(rp,InputRole.CATALYST)),(rwta,),geometry=PCRGeometry.unavailable(rwta.state,reason="The exact primer is printed in the 2015 Methods; this generic Smart-seq2-like cDNA state does not place its two binding sites."))
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
