"""Macaulay et al. G&T-seq separation and branch model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library, seg
from chemdraw import Scene, Segment
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="G&T-seq — separate genome and transcriptome from one cell"
NOTES="01_g-t-seq.html"
SOURCE=('Defining source: <a href="https://doi.org/10.1038/nmeth.3370">Macaulay et al., <i>Nature Methods</i> (2015)</a>; '
        '<a href="https://zenodo.org/records/895966">accepted author manuscript with Online Methods (Zenodo 895966)</a>; '
        '<a href="https://www.pure.ed.ac.uk/ws/files/181359881/Separation_and_parallel_sequencing_of_the_genomes_and_transcriptomes_of_single_cells_using_G_T_seq.pdf">authors\' detailed protocol, <i>Nature Protocols</i> (2016)</a>.')
SUMMARY="A biotinylated oligo(dT) captures a cell's poly(A) RNA on streptavidin beads, physically separating it from genomic DNA; the two fractions then undergo Smart-seq2-like WTA and MDA/PicoPLEX WGA, respectively."
CAVEAT="G&T-seq separates one cell into genome and transcriptome branches. The original 2015 Online Methods print the capture, template-switch and PCR oligos. Their TSO is AAGCAGTGGTATCAACGCAGAGTACrGrG+G (accepted-manuscript Online Methods, extracted lines 336–337); the authors' 2016 detailed protocol prints AAGCAGTGGTATCAACGCAGAGTACATrGrG+G (oligo list lines 652–653), adding AT immediately before rGrG+G. The model uses the 2015 value for this defining version; whether the 2016 value is an intentional revision or a source error is unresolved. Downstream selected Nextera XT index members and WGA kit-specific oligos still need source checks. INFERRED — the template-switch panel draws untemplated C additions pairing with the TSO rGrG+G end; that mechanism and its length are not printed. The drawn capture oligo and TSO residues are the printed 2015 oligos."
CAPTURE_HANDLE="AAGCAGTGGTATCAACGCAGAGTAC"
CAPTURE_POLYT="T"*30
TSO_DNA="AAGCAGTGGTATCAACGCAGAGTAC"
TSO_WRITTEN=TSO_DNA+"rGrG+G"
TSO_2016="AAGCAGTGGTATCAACGCAGAGTACATrGrG+G"
ISPCR="AAGCAGTGGTATCAACGCAGAGT"
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

def template_switch_local_scene():
    """Bead-captured RT and template switch drawn from the printed G&T oligos.

    The first strand is primed by the capture oligo itself. The untemplated C
    additions that pair with the TSO rGrG+G end are a mechanistic inference; their
    number is not printed by the source and the segment is marked inferred.
    """
    rna=[Segment("mRNA body","X"*30,placeholder=True),Segment("poly(A)","A"*30)]
    first=[Segment("capture handle",CAPTURE_HANDLE,"w1"),Segment("poly(dT)",CAPTURE_POLYT),
           Segment("3′ VN anchor","VN",placeholder=True),
           Segment("reverse-transcribed cDNA","X"*28,placeholder=True),
           Segment("untemplated C additions","CCC",inferred=True)]
    sc=Scene(); sc.strand("mRNA",rna,label="bead-captured mRNA")
    sc.anneal("first strand",first,to="mRNA",pair=("poly(dT)","poly(A)"),
              label="capture-oligo-primed cDNA",
              unpaired=("3′ VN anchor","reverse-transcribed cDNA","untemplated C additions"))
    sc.mark("first strand","capture handle","5′ biotin-TEG capture/RT oligo")
    sc.anneal("TSO",[Segment("TSO DNA handle",TSO_DNA,"w1"),Segment("3′ rGrG+G","GGG")],
              to="first strand",pair=("3′ rGrG+G","untemplated C additions"),label="2015 TSO")
    sc.mark("TSO","3′ rGrG+G","rGrG+G: two RNA G then one LNA G (printed residues)")
    return sc

def oligos():
    return [
        "Original 2015 Online Methods (accepted manuscript), extracted lines 319–342:",
        f"Biotinylated oligo-dT30VN: 5′-Biotin-TEG-{CAPTURE_HANDLE}{CAPTURE_POLYT}VN-3′",
        f"Template switching oligo: 5′-{TSO_WRITTEN}-3′",
        f"PCR primer: 5′-{ISPCR}-3′",
        f"2016 author protocol TSO: 5′-{TSO_2016}-3′ (extra AT before rGrG+G; model uses 2015 TSO; reason for change unresolved)",
    ]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 lysate=g.input("lysate","single-cell lysate containing genomic DNA and poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("single-cell lysate",(role("genomic DNA"),role("poly(A) RNA"))),primary=True)
 capture=g.input("capture_oligo","biotinylated oligo(dT) capture oligo",MaterialKind.OLIGO,
                 state=(Segment("5′ biotin-TEG capture handle",CAPTURE_HANDLE),
                        Segment("dT30",CAPTURE_POLYT),Segment("VN anchor","VN",placeholder=True)),
                 source=SOURCE+" (2015 accepted-manuscript Online Methods lines 319–320)")
 hybrid=g.product("hybrid","lysate with poly(A) RNA hybridized to biotin-oligo(dT)",MaterialKind.MOLECULAR,state=capture_scene())
 g.react("hybridize","Hybridize poly(A) RNA to biotinylated oligo(dT)",StepKind.HYBRIDIZATION,((lysate,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL)),(hybrid,))
 beads=g.input("streptavidin","streptavidin magnetic beads",MaterialKind.SELECTION)
 rna=g.product("rna_fraction","bead-bound poly(A) RNA",MaterialKind.MOLECULAR,state=molecule("RNA fraction",(role("poly(A) RNA"),role("biotin-oligo(dT)"))))
 dna=g.product("dna_fraction","mRNA-depleted genomic-DNA supernatant",MaterialKind.MOLECULAR,state=molecule("DNA fraction",(role("genomic DNA"),)))
 g.react("separate","Physically separate RNA and genomic-DNA fractions",StepKind.AFFINITY_CAPTURE,((hybrid,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(rna,dna))
 rt=g.input("rt","SuperScript II reverse transcriptase (2015 Methods lines 333–334)",MaterialKind.CATALYST)
 first=g.product("first_strand","oligo(dT)-primed transcript cDNA",MaterialKind.MOLECULAR,state=rna.state)
 g.react("reverse_transcribe","Reverse-transcribe bead-captured mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
 tso=g.input("tso","2015 G&T-seq template-switching oligo",MaterialKind.OLIGO,
             state=(Segment("TSO DNA handle",TSO_DNA),Segment("3′ rGrG+G","GGG")),
             source=SOURCE+" (2015 accepted-manuscript Online Methods lines 336–337)")
 switched=g.product("switched","template-switched full-length cDNA",MaterialKind.MOLECULAR,state=template_switch_local_scene())
 g.react("template_switch","INFERRED — template-switch: untemplated C additions (length not printed) pair with the TSO rGrG+G end",StepKind.TEMPLATE_SWITCHING,((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 rf=g.input("rna_wta_forward","2015 G&T-seq PCR primer, forward landing",MaterialKind.OLIGO,
            state=(Segment("G&T universal PCR primer",ISPCR),),
            source=SOURCE+" (2015 accepted-manuscript Online Methods lines 341–342)")
 rr=g.input("rna_wta_reverse","2015 G&T-seq PCR primer, opposite landing",MaterialKind.OLIGO,
            state=(Segment("G&T universal PCR primer",ISPCR),),
            source=SOURCE+" (2015 accepted-manuscript Online Methods lines 341–342)")
 rp=g.input("rna_wta_polymerase","KAPA HiFi HotStart ReadyMix, 18 cycles (2015 Methods lines 341–345)",MaterialKind.CATALYST)
 rwta=g.product("rna_wta","whole-transcriptome-amplified cDNA",MaterialKind.MOLECULAR,state=switched.state)
 g.react("amplify_rna","Amplify the transcriptome",StepKind.PCR,((switched,InputRole.TEMPLATE),(rf,InputRole.STRUCTURAL),(rr,InputRole.STRUCTURAL),(rp,InputRole.CATALYST)),(rwta,),geometry=PCRGeometry.unavailable(rwta.state,reason="The exact primer is printed in the 2015 Methods; this generic Smart-seq2-like cDNA state does not place its two binding sites."))
 random=g.input("wga_primer","MDA or PicoPLEX WGA priming reagent",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (alternative whole-genome-amplification branches)")
 wga_pol=g.input("wga_polymerase","GenomiPhi V2 (MDA) or PicoPlex kit amplification enzyme (2015 Methods lines 348–351)",MaterialKind.CATALYST)
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
