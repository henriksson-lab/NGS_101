"""GAGE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Workflow,feature
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_extensions import proximity_product
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind
CELL1=feature("gage_cell_round1","cell_barcode","combinatorial",group="gage_cell",part="ligation round 1",whitelist="96-well adapter set")
CELL2=feature("gage_cell_round2","cell_barcode","combinatorial",group="gage_cell",part="ligation round 2",whitelist="96-well adapter set")
def payload(label,rna=False):
 p=[seg("cell barcode round 1","B"*8,"cbc",placeholder=True,feature=CELL1),seg("cell barcode round 2","C"*8,"cbc",placeholder=True,feature=CELL2)]
 if rna:p.append(seg("RT UMI","U"*8,"umi",placeholder=True,feature=feature("gage_rna_umi","umi","random")))
 p.append(seg(label,"X"*42,placeholder=True));return p
DNA,DNA_P=nextera_library(payload("proximity-ligated genomic junction"),"GAGE-seq chromatin-contact library")
RNA,RNA_P=nextera_library(payload("cDNA insert",True),"GAGE-seq RNA library")
TITLE="GAGE-seq — gene expression and genome architecture"
NOTES="01_gage-seq.html"
SOURCE='<a href="https://doi.org/10.1101/2023.07.20.549578">Zhou et al., <i>bioRxiv</i> (2023)</a>.'
SUMMARY="In fixed cells, RNA is reverse-transcribed and genomic contacts are proximity-ligated. Two rounds of well-specific adapter ligation barcode both analytes before cDNA/DNA separation and library-specific tagmentation."
CAVEAT="The source establishes the two combinatorial barcode rounds and final Nextera/TruSeq PCR design. Representative barcode bases are typed placeholders because the full adapter table is supplementary. Strand-resolved barcode overhangs and exact branch-finishing PCR-primer landing geometries have not yet been verified from that table and remain explicitly unavailable."
FINAL_LIBRARIES=(("Chromatin-contact library",DNA,DNA_P,"Proximity-ligated DNA is preamplified, tagmented and indexed.","Paired reads traverse the two genomic sides of captured contacts."),("RNA library",RNA,RNA_P,"Biotinylated cDNA is recovered, dC-tailed, amplified and tagmented.","Cell barcode rounds and the RT UMI identify each transcript observation."))
READ_LENGTHS={"Chromatin-contact library":{"Read 1":150,"Read 2":150},"RNA library":{"Read 1":150,"Read 2":150}}
def barcode_geometry():return AdapterLigationGeometry.unavailable("The GAGE-seq preprint establishes analyte-specific two-round ligation adapters, but no strand-resolved adapter table has yet been verified from the supplementary source; their overhangs are not asserted.")
def sections():return [("Fix and prepare contacts",proximity_product(left="locus A",right="locus B",bridge="proximity-ligation junction").rows(),"Crosslinked chromatin is digested, filled and ligated in situ."),("Reverse-transcribe RNA",[Row(chunks=[("biotin-poly(dT) or random hexamer — UMI — cDNA","umi",False)])],"RNA receives an RT UMI and affinity handle."),("Ligate two cell indexes",barcode_geometry().rows()+[Row(chunks=[("round-1 barcode ** round-2 barcode ** DNA or cDNA","cbc",False)])],"Separate DNA- and cDNA-compatible adapter sets preserve modality while encoding cell identity; their overhang geometry remains explicitly unavailable."),("Separate and build libraries",[Row(chunks=[("streptavidin cDNA ──► RNA library\nsupernatant DNA  ──► contact library",None,False)])],"Each branch is preamplified, Nextera-tagmented and sequenced PE150.")]
def workflow():
 s=sections();w=Workflow(MolecularState(s[0][0],tuple(s[0][1])));w.react(s[1][0],s[1][1],note=s[1][2]);w.ligate_adapter(s[2][0],barcode_geometry(),s[2][1],note=s[2][2]);return w


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","fixed cells containing chromatin contacts and RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 nuclease=g.input("restriction_nuclease","chromatin restriction nuclease",MaterialKind.CATALYST)
 digested=g.product("digested","restriction-digested crosslinked chromatin with retained RNA",MaterialKind.MOLECULAR,state=cells.state)
 g.react("digest","Digest crosslinked chromatin in situ",StepKind.DIGESTION,((cells,InputRole.TEMPLATE),(nuclease,InputRole.CATALYST)),(digested,))
 fill_polymerase=g.input("fill_polymerase","restriction-end fill-in polymerase",MaterialKind.CATALYST)
 filled=g.product("filled","filled chromatin restriction ends with retained RNA",MaterialKind.MOLECULAR,state=digested.state)
 g.react("fill","Fill restriction-fragment ends",StepKind.EXTENSION,((digested,InputRole.TEMPLATE),(fill_polymerase,InputRole.CATALYST)),(filled,))
 contact_ligase=g.input("contact_ligase","proximity-ligation DNA ligase",MaterialKind.CATALYST)
 contacts=g.product("contacts","proximity-ligated chromatin contacts with retained RNA",MaterialKind.MOLECULAR,state=sections()[0][1])
 g.react("proximity_ligate","Ligate spatially adjacent chromatin ends in situ",StepKind.PROXIMITY_LIGATION,((filled,InputRole.TEMPLATE),(contact_ligase,InputRole.CATALYST)),(contacts,))

 rt_primer=g.input("rt_primer","biotinylated UMI-bearing poly(dT)/random-hexamer RT primers",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 copied=g.product("copied","cells containing UMI-labelled biotinylated cDNA and proximity-ligated DNA",MaterialKind.MOLECULAR,state=sections()[1][1])
 g.react("reverse_transcribe","Reverse-transcribe RNA with the affinity-tagged UMI primers",StepKind.REVERSE_TRANSCRIPTION,((contacts,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))

 barcode_ligase=g.input("barcode_ligase","cell-index DNA ligase",MaterialKind.CATALYST)
 current=copied
 reviewed=["cells","restriction_nuclease","digest","fill_polymerase","fill","contact_ligase","proximity_ligate","rt_primer","rt","reverse_transcribe","barcode_ligase"]
 for number in (1,2):
  adapter=g.input(f"barcode_{number}",f"round-{number} analyte-specific barcode adapter",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (strand-resolved supplementary adapter geometry not yet verified)")
  indexed=g.product(f"indexed_{number}",f"DNA and cDNA after cell-barcode round {number}",MaterialKind.MOLECULAR,state=sections()[2][1],inferred=True,source=SOURCE+" (barcode identity known; strand-resolved adapter geometry not yet verified)")
  g.react(f"ligate_{number}",f"Ligate cell-barcode round {number}",StepKind.LIGATION,((current,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(barcode_ligase,InputRole.CATALYST)),(indexed,),geometry=barcode_geometry())
  reviewed.extend((adapter.id,f"ligate_{number}"));current=indexed

 streptavidin=g.input("streptavidin","streptavidin affinity matrix",MaterialKind.SELECTION)
 rna=g.product("rna_branch","biotin-captured cDNA branch",MaterialKind.MOLECULAR,state=payload("cDNA insert",True))
 dna=g.product("dna_branch","proximity-ligated DNA flow-through branch",MaterialKind.MOLECULAR,state=payload("proximity-ligated genomic junction"))
 g.react("separate","Separate biotinylated cDNA from genomic contact DNA",StepKind.AFFINITY_CAPTURE,((current,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(rna,dna),changes_molecule=False)

 tdt=g.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
 rna_tailed=g.product("rna_tailed","poly(C)-tailed cell-indexed cDNA",MaterialKind.MOLECULAR,state=rna.state)
 g.react("tail_rna","Add a 3-prime poly(C) tail to cDNA",StepKind.TAILING,((rna,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(rna_tailed,))
 anchor=g.input("rna_anchor","G-rich cDNA amplification anchor",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact anchor sequence not yet transcribed)")
 anchor_polymerase=g.input("anchor_polymerase","anchor-extension polymerase",MaterialKind.CATALYST)
 rna_anchored=g.product("rna_anchored","anchor-extended cDNA",MaterialKind.MOLECULAR,state=rna_tailed.state,inferred=True,source=SOURCE+" (exact anchor sequence not yet transcribed)")
 g.react("extend_rna_anchor","Anneal and extend the cDNA amplification anchor",StepKind.EXTENSION,((rna_tailed,InputRole.TEMPLATE),(anchor,InputRole.STRUCTURAL),(anchor_polymerase,InputRole.CATALYST)),(rna_anchored,))

 dna_f=g.input("dna_preamp_forward","contact-branch preamplification forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not yet transcribed)")
 dna_r=g.input("dna_preamp_reverse","contact-branch preamplification reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (exact primer not yet transcribed)")
 dna_pol=g.input("dna_preamp_polymerase","contact-branch preamplification polymerase",MaterialKind.CATALYST)
 dna_preamp=g.product("dna_preamp","preamplified chromatin-contact products",MaterialKind.MOLECULAR,state=dna.state,inferred=True,source=SOURCE+" (exact preamplification geometry not yet transcribed)")
 g.react("preamp_dna","Preamplify the chromatin-contact branch",StepKind.PCR,((dna,InputRole.TEMPLATE),(dna_f,InputRole.STRUCTURAL),(dna_r,InputRole.STRUCTURAL),(dna_pol,InputRole.CATALYST)),(dna_preamp,),geometry=PCRGeometry.unavailable(DNA,reason="The preprint describes branch preamplification, but its exact primer landing geometry has not yet been verified from the supplement."))

 transposome=g.input("finishing_transposome","Nextera finishing transposome",MaterialKind.ADAPTER)
 tn5=g.input("finishing_tn5","Tn5 transposase",MaterialKind.CATALYST)
 branch_templates=[]
 for prefix,template,final in (("dna",dna_preamp,DNA),("rna",rna_anchored,RNA)):
  tagged=g.product(prefix+"_tagged",prefix.upper()+" branch after Nextera tagmentation",MaterialKind.MOLECULAR,state=final)
  waste=g.product(prefix+"_tn5_waste","non-library "+prefix.upper()+" tagmentation products",MaterialKind.WASTE)
  g.react("tagment_"+prefix,"Tagment the "+prefix.upper()+" branch for library completion",StepKind.TAGMENTATION,((template,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,waste))
  g.discard(waste,"does not retain the selected unlike-ended Nextera insert")
  reviewed.append("tagment_"+prefix);branch_templates.append(tagged)

 reviewed.extend(("streptavidin","separate","tdt","tail_rna","rna_anchor","anchor_polymerase","extend_rna_anchor","dna_preamp_forward","dna_preamp_reverse","dna_preamp_polymerase","preamp_dna","finishing_transposome","finishing_tn5"))
 finish_pcr_branch(g,"dna_final",branch_templates[0],DNA,DNA_P,reviewed,label="GAGE-seq chromatin-contact",read_lengths=READ_LENGTHS["Chromatin-contact library"],inferred_primers=True,limitation="Exact branch-finishing PCR primers have not yet been verified from the supplementary table.")
 finish_pcr_branch(g,"rna_final",branch_templates[1],RNA,RNA_P,reviewed,label="GAGE-seq RNA",read_lengths=READ_LENGTHS["RNA library"],inferred_primers=True,limitation="Exact branch-finishing PCR primers have not yet been verified from the supplementary table.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
