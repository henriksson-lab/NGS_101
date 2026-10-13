"""scCARE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg,truseq_library
from chemdraw import Row,Segment,feature,oligo
from multimodal_spatial import finish_pcr_branch
from protocol_extensions import proximity_product
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,StepKind
RT_HANDLE="ACGACGCTCTTCCGATCT"; RT1="CTTAGGAC"; UMI="N"*8; POLYT="T"*30+"VN"
CELL=feature("sccare_cell","cell_barcode","whitelist",whitelist="eight printed RT barcodes")
DNA,DNA_P=nextera_library([seg("cell barcode", "B"*8,"cbc",placeholder=True,feature=CELL),seg("contact junction","X"*46,placeholder=True)],"scCARE-seq DNA-contact library")
RNA,RNA_P=truseq_library([seg("cell barcode","B"*8,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*8,"umi",placeholder=True,feature=feature("sccare_umi","umi","random")),seg("cDNA insert","X"*42,placeholder=True)],"scCARE-seq RNA library")
TITLE="scCARE-seq — chromatin architecture and RNA expression"
NOTES="01_sccare-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41594-023-01066-9">Deng et al., <i>Nature Structural & Molecular Biology</i> (2023)</a>.'
SUMMARY="A biotinylated, well-barcoded oligo(dT) labels RNA while in-situ tagmentation and proximity ligation preserve chromatin contacts. DNA and RNA branches are separated and indexed for paired analysis."
CAVEAT="Eight exact RT barcodes are printed in Supplementary Table 1; the page shows one member and models the cell index as its whitelist identity. The exact terminal branch-PCR primer pairs have not yet been transcribed and remain explicitly unavailable."
FINAL_LIBRARIES=(("DNA-contact library",DNA,DNA_P,"Tn5-tagged, proximity-ligated genomic fragments.","Paired genomic reads report chromatin contacts; the cell index assigns the source cell."),("RNA library",RNA,RNA_P,"Well-barcoded UMI cDNA selectively amplified with its RNA primer.","The barcode and UMI precede the cDNA insert."))
def oligos():yield oligo("Oligo-dT-1",[Segment("RNA handle",RT_HANDLE,"r1"),Segment("UMI",UMI,"umi",placeholder=True,feature=feature("sccare_umi","umi","random")),Segment("cell barcode 1",RT1,"cbc",feature=CELL),Segment("oligo(dT)VN",POLYT,placeholder=True)],mods="one dT is biotinylated in the capture version")
def sections():return [("Fix and tagment nuclei",[Row(chunks=[("Tn5-N5 ** chromatin ** Tn5-N7","me",False)])],"Printed Nextera transposome oligos install DNA handles."),("Reverse-transcribe and label RNA",[Row(chunks=[("RNA handle — UMI — well barcode — oligo(dT) → poly(A) RNA","umi",False)])],"A biotinylated dT version permits physical recovery of the RNA branch."),("Generate contact junctions",proximity_product(left="tagged locus A",right="tagged locus B",bridge="proximity-ligation junction").rows(),"Crosslinked fragments are ligated before single-cell separation."),("Separate and amplify",[Row(chunks=[("DNA-contact PCR ← split → biotin-cDNA PCR",None,False)])],"Nextera primers recover DNA and a distinct RNA primer recovers cDNA.")]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nuclei=g.input("nuclei","fixed nuclei containing chromatin and poly(A) RNA",MaterialKind.BIOLOGICAL,state=sections()[0][1],primary=True)
 transposome=g.input("transposome","printed N5/N7-loaded Nextera transposome",MaterialKind.ADAPTER)
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","nuclei retaining Tn5-tagged chromatin and RNA",MaterialKind.MOLECULAR,state=sections()[0][1])
 twaste=g.product("tn5_waste","non-library chromatin tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment chromatin in fixed nuclei",StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste));g.discard(twaste,"does not retain a selected N5/N7 chromatin fragment")
 rt_primer=g.input("rt_primer","biotinylated well-barcoded UMI oligo(dT) primer",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 copied=g.product("copied","nuclei containing barcoded biotin-cDNA and tagged chromatin",MaterialKind.MOLECULAR,state=sections()[1][1])
 g.react("reverse_transcribe","Reverse-transcribe and label RNA",StepKind.REVERSE_TRANSCRIPTION,((tagged,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 ligase=g.input("contact_ligase","in-situ proximity ligase",MaterialKind.CATALYST)
 contacts=g.product("contacts","nuclei containing proximity-ligated chromatin and biotin-cDNA",MaterialKind.MOLECULAR,state=sections()[2][1])
 g.react("proximity_ligate","Generate chromatin-contact junctions",StepKind.PROXIMITY_LIGATION,((copied,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(contacts,))
 streptavidin=g.input("streptavidin","streptavidin affinity matrix",MaterialKind.SELECTION)
 rna=g.product("rna_branch","biotin-captured cDNA branch",MaterialKind.MOLECULAR,state=RNA)
 dna=g.product("dna_branch","DNA-contact flow-through branch",MaterialKind.MOLECULAR,state=DNA)
 g.react("separate","Separate biotin-cDNA from DNA-contact products",StepKind.AFFINITY_CAPTURE,((contacts,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),(rna,dna),changes_molecule=False)
 reviewed=["nuclei","transposome","tn5","tagment","rt_primer","rt","reverse_transcribe","contact_ligase","proximity_ligate","streptavidin","separate"]
 finish_pcr_branch(g,"dna",dna,DNA,DNA_P,reviewed,label="scCARE-seq DNA-contact",inferred_primers=True,limitation="The exact DNA-contact branch PCR-primer pair has not yet been transcribed from the supplementary table.")
 finish_pcr_branch(g,"rna",rna,RNA,RNA_P,reviewed,label="scCARE-seq RNA",inferred_primers=True,limitation="The exact RNA branch PCR-primer pair has not yet been transcribed from the supplementary table.")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
