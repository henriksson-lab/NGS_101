"""Slide-tags molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg,truseq_library
from chemdraw import Row,feature
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="Slide-tags — spatial barcodes transferred to intact nuclei"
NOTES="01_slide-tags.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-023-06837-4">Russell et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Photocleaved oligos from a sequenced bead array diffuse into overlying tissue nuclei; 10x single-nucleus RNA or Multiome processing then links spatial tags to each nucleus barcode."
CAVEAT="The public page shows the published TAGS/SLAC transfer logic and the multiome branch. Commercial 10x bead regions remain molecular-role placeholders."
CELL=feature("slidetags_cell","cell_barcode","whitelist",whitelist="10x Chromium bead barcode")
SPATIAL=feature("slidetags_spatial","spatial_barcode","whitelist",whitelist="sequenced Slide-tags puck")
UMI=feature("slidetags_umi","umi","random")
RNA,RNA_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("nuclear cDNA","X"*40,placeholder=True)],"Slide-tags RNA library")
ATAC,ATAC_P=nextera_library([seg("cell barcode association","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible chromatin","X"*40,placeholder=True)],"Slide-tags ATAC library")
TAG,TAG_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("tag UMI","U"*10,"umi",placeholder=True,feature=UMI),seg("spatial barcode","S"*20,"cbc",placeholder=True,feature=SPATIAL)],"Slide-tags spatial-tag library")
FINAL_LIBRARIES=(
 ("RNA library",RNA,RNA_P,"10x 3-prime or Multiome gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("ATAC library (Multiome branch)",ATAC,ATAC_P,"Tagmented accessibility library from the same spatially tagged nuclei.","Paired genomic reads use Nextera sites and retain shared cell identity."),
 ("Spatial-tag library",TAG,TAG_P,"Feature-barcode-style amplification links released puck tags to 10x nucleus barcodes.","The cell barcode groups spatial-tag UMIs; the tag sequence maps to bead coordinates."),
)
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library (Multiome branch)":{"Read 1":50,"Read 2":50},"Spatial-tag library":{"Read 1":28,"Read 2":50}}

def sections(): return [
 ("Place tissue on a spatially sequenced bead array",[Row(chunks=[("bead — photocleavable linker — spatial barcode — 10x-capture sequence","cbc",False)])],"Each puck barcode already has a measured x/y coordinate."),
 ("Photocleave spatial oligos into the overlying tissue",[Row(chunks=[("UV cleavage → local spatial tags enter nearby nuclei","w1",False)])],"Short diffusion transfers several local bead identities into each intact nucleus."),
 ("Dissociate tagged nuclei and partition with 10x beads",[Row(chunks=[("spatial tag + one GEM cell barcode", "cbc",False)])],"Feature-barcode capture associates puck tags with a nucleus barcode."),
 ("Prepare assay and spatial-tag libraries",modality_split("RNA library","ATAC library in Multiome experiments","spatial-tag library"),"The separate tag library supplies coordinates for the same cell barcodes found in molecular assay libraries."),
]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 array=g.input("array","spatially sequenced photocleavable bead array",MaterialKind.MOLECULAR,state=molecule("Slide-tags array",(role("mapped spatial barcode"),role("photocleavable linker"),role("10x capture sequence"))),primary=True)
 uv=g.input("uv","photocleavage condition",MaterialKind.CONDITION)
 tags=g.product("released_tags","photocleaved spatial barcode oligos",MaterialKind.MOLECULAR,state=molecule("released spatial tags",(role("spatial barcode"),role("10x capture sequence"))))
 g.react("photocleave","Photocleave spatial oligos from the mapped bead array",StepKind.CHEMICAL_CONVERSION,((array,InputRole.TEMPLATE),(uv,InputRole.CONDITION)),(tags,))
 nuclei=g.input("nuclei","intact tissue nuclei containing RNA and accessible chromatin",MaterialKind.BIOLOGICAL,state=molecule("tissue nuclei",(role("poly(A) RNA"),role("accessible chromatin"))),primary=True)
 tagged=g.product("tagged_nuclei","nuclei retaining local spatial tags",MaterialKind.MOLECULAR,state=molecule("spatially tagged nuclei",(role("spatial barcode tags"),role("poly(A) RNA"),role("accessible chromatin"))))
 unbound=g.product("unbound_tags","unbound or diffusely transferred spatial tags",MaterialKind.WASTE)
 g.react("bind_tags","Transfer local spatial tags into nearby nuclei",StepKind.BINDING,((nuclei,InputRole.TEMPLATE),(tags,InputRole.ANALYTE)),(tagged,unbound)); g.discard(unbound,"is not associated with a recovered nucleus")
 transposome=g.input("transposome","Multiome loaded transposome",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (complete commercial oligo sequence not disclosed)")
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagmented=g.product("tagmented","spatially tagged nuclei with tagmented accessible chromatin",MaterialKind.MOLECULAR,state=tagged.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin for the Multiome branch",StepKind.TAGMENTATION,((tagged,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagmented,twaste)); g.discard(twaste,"does not retain an amplifiable accessibility fragment")
 droplets=g.input("droplets","10x GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded spatially tagged nucleus products",MaterialKind.MOLECULAR,state=tagmented.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition tagged nuclei with barcoded gel beads",StepKind.PARTITION,((tagmented,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable nucleus barcode")
 rna=g.product("rna_branch","RNA branch",MaterialKind.MOLECULAR,state=occupied.state); atac=g.product("atac_branch","optional Multiome ATAC branch",MaterialKind.MOLECULAR,state=occupied.state); tag=g.product("tag_branch","spatial-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate RNA, ATAC and spatial-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,atac,tag),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded nuclear cDNA",MaterialKind.MOLECULAR,state=molecule("Slide-tags cDNA",(role("cell barcode"),role("RNA UMI"),role("nuclear cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe nuclear RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 capture=g.input("tag_capture","spatial-tag capture primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (complete commercial sequence not disclosed)")
 pol=g.input("tag_extension","spatial-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied_tags","cell-barcoded spatial-tag products",MaterialKind.MOLECULAR,state=molecule("Slide-tags tag precursor",(role("cell barcode"),role("tag UMI"),role("spatial barcode"))))
 g.react("copy_tags","Copy spatial tags with the nucleus barcode and UMI",StepKind.EXTENSION,((tag,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(copied,))
 reviewed=["array","uv","photocleave","nuclei","bind_tags","transposome","tn5","tagment","droplets","partition","split","bead_primer","rt","reverse_transcribe","tag_capture","tag_extension","copy_tags"]
 finish_pcr_branch(g,"rna",cdna,RNA,RNA_P,reviewed,label="Slide-tags RNA",read_lengths=READ_LENGTHS["RNA library"])
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="Slide-tags ATAC",read_lengths=READ_LENGTHS["ATAC library (Multiome branch)"])
 finish_pcr_branch(g,"tag",copied,TAG,TAG_P,reviewed,label="Slide-tags spatial tag",read_lengths=READ_LENGTHS["Spatial-tag library"])
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
