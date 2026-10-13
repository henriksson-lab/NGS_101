"""TEA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from antibody_tags import feature_library
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature
from tenx_rna import SEQ_PRIMERS,transcript_library
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role
CELL=feature("tea_cell","cell_barcode","whitelist",whitelist="10x Multiome bead barcode")
RNA=transcript_library("TEA-seq transcript library")
ATAC,ATAC_P=nextera_library([seg("shared cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("accessible DNA","X"*42,placeholder=True)],"TEA-seq ATAC library")
ADT,ADT_P=feature_library("TEA-seq antibody-tag library",barcode_id="tea_antibody",cell_id="tea_cell")
TITLE="TEA-seq — transcripts, epitopes and accessibility"
NOTES="01_tea-seq.html"
SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
SUMMARY="Oligo-tagged antibodies are combined with a 10x Single Cell Multiome workflow so RNA, surface epitopes and accessible chromatin share one cell identity but become three sequencing libraries."
CAVEAT="Commercial bead and TotalSeq-A sequences are represented structurally. The page does not claim proprietary bases."
FINAL_LIBRARIES=(("RNA library",RNA,SEQ_PRIMERS,"Poly(A) RNA captured by the Multiome bead.","Read 1 carries cell barcode and UMI; Read 2 reports cDNA."),("ATAC library",ATAC,ATAC_P,"Tagmented accessible chromatin carrying the same cell identity.","Paired genomic reads use Nextera sites."),("Antibody-tag library",ADT,ADT_P,"Polyadenylated antibody tags selectively amplified from the shared product.","Barcode and UMI identify the epitope observation."))
READ_LENGTHS={"RNA library":{"Read 1":28,"Read 2":90},"ATAC library":{"Read 1":51,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":51},"Antibody-tag library":{"Read 1":28,"Read 2":50}}
def sections():return [("Stain permeabilized cells",[Row(chunks=[("epitope — antibody — feature barcode — poly(A)","cbc",False)])],"TotalSeq-A oligos encode surface epitopes."),("Tagment and partition",[Row(chunks=[("accessible DNA + poly(A) RNA + antibody tags → one barcoded GEM","me",False)])],"The Multiome workflow associates all three analytes with the same cell."),("Split by analyte",[Row(chunks=[("shared product ├─ RNA\n               ├─ ATAC\n               └─ antibody tags",None,False)])],"Modality-specific PCRs create three libraries.")]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","permeabilized cells containing RNA, chromatin and surface epitopes",MaterialKind.BIOLOGICAL,state=molecule("TEA-seq cells",(role("poly(A) RNA"),role("accessible chromatin"),role("surface epitopes"))),primary=True)
 antibodies=g.input("antibodies","poly(A)-tailed TotalSeq-A antibodies",MaterialKind.SELECTION,inferred=True,source=SOURCE+" (commercial oligo bases proprietary)")
 stained=g.product("stained","cells retaining epitope-specific antibody tags",MaterialKind.MOLECULAR,state=molecule("stained cells",(role("poly(A) RNA"),role("accessible chromatin"),role("antibody tags"))))
 unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
 g.react("stain","Stain permeabilized cells with oligo-tagged antibodies",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound)); g.discard(unbound,"does not report a cell-associated epitope")
 transposome=g.input("transposome","Multiome loaded transposome",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (commercial oligo bases proprietary)")
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","stained cells retaining tagmented chromatin and RNA",MaterialKind.MOLECULAR,state=stained.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin",StepKind.TAGMENTATION,((stained,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not retain an amplifiable accessibility fragment")
 droplets=g.input("droplets","10x Multiome GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded RNA, ATAC and antibody-tag products",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 rna=g.product("rna_branch","RNA branch",MaterialKind.MOLECULAR,state=occupied.state); atac=g.product("atac_branch","ATAC branch",MaterialKind.MOLECULAR,state=occupied.state); adt=g.product("adt_branch","antibody-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate RNA, ATAC and antibody-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,atac,adt),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (commercial oligo bases proprietary)")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=molecule("TEA-seq cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe poly(A) RNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 tagprimer=g.input("tag_primer","antibody-tag capture primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (commercial oligo bases proprietary)")
 tagpol=g.input("tag_extension","antibody-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied_tags","cell-barcoded antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("TEA-seq ADT precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
 g.react("copy_tags","Transfer cell identity to antibody tags",StepKind.EXTENSION,((adt,InputRole.TEMPLATE),(tagprimer,InputRole.STRUCTURAL),(tagpol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","antibodies","stain","transposome","tn5","tagment","droplets","partition","split","bead_primer","rt","reverse_transcribe","tag_primer","tag_extension","copy_tags"]
 finish_pcr_branch(g,"rna",cdna,RNA,SEQ_PRIMERS,reviewed,label="TEA-seq RNA",read_lengths=READ_LENGTHS["RNA library"],inferred_primers=True)
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="TEA-seq ATAC",read_lengths=READ_LENGTHS["ATAC library"],inferred_primers=True)
 finish_pcr_branch(g,"adt",copied,ADT,ADT_P,reviewed,label="TEA-seq antibody tag",read_lengths=READ_LENGTHS["Antibody-tag library"],inferred_primers=True)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
