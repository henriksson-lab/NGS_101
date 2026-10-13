"""SEC-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library,transcript_library
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="SEC-seq — secretion-linked single-cell RNA sequencing"
NOTES="01_sec-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-39367-8">Cheng et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="A cell is held in an antibody-coated nanovial that captures its secreted protein. An oligo-barcoded detection antibody and the cell transcriptome are then read together with 10x Chromium."
CAVEAT="The publication uses commercial TotalSeq oligos; their feature identity and capture architecture are public, but the product sequence is proprietary."
RNA=transcript_library("SEC-seq gene-expression library")
SEC=poly_a_feature_library("secretion-antibody barcode",feature_len=15,feature_id="sec_antibody",name="SEC-seq secretion-tag library")
FINAL_LIBRARIES=(("Gene-expression library",RNA,SEQ_PRIMERS,"10x 3′ v3.1 transcript library from the nanovial-bound cell.","Read 1 records cell barcode and UMI; Read 2 reports cDNA."),("Secretion-tag library",SEC,SEQ_PRIMERS,"Feature-barcode library from the oligo-labeled antibody bound to captured secretion.","Read 1 records the same cell barcode and UMI; Read 2 identifies the secretion antibody."))
READ_LENGTHS={"Gene-expression library":{"Read 1":28,"Index 1 (i7)":10,"Read 2":90},"Secretion-tag library":{"Read 1":28,"Index 1 (i7)":10,"Read 2":90}}
def sections(): return [
 ("Capture one cell in a nanovial",[Row(chunks=[("antibody-coated cavity ← cell",None,False)])],"The cavity retains the cell and concentrates secreted protein on its surface."),
 ("Label captured secretion",[Row(chunks=[("nanovial — secreted IgG — detection antibody — DNA feature barcode","cbc",False)])],"An oligo-barcoded secondary antibody converts secretion into a sequenceable feature."),
 ("Partition cell-loaded nanovials",[Row(chunks=[("one GEM barcode → cellular mRNA + nanovial antibody tags","cbc",False)])],"The nanovial fits through the Chromium microfluidic channel and does not replace the standard bead chemistry."),
 ("Amplify two libraries",[Row(chunks=[("barcoded cDNA ├─ gene expression\n              └─ secretion feature tags",None,False)])],"The paper pools 80% RNA and 20% antibody-tag library for sequencing."),]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","cells containing poly(A) RNA and secreted protein analytes",MaterialKind.BIOLOGICAL,state=molecule("SEC-seq cells",(role("poly(A) RNA"),role("secreted proteins"))),primary=True)
 nanovials=g.input("nanovials","capture-antibody-coated nanovials",MaterialKind.SELECTION)
 loaded=g.product("loaded","cell-loaded nanovials retaining secreted protein",MaterialKind.MOLECULAR,state=molecule("loaded nanovial",(role("cell"),role("captured secretion"),role("poly(A) RNA"))))
 uncaptured=g.product("uncaptured","uncaptured secretion and empty nanovials",MaterialKind.WASTE)
 g.react("capture","Capture cells and their secretion in antibody-coated nanovials",StepKind.BINDING,((cells,InputRole.TEMPLATE),(nanovials,InputRole.SELECTION)),(loaded,uncaptured)); g.discard(uncaptured,"does not associate secretion with a retained cell-loaded nanovial")
 detection=g.input("detection_antibody","oligo-barcoded secretion-detection antibody",MaterialKind.SELECTION,inferred=True,source=SOURCE+" (commercial TotalSeq sequence proprietary)")
 labelled=g.product("labelled","nanovials bearing sequenceable secretion tags",MaterialKind.MOLECULAR,state=molecule("labelled nanovial",(role("cell"),role("poly(A) RNA"),role("secretion antibody tag"))))
 unbound=g.product("unbound","unbound detection-antibody tags",MaterialKind.WASTE)
 g.react("label","Label captured secretion with oligo-barcoded antibodies",StepKind.BINDING,((loaded,InputRole.TEMPLATE),(detection,InputRole.SELECTION)),(labelled,unbound)); g.discard(unbound,"does not report nanovial-captured secretion")
 droplets=g.input("droplets","10x Chromium GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded transcript and secretion-tag products",MaterialKind.MOLECULAR,state=labelled.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cell-loaded nanovials with barcoded beads",StepKind.PARTITION,((labelled,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable nanovial-associated cell barcode")
 rna=g.product("rna_branch","gene-expression branch",MaterialKind.MOLECULAR,state=occupied.state); sec=g.product("sec_branch","secretion-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate transcript and secretion-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,sec),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (commercial sequence proprietary)")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=molecule("SEC-seq cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe cellular mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 tagprimer=g.input("tag_primer","secretion-tag capture primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (commercial sequence proprietary)")
 tagpol=g.input("tag_extension","secretion-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied","cell-barcoded secretion-antibody tags",MaterialKind.MOLECULAR,state=molecule("SEC-seq tag precursor",(role("cell barcode"),role("UMI"),role("secretion antibody barcode"))))
 g.react("copy_tags","Transfer cell identity to secretion-antibody tags",StepKind.EXTENSION,((sec,InputRole.TEMPLATE),(tagprimer,InputRole.STRUCTURAL),(tagpol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","nanovials","capture","detection_antibody","label","droplets","partition","split","bead_primer","rt","reverse_transcribe","tag_primer","tag_extension","copy_tags"]
 finish_pcr_branch(g,"rna",cdna,RNA,SEQ_PRIMERS,reviewed,label="SEC-seq gene expression",read_lengths=READ_LENGTHS["Gene-expression library"],inferred_primers=True)
 finish_pcr_branch(g,"sec",copied,SEC,SEQ_PRIMERS,reviewed,label="SEC-seq secretion tag",read_lengths=READ_LENGTHS["Secretion-tag library"],inferred_primers=True)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
