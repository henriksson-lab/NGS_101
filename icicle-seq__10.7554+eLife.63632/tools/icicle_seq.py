"""ICICLE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from antibody_tags import feature_library
from chromatin_epigenetics import droplet_targeted_tagment_library
from chemdraw import Row, Segment, oligo
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role
ATAC,ATAC_P=droplet_targeted_tagment_library("accessible DNA",assay_name="poly(dT)/R2N capture branch",assay_barcode_len=8,assay_is_feature=False)
ADT,ADT_P=feature_library("ICICLE-seq antibody-tag library",barcode_id="icicle_antibody",cell_id="cell_barcode")
TITLE="ICICLE-seq — chromatin accessibility and epitopes"
NOTES="01_icicle-seq.html"
SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
SUMMARY="A custom Tn5 transposome with poly(dT)- or R2N-bearing ends makes ATAC fragments capturable by 10x 3′ gel beads alongside polyadenylated antibody tags."
CAVEAT="The displayed capture branch is structural: the paper's supplement defines several custom oligos, not a biological feature barcode."
FINAL_LIBRARIES=(("ATAC library",ATAC,ATAC_P,"Single-ended ATAC insert carrying a 10x cell barcode and UMI.","Read 1 is the 28-cycle cell-barcode/UMI read; Read 2 is the ATAC insert."),("Antibody-tag library",ADT,ADT_P,"ADT barcode from the same GEM is selectively amplified.","Read 1 carries cell barcode and UMI; Read 2 identifies the antibody."))
READ_LENGTHS={"ATAC library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100},"Antibody-tag library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":100}}
def oligos():
    """Representative ordered oligos from Supplementary Table 5."""
    yield oligo("poly(dT) Tn5 top", [Segment("poly(dT)30", "T"*30), Segment("degenerate VN", "VN", placeholder=True), Segment("mosaic end", "AGATGTGTATAAGAGACAG")])
    yield oligo("Tn5 bottom", [Segment("mosaic end", "CTGTCTCTTATACACATCT")], mods="/5Phos/")
    yield oligo("Tn5 S7 top", [Segment("Nextera S7 mosaic end", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")])
    yield oligo("ADT pre-index PCR reverse primer", [Segment("ADT reverse primer", "CCTTGGCACCCGAGAATTCC")])
    yield oligo("ADT i7 primer", [Segment("P7 and index flank", "CAAGCAGAAGACGGCATACGAGAT"), Segment("8-base i7 index", "NNNNNNNN", placeholder=True), Segment("ADT-binding arm", "GTGACTGGAGTTCCTTGGCACCCGAGAATTCC"), Segment("terminal phosphorothioate A", "*A", placeholder=True)])
def sections():return [("Stain cells with TotalSeq-A antibodies",[Row(chunks=[("antibody — ADT barcode — poly(A)","cbc",False)])],"Antibody-derived tags remain associated with each cell."),("Tagment with custom Tn5",[Row(chunks=[("poly(dT) or R2N capture end ** accessible DNA","me",False)])],"The custom transposome makes genomic fragments compatible with the 10x 3′ bead capture reaction."),("Barcode in GEMs",[Row(chunks=[("10x cell barcode + UMI → ATAC fragments and ADTs","cbc",False)])],"Both analytes receive one cell identity."),("Split and amplify",[Row(chunks=[("barcoded product ├─ ATAC PCR\n                 └─ ADT PCR",None,False)])],"The two libraries are sequenced separately.")]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","cells containing accessible chromatin and surface epitopes",MaterialKind.BIOLOGICAL,state=molecule("ICICLE-seq cells",(role("accessible chromatin"),role("surface epitopes"))),primary=True)
 antibodies=g.input("antibodies","poly(A)-tailed TotalSeq-A antibodies",MaterialKind.SELECTION)
 stained=g.product("stained","cells retaining antibody-derived tags",MaterialKind.MOLECULAR,state=molecule("stained cells",(role("accessible chromatin"),role("antibody tags"))))
 unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
 g.react("stain","Stain cells with TotalSeq-A antibodies",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound)); g.discard(unbound,"does not report a cell-associated epitope")
 transposome=g.input("transposome","poly(dT)- or R2N-ended custom transposome",MaterialKind.ADAPTER)
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","cells retaining capturable ATAC fragments and antibody tags",MaterialKind.MOLECULAR,state=stained.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin with custom capture ends",StepKind.TAGMENTATION,((stained,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not retain a capturable accessible-DNA fragment")
 droplets=g.input("droplets","10x 3-prime GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded ATAC and antibody-tag products",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 atac=g.product("atac_branch","ATAC branch",MaterialKind.MOLECULAR,state=occupied.state); adt=g.product("adt_branch","antibody-tag branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate ATAC and antibody-tag products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(atac,adt),changes_molecule=False)
 atacprimer=g.input("atac_capture","ATAC-to-bead capture primer",MaterialKind.OLIGO)
 atacpol=g.input("atac_extension","ATAC barcode-transfer polymerase",MaterialKind.CATALYST)
 copied_atac=g.product("copied_atac","cell-barcode/UMI-tagged ATAC fragments",MaterialKind.MOLECULAR,state=molecule("ICICLE ATAC precursor",(role("cell barcode"),role("UMI"),role("accessible DNA"))))
 g.react("copy_atac","Transfer bead barcode and UMI to custom ATAC fragments",StepKind.EXTENSION,((atac,InputRole.TEMPLATE),(atacprimer,InputRole.STRUCTURAL),(atacpol,InputRole.CATALYST)),(copied_atac,))
 tagprimer=g.input("tag_capture","antibody-tag capture primer",MaterialKind.OLIGO)
 tagpol=g.input("tag_extension","antibody-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied_tag=g.product("copied_tag","cell-barcode/UMI-tagged antibody tags",MaterialKind.MOLECULAR,state=molecule("ICICLE ADT precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
 g.react("copy_tags","Transfer bead barcode and UMI to antibody tags",StepKind.EXTENSION,((adt,InputRole.TEMPLATE),(tagprimer,InputRole.STRUCTURAL),(tagpol,InputRole.CATALYST)),(copied_tag,))
 reviewed=["cells","antibodies","stain","transposome","tn5","tagment","droplets","partition","split","atac_capture","atac_extension","copy_atac","tag_capture","tag_extension","copy_tags"]
 finish_pcr_branch(g,"atac",copied_atac,ATAC,ATAC_P,reviewed,label="ICICLE-seq ATAC",read_lengths=READ_LENGTHS["ATAC library"])
 finish_pcr_branch(g,"adt",copied_tag,ADT,ADT_P,reviewed,label="ICICLE-seq antibody tag",read_lengths=READ_LENGTHS["Antibody-tag library"])
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
