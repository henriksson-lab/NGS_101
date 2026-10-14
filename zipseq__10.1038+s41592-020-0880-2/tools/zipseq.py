"""ZipSeq photocaged spatial barcode model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row,Segment,feature,oligo
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library,transcript_library
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role
ZIP1="GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*28
RNA_LIB=transcript_library("ZipSeq gene-expression library")
ZIP_LIB=poly_a_feature_library("photoaddress zipcode",feature_len=8,feature_role="spatial_barcode",feature_id="zipseq_region",name="ZipSeq zipcode library")
FINAL_LIBRARIES=(("Gene-expression library",RNA_LIB,SEQ_PRIMERS,"Standard 10x transcript endpoint.","Read 1 recovers cell barcode and UMI; Read 2 enters cDNA."),("Zipcode library",ZIP_LIB,SEQ_PRIMERS,"Captured polyadenylated zipcode links an illuminated region to the same cell barcode.","Read 1 recovers cell barcode and UMI; Read 2 identifies the zipcode."))
TITLE="ZipSeq — sequential photo-uncaging of spatial zipcodes"
NOTES="01_zipseq.html"
SOURCE='<a href="https://doi.org/10.1038/s41592-020-0880-2">Hu et al., <i>Nature Methods</i> (2020)</a>.'
SUMMARY="A photocaged DNA anchor coats cells. Patterned illumination exposes an overhang, and a complementary polyadenylated zipcode is added after each illumination round. 10x RNA-seq reads the accumulated spatial tags with the transcriptome."
def oligos():
 yield oligo("antibody anchor",[Segment("anchor sequence","ATCGTTTTTTTTTTTTGATGTTGACGG")],mods="Amine-C6 linker")
 yield oligo("O1 blocking strand",[Segment("caged-overhang blocker","CGATCTGTGGTTGCTAC")])
 yield oligo("zipcode 1",[Segment("O1′ complement + read handle","GTAGCAACCACAGATCGCACCCGAGAATTC"),Segment("8-base region zipcode","CATGATGC","cbc",feature=feature("zipseq_region","spatial_barcode","whitelist",whitelist="paper zipcode list")),Segment("poly(A)28","A"*28)])
def sections(): return [
 ("Coat cells with photocaged anchors",[Row(chunks=[("cell membrane — anchor || caged strand [overhang hidden]",None,False)])],"The anchor may be antibody- or lipid-conjugated."),
 ("Illuminate a selected region",[Row(chunks=[("405-nm light → cage removed → overhang exposed", "w1",False)])],"Only cells in the projected region become competent for that round's zipcode."),
 ("Hybridize a poly(A)-tailed zipcode",[Row(chunks=[("exposed overhang || zipcode — poly(A)30", "cbc",False)])],"Repeated uncaging and hybridization give cells a history of illuminated regions."),
 ("Dissociate and run 10x 3′ RNA-seq",[Row(chunks=[("same cell barcode ├─ endogenous transcript\n                  └─ zipcode feature",None,False)])],"Zipcodes and mRNA are captured by poly(dT) in the same droplet."),]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","cells containing poly(A) RNA",MaterialKind.BIOLOGICAL,state=molecule("ZipSeq cells",(role("cell membrane"),role("poly(A) RNA"))),primary=True)
 anchors=g.input("anchors","membrane-conjugated photocaged DNA anchors",MaterialKind.MOLECULAR)
 coated=g.product("coated","cells coated with photocaged anchors",MaterialKind.MOLECULAR,state=molecule("anchor-coated cells",(role("poly(A) RNA"),role("caged anchor overhang"))))
 excess=g.product("excess","unbound photocaged anchors",MaterialKind.WASTE)
 g.react("coat","Coat cells with photocaged DNA anchors",StepKind.BINDING,((cells,InputRole.TEMPLATE),(anchors,InputRole.ANALYTE)),(coated,excess)); g.discard(excess,"is not attached to a retained cell")
 light=g.input("light","patterned 405-nm illumination condition",MaterialKind.CONDITION)
 uncaged=g.product("uncaged","region-selected cells with exposed anchor overhangs",MaterialKind.MOLECULAR,state=molecule("uncaged cells",(role("poly(A) RNA"),role("exposed anchor overhang"))))
 g.react("uncage","Illuminate a selected region to expose anchor overhangs",StepKind.CHEMICAL_CONVERSION,((coated,InputRole.TEMPLATE),(light,InputRole.CONDITION)),(uncaged,))
 zipcode=g.input("zipcode","complementary poly(A)-tailed zipcode oligo",MaterialKind.OLIGO,state=ZIP1)
 tagged=g.product("tagged","cells carrying the hybridized regional zipcode",MaterialKind.MOLECULAR,state=molecule("zipcode-tagged cells",(role("poly(A) RNA"),role("region zipcode"),role("zipcode poly(A)"))))
 g.react("hybridize","Hybridize the round-specific poly(A)-tailed zipcode",StepKind.HYBRIDIZATION,((uncaged,InputRole.TEMPLATE),(zipcode,InputRole.STRUCTURAL)),(tagged,))
 droplets=g.input("droplets","10x 3-prime GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","cell-barcoded transcript and zipcode products",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Dissociate and partition tagged cells with barcoded beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 rna=g.product("rna_branch","gene-expression branch",MaterialKind.MOLECULAR,state=occupied.state); zips=g.product("zip_branch","zipcode-feature branch",MaterialKind.MOLECULAR,state=occupied.state)
 g.react("split","Separate transcript and zipcode products",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(rna,zips),changes_molecule=False)
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("cdna","cell-barcoded transcript cDNA",MaterialKind.MOLECULAR,state=molecule("ZipSeq cDNA",(role("cell barcode"),role("UMI"),role("cDNA"))))
 g.react("reverse_transcribe","Capture and reverse-transcribe endogenous mRNA",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 zpol=g.input("zipcode_extension","zipcode barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied","cell-barcoded zipcode products",MaterialKind.MOLECULAR,state=molecule("ZipSeq zipcode precursor",(role("cell barcode"),role("UMI"),role("zipcode"))))
 g.react("copy_zipcodes","Extend the same bead primer across poly(A)-tailed zipcodes",StepKind.EXTENSION,((zips,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(zpol,InputRole.CATALYST)),(copied,))
 reviewed=["cells","anchors","coat","light","uncage","zipcode","hybridize","droplets","partition","split","bead_primer","rt","reverse_transcribe","zipcode_extension","copy_zipcodes"]
 finish_pcr_branch(g,"rna",cdna,RNA_LIB,SEQ_PRIMERS,reviewed,label="ZipSeq gene expression")
 finish_pcr_branch(g,"zip",copied,ZIP_LIB,SEQ_PRIMERS,reviewed,label="ZipSeq zipcode")
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
