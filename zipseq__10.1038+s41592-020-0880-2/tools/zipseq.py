"""ZipSeq photocaged spatial barcode model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row,Segment,oligo
from tenx_rna import SEQ_PRIMERS,poly_a_feature_library,transcript_library
ZIP1="GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*30
RNA_LIB=transcript_library("ZipSeq gene-expression library")
ZIP_LIB=poly_a_feature_library("photoaddress zipcode",feature_len=7,feature_role="spatial_barcode",feature_id="zipseq_region",name="ZipSeq zipcode library")
FINAL_LIBRARIES=(("Gene-expression library",RNA_LIB,SEQ_PRIMERS,"Standard 10x transcript endpoint.","Read 1 recovers cell barcode and UMI; Read 2 enters cDNA."),("Zipcode library",ZIP_LIB,SEQ_PRIMERS,"Captured polyadenylated zipcode links an illuminated region to the same cell barcode.","Read 1 recovers cell barcode and UMI; Read 2 identifies the zipcode."))
TITLE="ZipSeq — sequential photo-uncaging of spatial zipcodes"
NOTES="01_zipseq.html"
SOURCE='<a href="https://doi.org/10.1038/s41592-020-0880-2">Hu et al., <i>Nature Methods</i> (2020)</a>.'
SUMMARY="A photocaged DNA anchor coats cells. Patterned illumination exposes an overhang, and a complementary polyadenylated zipcode is added after each illumination round. 10x RNA-seq reads the accumulated spatial tags with the transcriptome."
def oligos(): yield oligo("zipcode 1",[Segment("uncaged-overhang complement + handle + zipcode","GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC","cbc"),Segment("poly(A)30","A"*30)])
def sections(): return [
 ("Coat cells with photocaged anchors",[Row(chunks=[("cell membrane — anchor || caged strand [overhang hidden]",None,False)])],"The anchor may be antibody- or lipid-conjugated."),
 ("Illuminate a selected region",[Row(chunks=[("405-nm light → cage removed → overhang exposed", "w1",False)])],"Only cells in the projected region become competent for that round's zipcode."),
 ("Hybridize a poly(A)-tailed zipcode",[Row(chunks=[("exposed overhang || zipcode — poly(A)30", "cbc",False)])],"Repeated uncaging and hybridization give cells a history of illuminated regions."),
 ("Dissociate and run 10x 3′ RNA-seq",[Row(chunks=[("same cell barcode ├─ endogenous transcript\n                  └─ zipcode feature",None,False)])],"Zipcodes and mRNA are captured by poly(dT) in the same droplet."),]
