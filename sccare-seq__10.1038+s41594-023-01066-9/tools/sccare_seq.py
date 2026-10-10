"""scCARE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg,truseq_library
from chemdraw import Row,Segment,feature,oligo
RT_HANDLE="ACGACGCTCTTCCGATCT"; RT1="CTTAGGAC"; UMI="N"*8; POLYT="T"*30+"VN"
CELL=feature("sccare_cell","cell_barcode","whitelist",whitelist="eight printed RT barcodes")
DNA,DNA_P=nextera_library([seg("cell barcode", "B"*8,"cbc",placeholder=True,feature=CELL),seg("contact junction","X"*46,placeholder=True)],"scCARE-seq DNA-contact library")
RNA,RNA_P=truseq_library([seg("cell barcode","B"*8,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*8,"umi",placeholder=True,feature=feature("sccare_umi","umi","random")),seg("cDNA insert","X"*42,placeholder=True)],"scCARE-seq RNA library")
TITLE="scCARE-seq — chromatin architecture and RNA expression"
NOTES="01_sccare-seq.html"
SOURCE='<a href="https://doi.org/10.1038/s41594-023-01066-9">Deng et al., <i>Nature Structural & Molecular Biology</i> (2023)</a>.'
SUMMARY="A biotinylated, well-barcoded oligo(dT) labels RNA while in-situ tagmentation and proximity ligation preserve chromatin contacts. DNA and RNA branches are separated and indexed for paired analysis."
CAVEAT="Eight exact RT barcodes are printed in Supplementary Table 1; the page shows one member and models the cell index as its whitelist identity."
FINAL_LIBRARIES=(("DNA-contact library",DNA,DNA_P,"Tn5-tagged, proximity-ligated genomic fragments.","Paired genomic reads report chromatin contacts; the cell index assigns the source cell."),("RNA library",RNA,RNA_P,"Well-barcoded UMI cDNA selectively amplified with its RNA primer.","The barcode and UMI precede the cDNA insert."))
def oligos():yield oligo("Oligo-dT-1",[Segment("RNA handle",RT_HANDLE,"r1"),Segment("UMI",UMI,"umi",placeholder=True,feature=feature("sccare_umi","umi","random")),Segment("cell barcode 1",RT1,"cbc",feature=CELL),Segment("oligo(dT)VN",POLYT,placeholder=True)],mods="one dT is biotinylated in the capture version")
def sections():return [("Fix and tagment nuclei",[Row(chunks=[("Tn5-N5 ** chromatin ** Tn5-N7","me",False)])],"Printed Nextera transposome oligos install DNA handles."),("Reverse-transcribe and label RNA",[Row(chunks=[("RNA handle — UMI — well barcode — oligo(dT) → poly(A) RNA","umi",False)])],"A biotinylated dT version permits physical recovery of the RNA branch."),("Generate contact junctions",[Row(chunks=[("tagged locus A ** proximity ligation ** tagged locus B",None,False)])],"Crosslinked fragments are ligated before single-cell separation."),("Separate and amplify",[Row(chunks=[("DNA-contact PCR ← split → biotin-cDNA PCR",None,False)])],"Nextera primers recover DNA and a distinct RNA primer recovers cDNA.")]
