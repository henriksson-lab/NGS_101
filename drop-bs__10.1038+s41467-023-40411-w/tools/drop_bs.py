"""Drop-BS molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from base_conversion import bisulfite_path
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature,oligo

CELL_LEN=15
RANDOM_PRIMER="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"+"N"*9
FINAL_LIBRARY,SEQ_PRIMERS=nextera_library([
 seg("bisulfite-converted genomic insert","X"*38,placeholder=True),
 seg("bead-ligation handle","X"*20,placeholder=True),
 seg("15-nt bead cell barcode","D"*CELL_LEN,"cbc",placeholder=True,
     feature=feature("dropbs_cell","cell_barcode","whitelist",whitelist="Drop-BS bead barcode set")),
],"Drop-BS library")
TITLE="Drop-BS — droplet barcoding before bisulfite conversion"
NOTES="01_drop-bs.html"
SOURCE='<a href="https://doi.org/10.1038/s41467-023-40411-w">Fu et al., <i>Nature Communications</i> (2023)</a>.'
SUMMARY="MNase-fragmented DNA and a photocleavable barcode bead meet in one droplet. Ligation attaches cell identity before a second droplet conversion, and random priming plus indexed PCR completes the library."
FINAL_CAPTION="Nextera-ended Drop-BS library carrying the bead-derived 15-nt cell barcode on the converted genomic molecule."
SEQUENCING_INTRO="Paired-end genomic reads retain the bead barcode; i5 and i7 are added during the final indexing PCR."
def oligos(): yield oligo("random primer",[seg("S5+ME","TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG","me"),seg("N9","N"*9,placeholder=True)])
def sections(): return [
 ("Lyse and fragment one nucleus in a droplet",[Row(chunks=[("nucleus + MNase → genomic fragments",None,False)])],"MNase fragmentation is optimized independently of conversion."),
 ("Fuse with one barcode-bead droplet",[Row(chunks=[("UV release: bead barcode ** genomic fragment", "cbc",False)])],"End repair and ligation append the photocleaved 15-nt bead barcode."),
 ("Bisulfite-convert barcoded DNA in droplets",[Row(chunks=[(bisulfite_path(protected=False).text(),None,False)]),Row(chunks=[(bisulfite_path(protected=True).text(),"w1",False)])],"Cell identity is already covalent before droplets are broken."),
 ("Random-prime and index",[Row(chunks=[("S5—ME—N9 → extension; i5/i7 PCR completes P5/P7",None,False)])],"The printed N9 primer installs the second amplification handle."),]
