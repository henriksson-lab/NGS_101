"""GAGE-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature
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
CAVEAT="The source establishes the two combinatorial barcode rounds and final Nextera/TruSeq PCR design. Representative barcode bases are typed placeholders because the full adapter table is supplementary."
FINAL_LIBRARIES=(("Chromatin-contact library",DNA,DNA_P,"Proximity-ligated DNA is preamplified, tagmented and indexed.","Paired reads traverse the two genomic sides of captured contacts."),("RNA library",RNA,RNA_P,"Biotinylated cDNA is recovered, dC-tailed, amplified and tagmented.","Cell barcode rounds and the RT UMI identify each transcript observation."))
READ_LENGTHS={"Chromatin-contact library":{"Read 1":150,"Read 2":150},"RNA library":{"Read 1":150,"Read 2":150}}
def sections():return [("Fix and prepare contacts",[Row(chunks=[("locus A ** proximity-ligation junction ** locus B",None,False)])],"Crosslinked chromatin is digested, filled and ligated in situ."),("Reverse-transcribe RNA",[Row(chunks=[("biotin-poly(dT) or random hexamer — UMI — cDNA","umi",False)])],"RNA receives an RT UMI and affinity handle."),("Ligate two cell indexes",[Row(chunks=[("round-1 barcode ** round-2 barcode ** DNA or cDNA","cbc",False)])],"Separate DNA- and cDNA-compatible adapter sets preserve modality while encoding cell identity."),("Separate and build libraries",[Row(chunks=[("streptavidin cDNA ──► RNA library\nsupernatant DNA  ──► contact library",None,False)])],"Each branch is preamplified, Nextera-tagmented and sequenced PE150.")]
