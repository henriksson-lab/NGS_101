"""scMT-seq joint RNA and methylome workflow."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scMT-seq — transcriptome and DNA methylome"
NOTES="01_scmt-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1186/s13059-016-0950-z">Hu et al., <i>Genome Biology</i> (2016)</a>.'
SUMMARY="Cytoplasmic RNA and nuclear genomic DNA from one cell are physically separated for transcriptome sequencing and single-cell methylome sequencing."
RNA,RP=illumina_branch("scMT-seq RNA library","single-cell cDNA",architecture="nextera")
DNA,DP=illumina_branch("scMT-seq methylome library","bisulfite-converted genomic DNA")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Amplified single-cell cDNA is converted to an Illumina library.","Declared primers read the RNA branch."),("DNA methylome library",DNA,DP,"Bisulfite-converted nuclear DNA yields the methylome library.","Declared primers read the DNA branch."))
def sections(): return [
 ("Lyse one cell and retain its nucleus",[Row(chunks=[("cell → cytoplasmic RNA + nuclear genomic DNA","r1",False)])],"The fractions are separated before either amplification."),
 ("Separate the matched analytes",analyte_split("cytoplasmic RNA → transcriptome","nuclear DNA → methylome"),"Cell identity is maintained by plate position."),
 ("Convert the RNA fraction",[Row(chunks=[("reverse transcription → cDNA amplification → Illumina library","tso",False)])],"The cytoplasmic branch reports expression."),
 ("Convert the DNA fraction",[Row(chunks=[("bisulfite conversion → random priming → indexed library","r2",False)])],"The nuclear branch reports methylated cytosines."),]
