"""scCAT-seq physical RNA/nucleus separation."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import analyte_split
TITLE="scCAT-seq — transcriptome and chromatin accessibility"
NOTES="01_sccat-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-018-08205-7">Liu et al., <i>Nature Communications</i> (2019)</a>.'
SUMMARY="Gentle lysis separates cytoplasmic RNA from the nucleus of one cell; Smart-seq2 and scATAC-seq then process the matched fractions."
RNA,RP=illumina_branch("scCAT-seq RNA library","Smart-seq2 cDNA",architecture="nextera")
ATAC,AP=illumina_branch("scCAT-seq ATAC library","accessible chromatin fragment",architecture="nextera")
FINAL_LIBRARIES=(("Transcriptome library",RNA,RP,"Smart-seq2 cDNA after Nextera library construction.","Nextera-compatible primers read the RNA branch."),("Chromatin-accessibility library",ATAC,AP,"The matched nucleus yields a tagmented accessibility library.","Nextera-compatible primers read the ATAC branch."))
def sections(): return [
 ("Lyse one cell without disrupting its nucleus",[Row(chunks=[("cell → cytoplasmic RNA + intact nucleus","r1",False)])],"The nucleus is physically retained while cytoplasmic RNA enters the supernatant."),
 ("Separate the matched analytes",analyte_split("cytoplasmic RNA → Smart-seq2","nucleus → scATAC-seq"),"The two branches retain common cell provenance, not a covalent barcode."),
 ("Build the RNA branch",[Row(chunks=[("oligo(dT) RT → template switching → cDNA PCR → tagmentation","tso",False)])],"Smart-seq2 produces the transcriptome library."),
 ("Build the chromatin branch",[Row(chunks=[("nuclear chromatin → Tn5 tagmentation → indexed PCR","me",False)])],"The intact matched nucleus produces the accessibility library."),]
