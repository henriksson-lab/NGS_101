"""Paired-seq joint RNA and chromatin accessibility."""
from chemdraw import Row
from multibranch import illumina_branch
TITLE="Paired-seq — combinatorial RNA and chromatin accessibility"
NOTES="01_paired-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41594-019-0323-x">Zhu et al., <i>Nature Structural &amp; Molecular Biology</i> (2019)</a>.'
SUMMARY="Permeabilized nuclei undergo RNA reverse transcription and chromatin tagmentation before four rounds of split-pool cellular indexing."
CAVEAT="Barcode rounds and branch chemistry are source-established. The compact final Illumina shells are marked inferred because the full protocol-specific molecules have not yet been transcribed base by base."
BC=(("round 1 barcode",3),("round 2 barcode",7),("round 3 barcode",7),("round 4 barcode",7))
RNA,RP=illumina_branch("Paired-seq RNA library","cDNA fragment",cell_parts=BC,umi=10,inferred_adapters=True)
ATAC,AP=illumina_branch("Paired-seq ATAC library","accessible chromatin fragment",architecture="nextera",cell_parts=BC,inferred_adapters=True)
FINAL_LIBRARIES=(("RNA library",RNA,RP,"cDNA reads carry the four-part cell barcode and UMI.","Illumina primers read the RNA branch."),("Chromatin-accessibility library",ATAC,AP,"Tagmented DNA reads carry the same four-part cell barcode.","Nextera-compatible primers read the ATAC branch."))
def sections(): return [
 ("Reverse-transcribe RNA and tagment chromatin in nuclei",[Row(chunks=[("RNA → cDNA     accessible DNA ** Tn5 ends","me",False)])],"Both analytes remain inside the same permeabilized nucleus."),
 ("Apply four split-pool barcode rounds",[Row(chunks=[("BC1 ** BC2 ** BC3 ** BC4 ** RNA/ATAC molecule","cbc",False)])],"Successive ligations create the shared combinatorial cell identity."),
 ("Add a terminal homopolymer and anchor",[Row(chunks=[("barcoded molecule—poly(C) ← G-anchor extension","tso",False)])],"Terminal transferase and anchored amplification make indexed material amplifiable."),
 ("INFERRED — complete RNA and ATAC sequencing shells",[Row(chunks=[("restriction split → RNA library | ATAC library","r2",False)])],"Branch-selective restriction is source-established; the compact outer adapter completion is inferred."),]
