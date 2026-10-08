"""EMP ITS1f–ITS2 fusion-primer library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
from batch_ngs import seg
from chemdraw import Construct,Row,revcomp

TITLE="Earth Microbiome Project ITS amplicon sequencing"
NOTES="01_emp-its.html"
SOURCE='<a href="https://earthmicrobiome.org/protocols-and-standards/its/">Earth Microbiome Project ITS Illumina Amplicon Protocol (EMP.ITSkabir)</a>.'
SUMMARY="Fusion primers amplify fungal ITS1 while directly adding an Illumina arm and a ten-base Golay sample barcode; extended custom read primers interrogate the low-complexity amplicon."
CAVEAT="The maintained protocol page gives the ordered fusion primers and the custom-primer extensions, but not every complete sequencing-primer oligo. Binding roles are therefore explicit while unpublished bases remain unknown."
SEQ_PRIMERS=()
FWD_ARM="AATGATACGGCGACCACCGAGATCTACAC"
FWD_LINK="GG"
FWD_LOCUS="CTTGGTCATTTAGAGGAAGTAA"
REV_ARM="CAAGCAGAAGACGGCATACGAGAT"
REV_LINK="CG"
REV_LOCUS="GCTGCGTTCTTCATCGATGC"
FINAL_LIBRARY=Construct([
 seg("P5-side fusion arm",FWD_ARM,"p5"),seg("forward linker",FWD_LINK),
 seg("ITS1f site",FWD_LOCUS,"r1"),seg("fungal ITS1 insert","N"*34,placeholder=True),
 seg("ITS2 site reverse complement",revcomp(REV_LOCUS),"r2"),seg("reverse linker",revcomp(REV_LINK)),
 seg("10-nt Golay barcode","B"*10,"cbc",placeholder=True),
 seg("P7-side fusion arm reverse complement",revcomp(REV_ARM),"p7")],name="EMP ITS fusion-primer library")
FINAL_CAPTION="EMP.ITSkabir fusion-primer product; orientation follows the published ordered primers and exposes the inline Golay barcode."
SEQUENCING_ENDING="Custom forward and reverse primers anneal across the respective Illumina-arm/locus boundaries and extend 19 nt and 15 nt farther at their 3′ ends than the amplification-primer boundaries. The Golay barcode is read with the corresponding custom index primer. Complete oligo bases are not printed on the maintained protocol page, so the page does not fabricate exact placements."

def sections():
 return [
  ("One PCR installs the complete fusion-primer construct",[
   Row(chunks=[("P5-side arm—GG—ITS1f → fungal ITS1 ← ITS2—CG—Golay10—P7-side arm", "cbc", False)]),
  ],"Both locus selection and platform-arm installation occur in the same amplification."),
  ("Custom sequencing primers cross into the locus",[
   Row(chunks=[("forward custom primer: Illumina landing region + 19 additional 3' bases →", "r1", False)]),
   Row(chunks=[("← reverse custom primer: Illumina landing region + 15 additional 3' bases", "r2", False)]),
  ],"The extra locus-matching bases raise primer melting temperature; standard sequencing primers are not interchangeable here."),
 ]
