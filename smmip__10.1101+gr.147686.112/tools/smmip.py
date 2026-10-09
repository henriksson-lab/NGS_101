"""Single-molecule molecular inversion probe model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
from targeted_ngs import smmip_rows

TITLE="Single-molecule molecular inversion probes (smMIPs)"
NOTES="01_smmip.html"
SOURCE='Hiatt et al. 2013, <a href="https://doi.org/10.1101/gr.147686.112">doi:10.1101/gr.147686.112</a>; detailed MIP protocol: <a href="https://doi.org/10.1007/978-1-4939-6442-0_6">O\'Roak et al. 2017</a>.'
SUMMARY="A 12-nt-tagged inversion probe copies a target gap, circularizes, survives exonuclease and is universally amplified so reads can be collapsed by original molecule."
CAVEAT="Probe arms and captured gaps are target-specific. The model uses the published invariant backbone/run primers and leaves only the arms, target, twelve-base molecular tag and eight-base sample index variable."
MOLECULAR_TAG=feature("molecular_tag","umi","random")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")

FORWARD_HANDLE="ATACGAGATCCGTAATCGGGAAGCTGAAG"
REVERSE_HANDLE_TOP="ACACTACCGTCGGATCGTGCGTGT"
READ1_PRIMER="CATACGAGATCCGTAATCGGGAAGCTGAAG"
READ2_PRIMER="ACACGCACGATCCGACGGTAGTGT"
INDEX1_PRIMER=REVERSE_HANDLE_TOP

FINAL_LIBRARY=Construct([
 seg("P5",il.P5,"p5"),
 seg("smMIP forward handle",FORWARD_HANDLE,"r1"),
 seg("extension targeting arm","X"*18,placeholder=True),
 seg("captured target","N"*34,placeholder=True),
 seg("ligation targeting arm","Y"*22,placeholder=True),
 seg("12-nt single-molecule tag","U"*12,"umi",placeholder=True,feature=MOLECULAR_TAG),
 seg("smMIP reverse handle",REVERSE_HANDLE_TOP,"r2"),
 seg("eight-base i7 index reverse complement","I"*8,"cbc",placeholder=True,feature=I7_FEATURE),
 seg("P7 reverse complement",il.P7_RC,"p7"),
],name="single-index smMIP capture product")
SEQ_PRIMERS=(
 sp.custom("Read 1","MIP forward sequencing primer",READ1_PRIMER,
           "O'Roak et al. 2017 detailed protocol",
           "The 5'-terminal C is an unpaired flap; the 31-nt 3' end anneals exactly."),
 sp.custom("Index 1 (i7)","MIP sample-index sequencing primer",INDEX1_PRIMER,
           "O'Roak et al. 2017 detailed protocol"),
 sp.custom("Read 2","MIP reverse sequencing primer",READ2_PRIMER,
           "O'Roak et al. 2017 detailed protocol"),
)
if problems:=sp.verify(FINAL_LIBRARY,SEQ_PRIMERS,
                        required_roles=("Read 1","Index 1 (i7)","Read 2")):
 raise ValueError("smMIP final library: "+"; ".join(problems))

FINAL_CAPTION="Single-index smMIP amplicon. Read 2 encounters the twelve-base molecular tag first; the separate eight-base index read identifies the sample."
SEQUENCING_INTRO="The original assay collects paired reads plus one eight-base sample-index read using three custom primers. It has no i5/Index 2 read."

def sections():
 return [
  ("Hybridize, fill and close one probe",smmip_rows(),
   "Both targeting arms must bind one genomic molecule. Extension across the gap and ligation convert only successful captures into closed circles."),
  ("Digest everything linear",[
   Row(chunks=[("unreacted probe + genomic DNA + incompletely closed products",None,False)]),
   Row(chunks=[("                         ↓ exonuclease",None,False)]),
   Row(chunks=[("closed smMIP circles survive as the selected template pool", "umi", False)]),
  ],"Circular topology, rather than bead affinity, performs the enrichment."),
  ("Amplify from the invariant backbone",[
   Row(chunks=[("circle—[12-nt tag]—captured gap → universal PCR → indexed library", "umi", False)]),
   Row(chunks=[("same tag across PCR descendants → single-molecule consensus", "cbc", False)]),
  ],"The molecular tag is already covalently attached before PCR, allowing family consensus and molecule counting."),
 ]
