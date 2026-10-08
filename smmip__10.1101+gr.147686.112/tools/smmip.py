"""Single-molecule molecular inversion probe model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from targeted_ngs import smmip_library,smmip_rows

TITLE="Single-molecule molecular inversion probes (smMIPs)"
NOTES="01_smmip.html"
SOURCE='Hiatt et al. 2013, <a href="https://doi.org/10.1101/gr.147686.112">doi:10.1101/gr.147686.112</a>.'
SUMMARY="A 12-nt-tagged inversion probe copies a target gap, circularizes, survives exonuclease and is universally amplified so reads can be collapsed by original molecule."
CAVEAT="Probe arms and captured gaps are target-specific. The invariant 12-nt tag and capture topology are published; the final diagram uses generic target bases and universal Illumina roles."
FINAL_LIBRARY,SEQ_PRIMERS=smmip_library()
FINAL_CAPTION="Universal-PCR product from one captured smMIP circle. The twelve-base tag precedes amplification and identifies the original captured molecule."
SEQUENCING_INTRO="Paired-end reads traverse probe backbone/tag and the captured target; sample indexes are introduced by universal PCR."

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
