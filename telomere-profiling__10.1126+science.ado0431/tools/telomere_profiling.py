"""Telomere Profiling enrichment workflow (Karimian et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Row,Segment,feature,feature_rows,strand_row
from telomere import AffinityCapture,PhasedTelomereAdapter,TelomereEnd

TELOTAG_BARCODE=Segment(
    "24-nt TeloTag sample barcode", "N"*24, "cbc", placeholder=True,
    feature=feature("telotag_sample_index", "sample_index", "whitelist",
                    whitelist="Karimian et al. Extended Data Table 2"))
TELOTAG_CORE=Segment("TeloTag release-site flank", "N"*4, "r3", placeholder=True)
TELOTAG=PhasedTelomereAdapter(
    TELOTAG_BARCODE.top+TELOTAG_CORE.top,
    core_segments=(TELOTAG_BARCODE,TELOTAG_CORE))
TELOTAG_CONSTRUCT=Construct([TELOTAG_BARCODE,TELOTAG_CORE],name="barcoded TeloTag core")
CAPTURE=AffinityCapture("biotin","adapter-encoded restriction site")

def panels(): return [("Digest high-molecular-weight genomic DNA away from telomere ends",[Row(chunks=[("internal restriction cut | subtelomere -- telomere -- native 3′ overhang",None,False)])],"Restriction reduces non-telomeric length while leaving the native chromosome terminus available for tagging."),("Anneal six telomeric splints and ligate the biotinylated TeloTag",[*TELOTAG.annealed_scene(TelomereEnd()).rows(),strand_row(TELOTAG_CONSTRUCT),*feature_rows(TELOTAG_CONSTRUCT,prefix_width=4)],"The six splints cover all telomere-repeat phases. The 24-base barcode demultiplexes samples; ** marks terminal ligation."),("Select tagged ends on streptavidin and release at the adapter restriction site",[Row(chunks=[("biotin-TeloTag -> streptavidin capture -> wash -> restriction release",None,False)])],"Affinity selection removes most untagged genomic fragments; release retains the tagged telomere molecule."),("Attach Oxford Nanopore sequencing adapters",[Row(chunks=[("[motor-loaded ONT adapter] ** [enriched TeloTag–telomere fragment]",None,False)])],"The selected material enters the standard ligation library workflow; ** marks adapter ligation.")]
