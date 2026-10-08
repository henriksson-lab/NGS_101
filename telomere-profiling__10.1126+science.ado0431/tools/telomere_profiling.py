"""Telomere Profiling enrichment workflow (Karimian et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from telomere import AffinityCapture,PhasedTelomereAdapter,TelomereEnd

TELOTAG=PhasedTelomereAdapter("N"*28)
CAPTURE=AffinityCapture("biotin","adapter-encoded restriction site")

def panels(): return [("Digest high-molecular-weight genomic DNA away from telomere ends",[Row(chunks=[("internal restriction cut | subtelomere -- telomere -- native 3′ overhang",None,False)])],"Restriction reduces non-telomeric length while leaving the native chromosome terminus available for tagging."),("Anneal six telomeric splints and ligate the biotinylated TeloTag",TELOTAG.annealed_scene(TelomereEnd()),"The six splints cover all telomere-repeat phases. ** marks the terminal ligation."),("Select tagged ends on streptavidin and release at the adapter restriction site",[Row(chunks=[("biotin-TeloTag -> streptavidin capture -> wash -> restriction release",None,False)])],"Affinity selection removes most untagged genomic fragments; release retains the tagged telomere molecule."),("Attach Oxford Nanopore sequencing adapters",[Row(chunks=[("[motor-loaded ONT adapter] ** [enriched TeloTag–telomere fragment]",None,False)])],"The selected material enters the standard ligation library workflow; ** marks adapter ligation.")]
