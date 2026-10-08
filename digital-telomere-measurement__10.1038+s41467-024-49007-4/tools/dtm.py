"""Digital Telomere Measurement workflow (Sanchez et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from telomere import PhasedTelomereAdapter,TelomereEnd

# The main article defines the duplex architecture but places ordered bases in its
# supplementary oligo table. Unknown bases stay N rather than becoming guessed sequence.
CAPTURE=PhasedTelomereAdapter("N"*24)

def panels(): return [("Duplex each barcoded telomere-capture oligo with its sequencing tether",[Row(chunks=[("5′-[barcoded capture oligo]-[telomere arm]",None,False)]),Row(chunks=[("3′-[sequencing tether]",None,False)])],"The tether is pre-annealed before genomic DNA is added; undisclosed here means explicit N, not invented bases."),("Anneal the phased capture end and ligate it to native telomeres",CAPTURE.annealed_scene(TelomereEnd()),"A capture family recognizes every human telomere-repeat phase. ** marks the terminal ligation."),("Fill any gap with Sulfolobus DNA polymerase IV, then mechanically fragment",[Row(chunks=[("capture ** filled telomere -- subtelomere -- | random mechanical break",None,False)])],"Gap fill completes the tagged end before random fragmentation."),("Ligate the Oxford Nanopore sequencing adapter to the prepared library",[Row(chunks=[("[motor-loaded ONT adapter] ** [prepared telomere fragment]",None,False)])],"The completed capture construct is converted to a nanopore library; ** marks adapter ligation.")]
