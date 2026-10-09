"""Digital Telomere Measurement workflow (Sanchez et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Row,Segment,feature,feature_rows,strand_row
from telomere import PhasedTelomereAdapter,TelomereEnd

# The main article defines the duplex architecture but places ordered bases in its
# supplementary oligo table. Unknown bases stay N rather than becoming guessed sequence.
CAPTURE_BARCODE=Segment(
    "24-nt ONT capture barcode", "N"*24, "cbc", placeholder=True,
    feature=feature("capture_sample_index", "sample_index", "whitelist",
                    whitelist="Oxford Nanopore native barcode set"))
CAPTURE=PhasedTelomereAdapter(CAPTURE_BARCODE.top,
                              core_segments=(CAPTURE_BARCODE,))
CAPTURE_CORE=Construct([CAPTURE_BARCODE],name="DTM capture barcode")

def panels(): return [("Duplex each barcoded telomere-capture oligo with its sequencing tether",[strand_row(CAPTURE_CORE),*feature_rows(CAPTURE_CORE,prefix_width=4),Row(chunks=[("3′-[sequencing tether]",None,False)])],"The 24-nt ONT barcode identifies the sample. The tether is pre-annealed before genomic DNA is added; undisclosed bases remain explicit N."),("Anneal the phased capture end and ligate it to native telomeres",CAPTURE.annealed_scene(TelomereEnd()),"A capture family recognizes every human telomere-repeat phase. ** marks the terminal ligation."),("Fill any gap with Sulfolobus DNA polymerase IV, then mechanically fragment",[Row(chunks=[("capture ** filled telomere -- subtelomere -- | random mechanical break",None,False)])],"Gap fill completes the tagged end before random fragmentation."),("Ligate the Oxford Nanopore sequencing adapter to the prepared library",[Row(chunks=[("[motor-loaded ONT adapter] ** [prepared telomere fragment]",None,False)])],"The completed capture construct is converted to a nanopore library; ** marks adapter ligation.")]
