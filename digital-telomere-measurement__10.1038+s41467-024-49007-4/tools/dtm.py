"""Digital Telomere Measurement workflow (Sanchez et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row, Scene, Segment, oligo
from telomere import TelomereEnd

SEQ_TETHER="AACCTTGGAGATGCACGGAGCAAGCAAT"
CAPTURE="TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC"
# Supplementary Data 3 prints the oligos unpartitioned.  The partition is derived
# within that table: all eight canonical capture oligos share their first 22 nt
# (the seqTether-complementary arm) and last 7 nt (CCTAACC); the seven TSQ variants
# share the first 46 nt with p-t3-nb01 and differ only in the last 7 nt.  So the
# barcode is 24 nt and the telomere-targeting arm 7 nt.
CAPTURE_TETHER_ARM=CAPTURE[:22]
CAPTURE_BARCODE=CAPTURE[22:46]
CAPTURE_TELOMERE_ARM=CAPTURE[46:]

def oligos():
    return (
        oligo("seqTether",(Segment("seqTether",SEQ_TETHER),)),
        oligo("p-t3-nb01 capture oligo",(
            Segment("seqTether-complementary arm",CAPTURE_TETHER_ARM),
            Segment("capture barcode",CAPTURE_BARCODE,"cbc"),
            Segment("telomere arm",CAPTURE_TELOMERE_ARM)),mods="/5Phos/"),
    )

def _duplex():
    sc=Scene()
    sc.strand("capture",[Segment("seqTether-complementary arm",CAPTURE_TETHER_ARM),
                          Segment("capture barcode",CAPTURE_BARCODE,"cbc"),
                          Segment("telomere arm",CAPTURE_TELOMERE_ARM)],
              label="p-t3-nb01 capture",mod5="p")
    sc.anneal("tether",[Segment("paired seqTether",SEQ_TETHER[:22]),
                         Segment("seqTether 3-prime tail",SEQ_TETHER[22:])],
              to="capture",pair=("paired seqTether","seqTether-complementary arm"),
              unpaired=("seqTether 3-prime tail",),label="seqTether")
    sc.mark("capture","capture barcode",f"{len(CAPTURE_BARCODE)}-nt capture barcode (variable region of Supplementary Data 3)")
    return sc

def _telomere_capture():
    end=TelomereEnd()
    sc=end.scene()
    sc.anneal("capture oligo",[Segment("capture tether arm",CAPTURE_TETHER_ARM),
                                Segment("capture barcode",CAPTURE_BARCODE,"cbc"),
                                Segment("CCTAACC telomere arm",CAPTURE_TELOMERE_ARM)],
              to="G-rich strand",pair=("CCTAACC telomere arm","G-rich 3-prime overhang"),
              shift=4,unpaired=("capture tether arm","capture barcode"),
              label="capture oligo",mod5="p")
    sc.anneal("seqTether",[Segment("paired seqTether",SEQ_TETHER[:22]),
                           Segment("seqTether 3-prime tail",SEQ_TETHER[22:])],
              to="capture oligo",pair=("paired seqTether","capture tether arm"),
              unpaired=("seqTether 3-prime tail",),label="seqTether")
    sc.footer("the 5′-phosphorylated capture end is distal (seqTether-paired); its 3′ end faces the C-rich 5′ terminus across a potential gap. Which nick T4 ligase seals is not stated in the Methods.",
              "C-rich strand")
    return sc

def panels():
    return [
        ("Anneal each barcoded capture oligo to seqTether",_duplex(),
         "The representative p-t3-nb01 sequence and seqTether are printed in Supplementary Data 3. The 22-nt tether-complementary arm, 24-nt barcode and 7-nt CCTAACC telomere arm are derived by comparing the printed canonical and variant capture oligos."),
        ("Anneal and ligate the capture duplex at native telomeres",_telomere_capture(),
         "CCTAACC is shown in one representative register on the terminal G-rich overhang; the native terminal phase is not drawn from the source. The Methods ligate the duplex overnight with T4 DNA ligase at 37 °C and heat-inactivate at 65 °C, but do not identify the sealed junction."),
        ("Fill the capture–telomere gap with Sulfolobus DNA polymerase IV",
         [Row(chunks=[("capture/tether ** gap-filled native telomere -- subtelomere",None,False)])],
         "Sulfolobus DNA polymerase IV (56 °C 2 min, 72 °C 15 min) fills potential gaps between the capture oligo and the double-strand/single-strand junction; ProNex beads (1.6×) then recover the post-capture library. The published core protocol does not mechanically fragment this product."),
        ("Ligate the Oxford Nanopore sequencing adapter",
         [Row(chunks=[("[motor-loaded ONT adapter] ** [gap-filled captured chromosome end]",None,False)])],
         "AMII (R9.4.1) or NA (R10.4) adapter with Quick T4 DNA ligase, 1 h at room temperature. The proprietary motor-adapter terminal bases are unavailable; the published capture and tether bases remain explicit."),
    ]
