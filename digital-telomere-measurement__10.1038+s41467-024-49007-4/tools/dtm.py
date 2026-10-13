"""Digital Telomere Measurement workflow (Sanchez et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row, Scene, Segment, oligo
from telomere import TelomereEnd

SEQ_TETHER="AACCTTGGAGATGCACGGAGCAAGCAAT"
CAPTURE="TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC"
CAPTURE_TETHER_ARM=CAPTURE[:22]
CAPTURE_BARCODE=CAPTURE[22:-6]
CAPTURE_TELOMERE_ARM=CAPTURE[-6:]

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
    sc.mark("capture","capture barcode","published 25-nt capture barcode")
    return sc

def _telomere_capture():
    end=TelomereEnd()
    sc=end.scene()
    sc.anneal("capture oligo",[Segment("capture/tether body",CAPTURE[:-6],"cbc"),
                                Segment("CTAACC telomere arm",CAPTURE_TELOMERE_ARM)],
              to="G-rich strand",pair=("CTAACC telomere arm","G-rich 3-prime overhang"),
              shift=4,unpaired=("capture/tether body",),
              label="capture duplex",mod5="p")
    sc.ligation_end("capture oligo","capture/tether body",kind="cohesive",
                    polarity="5-prime",owner="adapter",text="terminal capture ligation")
    sc.footer("** ligase joins the capture 5′ phosphate at the native C-rich terminus",
              "capture oligo","capture/tether body")
    return sc

def panels():
    return [
        ("Anneal each barcoded capture oligo to seqTether",_duplex(),
         "The representative p-t3-nb01 sequence is printed in Supplementary Data 5: a 22-nt tether-complementary arm, 25-nt barcode and 6-nt telomere arm."),
        ("Anneal and ligate the capture duplex at native telomeres",_telomere_capture(),
         "CTAACC is shown in one representative register on the terminal G-rich overhang; the capture oligo is covalently sealed at the chromosome end."),
        ("Fill the capture–telomere gap with Sulfolobus DNA polymerase IV",
         [Row(chunks=[("capture/tether ** gap-filled native telomere -- subtelomere",None,False)])],
         "Polymerase IV completes the tagged telomere substrate. The published core protocol does not mechanically fragment this product."),
        ("Ligate the Oxford Nanopore sequencing adapter",
         [Row(chunks=[("[motor-loaded ONT adapter] ** [gap-filled captured chromosome end]",None,False)])],
         "The proprietary motor-adapter terminal bases are unavailable; the published capture and tether bases remain explicit."),
    ]
