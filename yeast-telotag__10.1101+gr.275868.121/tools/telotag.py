"""Yeast TeloTag nanopore workflow (Sholes et al. 2022)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Scene,Segment,complement_segments,Row
from telomere import TelomereEnd

END=TelomereEnd(repeat="TGGTGT",duplex_copies=3,overhang_copies=2)
POLYA="A"*18; TELOTAG_CORE="N"*24

def tail_scene():
    top=[Segment("subtelomere","X"*16,placeholder=True),Segment("telomere",END.duplex,"r1"),Segment("native overhang",END.overhang,"r2"),Segment("poly(A) tail",POLYA,"r3")]
    sc=Scene(); sc.strand("G-rich strand",top,label="chromosome")
    sc.anneal("C-rich strand",complement_segments(top[:2]),to="G-rich strand",pair=("telomere'","telomere"),label="chromosome")
    tag=[Segment("TeloTag core",TELOTAG_CORE,"r3",placeholder=True),Segment("oligo(dT)","T"*18,"r2")]
    sc.anneal("TeloTag",tag,to="G-rich strand",pair=("oligo(dT)","poly(A) tail"),label="TeloTag",unpaired=("TeloTag core",))
    sc.footer("** T4 ligase seals the filled TeloTag–chromosome nick","TeloTag","TeloTag core")
    return sc

def panels(): return [("Terminal transferase adds poly(A) at the native 3′ chromosome end",tail_scene(),"The poly(A) tail preserves a recognizable terminal landmark; the paper uses five TeloTag primers."),("Sulfolobus DNA polymerase IV fills the end and T4 ligase seals the nick",[Row(chunks=[("[duplex chromosome][telomere][poly(A):oligo(dT)][TeloTag] **",None,False)])],"Fill-in makes a blunt tagged end; ** identifies the sealing ligation."),("Attach Oxford Nanopore library adapters",[Row(chunks=[("[motor-loaded ONT adapter] ** [tagged chromosome]",None,False)])],"The manufacturer’s nanopore adapter is ligated after TeloTag completion.")]
