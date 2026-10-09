"""Yeast TeloTag nanopore workflow (Sholes et al. 2022)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import (Construct,Row,Scene,Segment,complement_segments,feature,
                      feature_rows)
from telomere import TelomereEnd

END=TelomereEnd(repeat="TGGTGT",duplex_copies=3,overhang_copies=2)
POLYA="A"*18; TELOTAG_CORE="N"*24
SAMPLE_INDEX=feature("ont_native_sample_index", "sample_index", "whitelist",
                     whitelist="Oxford Nanopore EXP-NBD104 native barcode set")

def tail_scene():
    top=[Segment("subtelomere","X"*16,placeholder=True),Segment("telomere",END.duplex,"r1"),Segment("native overhang",END.overhang,"r2"),Segment("poly(A) tail",POLYA,"r3")]
    sc=Scene(); sc.strand("G-rich strand",top,label="chromosome")
    sc.anneal("C-rich strand",complement_segments(top[:2]),to="G-rich strand",pair=("telomere'","telomere"),label="chromosome")
    tag=[Segment("TeloTag core",TELOTAG_CORE,"r3",placeholder=True),Segment("oligo(dT)","T"*18,"r2")]
    sc.anneal("TeloTag",tag,to="G-rich strand",pair=("oligo(dT)","poly(A) tail"),label="TeloTag",unpaired=("TeloTag core",))
    sc.footer("** T4 ligase seals the filled TeloTag–chromosome nick","TeloTag","TeloTag core")
    return sc

def nanopore_library() -> Construct:
    return Construct([
        Segment("motor-loaded ONT adapter", "[motor adapter]", "r2", placeholder=True),
        Segment("left ONT native barcode", "[ONT barcode]", "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
        Segment("TeloTag-labelled chromosome", "[tagged chromosome]", placeholder=True),
        Segment("right ONT native barcode", "[ONT barcode]", "cbc", placeholder=True,
                feature=SAMPLE_INDEX),
    ], name="barcoded TeloTag nanopore library")

def nanopore_library_rows():
    lib=nanopore_library(); sc=Scene(); sc.strand("library",list(lib),label="nanopore library")
    sc.junction("library","motor-loaded ONT adapter","left ONT native barcode","adapter assembly")
    prefix=len("nanopore library") + 1 + len("5'- ")
    return [*sc.rows(),*feature_rows(lib,prefix_width=prefix)]

def panels():
    return [
        ("Terminal transferase adds poly(A) at the native 3′ chromosome end",
         tail_scene(),
         "The poly(A) tail preserves a recognizable terminal landmark; the paper uses five TeloTag primers."),
        ("Sulfolobus DNA polymerase IV fills the end and T4 ligase seals the nick",
         [Row(chunks=[("[duplex chromosome][telomere][poly(A):oligo(dT)][TeloTag] **",None,False)])],
         "Fill-in makes a blunt tagged end; ** identifies the sealing ligation."),
        ("Attach native barcodes and the Oxford Nanopore sequencing adapter",
         nanopore_library_rows(),
         "EXP-NBD104 supplies the sample barcode used at both molecule ends; the motor-loaded sequencing adapter is assembled after barcoding."),
    ]
