"""Construct model for the documented BD Rhapsody WTA alpha workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Scene,Segment

# Printed by BD in Whole Transcriptome Analysis Alpha Protocol 23-21179-00 (12/2018).
RANDOMER_HANDLE="TCAGACGTGTGCTCTTCCGATCT"
RANDOM_NT=9
DRAWN_POLYT=15  # visual tract only; BD does not publish the bead oligo's exact sequence

def _s(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def bead_capture_segments():
    """Vendor-published architecture, with proprietary bases represented symbolically."""
    return [_s("bead universal region","XXXXXXXXXXXX",inferred=True,placeholder=True),
            _s("cell label","NNNNNNNNN","cbc",inferred=True,placeholder=True),
            _s("UMI","NNNNNNNN","umi",inferred=True,placeholder=True),
            _s("poly(T)","T"*DRAWN_POLYT,inferred=True)]

def randomer_segments():
    return [_s("published R2-side handle",RANDOMER_HANDLE,"r2"),
            _s("random N9","N"*RANDOM_NT,placeholder=True)]

def capture_scene():
    mrna=[_s("mRNA body","XXXXXXXXXXXX",placeholder=True),
          _s("poly(A)","A"*DRAWN_POLYT)]
    sc=Scene();sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead oligo",bead_capture_segments(),to="mRNA",
              pair=("poly(T)","poly(A)"),label="bead oligo")
    sc.arrow("bead oligo","reverse transcriptase")
    return sc

def bead_cdna():
    return Construct([*bead_capture_segments(),
                      _s("first-strand cDNA","XXXXXXXX...XXXXXXXX",placeholder=True)],
                     name="bead-bound first-strand cDNA")

def random_priming_scene():
    template=[_s("bead-proximal cDNA","XXXXXXXX...XXXXXXXX",placeholder=True),
              _s("random site","XXXXXXXXX",placeholder=True)]
    sc=Scene();sc.strand("first strand",template,label="bead-bound cDNA")
    sc.anneal("Randomer",randomer_segments(),to="first strand",
              pair=("random N9","random site"),label="Randomer")
    sc.arrow("Randomer","Klenow exo-")
    return sc

def rpe_product():
    return Construct([*randomer_segments()[:-1],
                      _s("random-primed cDNA","XXXXXXXX...XXXXXXXX",placeholder=True),
                      _s("poly(A) copy","AAA",placeholder=True),
                      _s("UMI'","N"*8,"umi",inferred=True,placeholder=True),
                      _s("cell label'","N"*9,"cbc",inferred=True,placeholder=True),
                      _s("bead universal complement","XXXXXXXXXXXX",inferred=True,
                         placeholder=True)],name="released RPE product")

def final_library():
    return Construct([
        _s("unpublished P5/Read-1 arm","[P5 / READ 1 ARM]",inferred=True,placeholder=True),
        _s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
        _s("UMI","[UMI]","umi",inferred=True,placeholder=True),
        _s("poly(T)/cDNA","TTTXXXXXXXX...XXXXXXXX",placeholder=True),
        _s("published Randomer handle",RANDOMER_HANDLE,"r2"),
        _s("unpublished index/P7 arm","[i7 / P7 ARM]",inferred=True,placeholder=True),
    ],name="BD Rhapsody WTA library boundary")

def _validate():
    if "".join(s.top for s in randomer_segments()) != RANDOMER_HANDLE+"N"*RANDOM_NT:
        raise ValueError("Randomer no longer matches the BD alpha protocol")
    capture_scene().rows();random_priming_scene().rows()
_validate()
