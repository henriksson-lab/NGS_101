"""Construct model for the individual dscATAC-seq and dsciATAC-seq pages."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
from chemdraw import Construct, Scene, Segment

C1="TATGCATGAC"; C2="AGTCACTGAG"; BC_NT=7; TN5_BC_NT=6
REP_TN5_BC="AAAGAA"

def _s(n,t,g=None,**kw): return Segment(n,t,g,**kw)

def bead_published_segments():
    return [_s("BC1","N"*BC_NT,"cbc",placeholder=True),
            _s("phase block","NN",placeholder=True),_s("constant 1",C1),
            _s("BC2","N"*BC_NT,"cbc",placeholder=True),_s("constant 2",C2),
            _s("BC3","N"*BC_NT,"cbc",placeholder=True),_s("s5 primer",nx.S5,"s5")]

def tn5_r1(protocol):
    if protocol=="dscATAC-seq": return [_s("s5",nx.S5,"s5"),_s("ME",nx.ME,"me")]
    if protocol=="dsciATAC-seq": return [_s("s5",nx.S5,"s5"),
        _s("Tn5 barcode",REP_TN5_BC,"cbc"),_s("ME",nx.ME,"me")]
    raise ValueError(protocol)

def bead_priming_scene():
    template=[_s("s5'",nx.S5_RC,"s5"),_s("ME'",nx.ME_RC,"me"),
              _s("insert","XXXXXXXX...",placeholder=True)]
    bead=[_s("unpublished bead 5' arm","XXXXXXXXXXXX",inferred=True,placeholder=True),
          *bead_published_segments()]
    sc=Scene(); sc.strand("tagmented strand",template,label="fragment")
    sc.anneal("released bead oligo",bead,to="tagmented strand",
              pair=("s5 primer","s5'"),label="bead primer",
              unpaired=("unpublished bead 5' arm","BC1","phase block","constant 1",
                        "BC2","constant 2","BC3"))
    sc.arrow("released bead oligo","polymerase")
    return sc

def filled_fragment(protocol):
    return Construct([*tn5_r1(protocol),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
                      _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7")],
                     name=f"{protocol} gap-filled fragment")

def final_library(protocol):
    inner=[*bead_published_segments()]
    if protocol=="dsciATAC-seq": inner.append(_s("Tn5 barcode", "N"*TN5_BC_NT,"cbc",placeholder=True))
    elif protocol!="dscATAC-seq": raise ValueError(protocol)
    inner += [_s("ME",nx.ME,"me"),_s("genomic insert","XXXXXXXX...XXXXXXXX",placeholder=True),
              _s("ME'",nx.ME_RC,"me"),_s("s7'",nx.S7_RC,"s7"),
              _s("i7'","N"*8,"cbc",placeholder=True),_s("P7'",il.P7_RC,"p7")]
    return Construct([_s("unpublished bead/P5 arm","[BEAD RELEASE / P5 / R1 ARM]",inferred=True,
                         placeholder=True),*inner],name=f"{protocol} final library")

def _validate():
    if "".join(s.top for s in tn5_r1("dsciATAC-seq")) != nx.S5+REP_TN5_BC+nx.ME:
        raise ValueError("barcoded Tn5 adaptor no longer matches Supplementary Table 6")
    bead_priming_scene().rows()
_validate()
