"""Compatibility entry point for the original2014MARS page."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import mars_graph as P
from batch_page import render

def render_page(protocol="MARS-seq",**kwargs):
 if protocol!="MARS-seq": raise ValueError("2019MARS-seq2.0 has its own native renderer")
 return render(P)
