"""Telobait enrichment for PacBio HiFi (Tham et al. 2023)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Row,Segment
from dumbbell import Dumbbell,dumbbell_rows
from telomere import AffinityCapture,PhasedTelomereAdapter,TelomereEnd
from restriction import RSAI,HINFI,ECORI

END=TelomereEnd(); F01="ATCGACGGTTCAAGATGCCAGATGCACGGAGCAGAATTCTTTTAC"
R_PHASES=("CCCTAA","ACCCTA","AACCCT","TAACCC","CTAACC","CCTAAC")
TELOBAIT=PhasedTelomereAdapter(F01,arm_copies=1)
CAPTURE=AffinityCapture("3′ biotin",ECORI.name)
INSERT=Construct([Segment("subtelomere","X"*30,placeholder=True),Segment("telomere","TTAGGG"*4,"r1"),Segment("telobait barcode + EcoRI","X"*24,"r3",placeholder=True)],name="captured telomere fragment")
SMRTBELL=Dumbbell(INSERT,"N"*16,"N"*16,name="telobait SMRTbell")

def panels(): return [("Trim subtelomeric DNA with RsaI and HinfI",[Row(chunks=[(f"subtelomere -- {RSAI.name}/{HINFI.name} cuts -- telomere -- native 3′ overhang",None,False)])],"Neither enzyme cuts the canonical TTAGGG repeat; the telomere end remains intact."),("Anneal and ligate one of six phased, barcoded 3′-biotin telobaits",TELOBAIT.annealed_scene(END),"The six ends cover all repeat phases. ** marks the terminal ligation; displayed F01 is one of the published barcoded sets."),("Capture on streptavidin, wash, and release with EcoRI",[Row(chunks=[("3′ biotin -> streptavidin capture -> wash -> EcoRI release",None,False)])],"The EcoRI site built into the telobait releases enriched telomere fragments."),("Build a PacBio HiFi SMRTbell",dumbbell_rows(SMRTBELL),"Hairpins close the selected fragment; ** marks each sealed insert–adapter junction.")]
