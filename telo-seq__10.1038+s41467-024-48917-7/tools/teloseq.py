"""Telo-seq targeted Oxford Nanopore workflow (Schmidt et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from telomere import PhasedTelomereAdapter,TelomereEnd
from restriction import ECORV

CORE="AGCAATACGTAACTGAACGAAGT"
TELORETTES=(
"AGCAATACGTAACTGAACGAAGTCCCTAACCCTAACCCTAA",
"AGCAATACGTAACTGAACGAAGTTAACCCTAACCCTAACCC",
"AGCAATACGTAACTGAACGAAGTCTAACCCTAACCCTAACC",
"AGCAATACGTAACTGAACGAAGTCCTAACCCTAACCCTAAC",
"AGCAATACGTAACTGAACGAAGTAACCCTAACCCTAACCCT",
"AGCAATACGTAACTGAACGAAGTACCCTAACCCTAACCCTA")
S1="ACTTCGTTCAGTTACGTATTGCTAGCAAT"
ADAPTER=PhasedTelomereAdapter(CORE)

def panels(): return [("Anneal six phased telorettes and ligate the native telomere terminus",ADAPTER.annealed_scene(TelomereEnd()),"All oligos are 5′ phosphorylated. ** marks ligation to the native C-rich strand."),("EcoRV cuts distal genomic DNA; Klenow exo− dA-tails only those new blunt ends",[Row(chunks=[(f"telorette ** telomere -- subtelomere -- {ECORV.name}|A     (internal cut end only)",None,False)])],"The terminal telorette end is preserved while internal EcoRV ends receive dA, suppressing wrong-end adapter ligation."),("S1 splint recruits the Oxford Nanopore sequencing adapter",[Row(chunks=[("[motor-loaded ONT adapter] ** [S1:telorette core] ** telomere -- subtelomere -- EcoRV|A",None,False)])],"S1 bridges the disclosed telorette core to the AMII adapter. ** marks the ligated boundaries.")]
