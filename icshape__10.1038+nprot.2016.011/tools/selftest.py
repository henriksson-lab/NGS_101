#!/usr/bin/env python3
"""Source-boundary checks for icSHAPE reaction identity and operation order."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from checks import Check, Source, have
from storage import sources_root
import chem

check = Check()
paper = Source(sources_root() / HERE.parent.name / 'nprot.2016.011_PMC4896316.xml.txt',
               'Flynn2016 original protocol PMC4896316, procedure26/40–50')
if have(check, paper, label='Original icSHAPE procedure'):
    text = paper.path.read_text()
    graph = chem.protocol_graph()
    order = list(graph.steps)
    check('published ligation enzyme is T4 RNA ligase1',
          'T4 RNA ligase 1, high concentration' in text and
          graph.materials['t4_rnl1'].label == 'T4 RNA ligase 1')
    check('source RT precedes capture and RNA digestion/cDNA release',
          text.index('Resuspend RNA samples from Step 39') <
          text.index('Add 10 µl of streptavidin bead solution') and
          order.index('reverse_transcription') < order.index('streptavidin_capture') <
          order.index('release_cdna') < order.index('circularize'))
check.report()
