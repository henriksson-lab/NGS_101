#!/usr/bin/env python3
"""Original CAGEscan table transcription and source-branch regressions."""
import os
import re
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check, Source, have
from protocol_graph import Disposition, StepKind
import chem

check = Check()
bundle = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name
table = Source(bundle / "nmeth.1470_41592_2010_BFnmeth1470_MOESM315_ESM.pdf.txt", "Original Supplementary Table8, PDFp26")
if have(check, table, label="Original CAGEscan oligo table"):
    text = table.path.read_text().split("CAGEscan TS primer", 1)[1]
    values = re.findall(r"\n([ACGT]{10,}(?:XnrGrGrG|N6)?)\n", text)
    check("all eight displayed oligos agree with original Table8, including RNA rG notation",
          len(values) == 8 and all(value in shown for value, shown in zip(values, chem.CAGE_OLIGOS)))
    check("source RT/semi-reverse difference is preserved rather than silently normalized",
          values[1].startswith("TAGTCGAACGAAG") and "TAGTCGAACTGAAG" in values[3])
methods = Source(bundle / "nmeth.1470_supp_NIHMS207200-supplement-supplement_1.pdf.txt", "Original OnlineMethods CAGEscan section66–118")
if have(check, methods, label="Original CAGEscan Methods"):
    text = methods.path.read_text().split("CAGEscan reverse-transcription", 1)[1].split("Hep G2 total RNA", 1)[0]
    graph = chem.protocol_graph()
    check("both source pilot/production PCR pairs survive graph assembly",
          sum(s.kind is StepKind.PCR for s in graph.steps.values()), 2 * len(re.findall("quantitative PCR", text)))
    check("source RNAClean and both AMPure cleanups are retained",
          sum(s.kind is StepKind.CLEANUP for s in graph.steps.values()), text.count("RNAClean XP") + text.count("AMPure XP"))
    check("first-strand remainder receives an unresolved rather than invented physical fate",
          graph.terminals["remainder1"].disposition, Disposition.UNRESOLVED)
    notice = next(a for a in graph.annotations if a.title == "Apparent primer concentration unit issue in the source")
    check("source concentration-unit issue remains displayed without a replacement molarity",
          "100 mM" in text and "200 mM" in text and "100 mM" in notice.claims[0].claim and
          "200 mM" in notice.claims[1].claim and "without assigning a corrected" in notice.model_choice)
check.report()
