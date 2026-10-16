#!/usr/bin/env python3
"""Original2008 SOM source-boundary regressions; no later recipe substitution."""
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
from checks import Check, Source, have
from protocol_graph import InputRole, StepKind
import chem

check = Check()
root = data_root() / "sources"
source = Source(root / HERE.parent.name / "594009_9144e1ea58ad4573871ec024827290f3.pdf.txt",
                "Original Core2008 Science1162228 Supporting Online Material")
if have(check, source, label="Original2008 SOM"):
    text = source.path.read_text()
    graph = chem.protocol_graph()
    check("recovered source is the original2008 SOM, not the later implementation",
          "DOI: 10.1126/science.1162228" in text and "Published 4 December 2008" in text)
    lanes = text.split("Lanes: 1) Input", 1)[1].split("Figure S16", 1)[0]
    count = len(set(re.findall(r"Elution[- ](\d)", lanes)))
    check("every original FigureS15 enrichment/elution is represented",
          sum(s.kind is StepKind.AFFINITY_CAPTURE for s in graph.steps.values()), count)
    steps = list(graph.steps)
    check("source FigureS1 5′-then-3′ adapter order remains intact",
          "5’ end follwed by the 3’end" in text and steps.index("ligate5") < steps.index("ligate3"))
    check("each adapter step precedes its own source antibody capture",
          steps.index("ligate5") < steps.index("second_capture") < steps.index("ligate3") < steps.index("third_capture"))
    check("two original DNase treatments remain on the molecular route",
          "another digestion with RNase-free DNaseI" in text and
          steps.index("digest_dna_first") < steps.index("hydrolyse") < steps.index("digest_dna_second"))
    check("source RNA digestion precedes PCR rather than becoming a physical discard",
          "RNAse cocktail" in text and "RNAse H" in text and
          any(i.material == "pcr_template" and i.role is InputRole.TEMPLATE for i in graph.steps["pcr"].inputs))
    reads = int(re.search(r"All reads were (\d+) bases long", text)[1])
    check("original read length retained without selecting a modern primer site",
          chem.READ_LENGTHS == {"Read 1": reads} and
          chem.READ_SPECS[0].cycles == str(reads) and not chem.SEQ_PRIMERS)
check.report()
