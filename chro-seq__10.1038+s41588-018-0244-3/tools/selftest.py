#!/usr/bin/env python3
"""Check transcribed predecessor oligos and original ChRO-specific source boundaries."""
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
from checks import Check, Source, have
from protocol_graph import StepKind
import chem

check = Check()
base = data_root() / "sources"
bundle = base / HERE.parent.name
original = bundle / "PMC5502525_PMC5502525.xml.txt"
if not original.exists():
    original = base / "doi_10.1038+nprot.2016.086" / "nprot.2016.086_PMC5502525.xml.txt"
predecessor = Source(original, "Original Mahat2016 PRO-seq Table1 (separate from PRO-cap)")
if have(check, predecessor, label="Original PRO-seq Table1"):
    text = predecessor.path.read_text()
    table = text.split("Oligonucleotides required for PRO-seq and PRO-cap.", 1)[1]
    values = (
        re.search(r"VRA3\s+([ACGU]+-/Inverted dT/)", table)[1],
        re.search(r"VRA5\s+([ACGU]+)", table)[1],
        re.search(r"RP1\s+([ACGT]+)", table)[1],
        re.sub(r"\s+", "", re.search(r"RPI-n\s+([ACGTN\s]+)\s+DNA oligo", table)[1]),
    )
    for name, value, displayed in zip(("VRA3", "VRA5", "RP1", "RPI-n"), values, chem.OLIGOS):
        check("predecessor " + name + " agrees with original ordered sequence/modifications", value in displayed)
    check("VRA3 source end protection/phosphate remains explicit",
          "The 5′ end is phosphorylated" in table and "p" + values[0] in chem.OLIGOS[0])

article = Source(bundle / "s41588-018-0244-3_PMC6204104.xml.txt", "Original ChRO2018 Methods154–158")
if have(check, article, label="Original ChRO reaction order"):
    recipe = article.path.read_text().split("Chromatin Run-On and sequencing (ChRO-seq) library preparation.", 1)[1].split("Mapping ChRO-seq", 1)[0]
    source_captures = ("Nascent RNA was purified by binding streptavidin beads", "A second bead binding", "A third bead binding")
    graph = chem.protocol_graph()
    check("all three source captures survive graph assembly",
          sum(fragment in recipe for fragment in source_captures),
          sum(step.kind is StepKind.AFFINITY_CAPTURE for step in graph.steps.values()))

sheet = Source(bundle / "s41588-018-0244-3_41588_2018_244_MOESM3_ESM.xlsx.txt", "Original SupplementaryTable1 sample UMI column")
if have(check, sheet, label="Original sample-specific UMI assignments"):
    rows = sheet.path.read_text().split("## sheet Supplementary Table 2", 1)[0].splitlines()[2:]
    marked = {row.split('\t')[0] for row in rows if 'X' in row.split('\t')[6:7]}
    annotation = next(a.explanation for a in chem.protocol_graph().annotations if a.title == "Sample-specific six-base molecule identifiers")
    examples = ("UMU94-13_leChr", "GBM-88-04_Chr", "PDX-88-04_P57")
    check("public UMI examples are actually marked in original Table1",
          all(name in marked and name in annotation for name in examples))
check.report()
