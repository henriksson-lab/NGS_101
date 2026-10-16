#!/usr/bin/env python3
"""Source boundaries and workflow regressions, without selecting an unknown oligo kit."""
import os
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check, Source, have
from protocol_graph import InputRole
import chem

check = Check()
data = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name
article = Source(data / "s13059-017-1247-6_PMC5488340.xml.txt", "Ponnaluri2017 original Methods, NicE-seq library construction")
if have(check, article, label="Original NicE-seq library method"):
    text = article.path.read_text()
    recipe = text.split("NicE-seq library construction\n", 1)[1].split("ATAC-seq library construction", 1)[0]
    graph = chem.protocol_graph()
    capture_input = next(i.material for i in graph.steps["capture"].inputs if i.role is InputRole.TEMPLATE)
    check("source direct ligation-to-capture route has no inserted kit cleanup",
          "Without further purification" in recipe and
          capture_input in graph.steps["adapter_ligation"].outputs)
    pcr_input = next(i.material for i in graph.steps["library_pcr"].inputs if i.role is InputRole.TEMPLATE)
    check("source bead-bound PCR follows capture rather than a free-library substitution",
          "beads were resuspended" in recipe.lower() and
          pcr_input in graph.steps["sample_beads"].outputs and
          any(i.material in graph.steps["capture"].outputs for i in graph.steps["sample_beads"].inputs))
    resuspension, aliquot = re.search(r"resuspended in (\d+)\s*μL.*?and (\d+)\s*μL was used", recipe).groups()
    check("source aliquot volumes and unresolved remainder are preserved",
          (chem.BEAD_RESUSPENSION_UL, chem.PCR_ALIQUOT_UL, chem.BEAD_REMAINDER_UL,
           graph.terminals["bead_remainder"].disposition.value),
          (int(resuspension), int(aliquot), int(resuspension) - int(aliquot), "unresolved"))


geo = Source(data / "geo-GSE97889-all-full.txt", "Original GSE97889 family sample extraction metadata")
if have(check, article, geo, label="Original article/deposit differences"):
    a = article.path.read_text()
    g = geo.path.read_text()
    pairs = (
        (re.search(r"(\d+) U of DNA polymerase I", a)[1], re.search(r"(\d+) U of DNA polymerase I", g)[1]),
        ("overnight", "1 hour"),
        (re.search(r"B&W buffer plus ([\d.]+%)", a)[1], re.search(r"B&W buffer plus ([\d.]+%)", g)[1]),
    )
    check("page retains all three independently located article/deposit differences",
          all(value in chem.CAVEAT for pair in pairs for value in pair))

check.section("Unassigned study-specific oligo identity")
library, primers = chem.nice_library()
check("unassigned oligo identity cannot introduce canonical adapter bases or index cycles",
      all(segment.is_role_token() for segment in library) and not primers and chem.READ_LENGTHS is None)
check.report()
