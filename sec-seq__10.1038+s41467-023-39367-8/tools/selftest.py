#!/usr/bin/env python3
"""SEC-seq source transcription checks; absent source originals skip cleanly."""
from pathlib import Path
import os
import re
import sys
H = Path(__file__).resolve().parent
ROOT = H.parents[1]
sys.path[:0] = [str(H), str(ROOT / "lib")]
import sec_seq as M
from checks import Check, Source, have

c = Check()
bundle = Path(os.environ.get("CHEM_DATA", ROOT / "_data")) / "sources" / "sec-seq__10.1101+2022.08.25.505190"
guide = Source(bundle / "10x_CG000317_ChromiumNextGEMSingleCell3_v3.1_CellSurfaceProtein_RevE.pdf.txt", M.GUIDE)
article = Source(bundle / "s41467-023-39367-8_PMC10272111.xml.txt", M.SOURCE)
if have(c, guide, label="CG000317 complete final-product transcription"):
    text = re.sub(r"\s+", "", guide.path.read_text())
    products = re.findall(r"5[′’']-([A-Za-z0-9_-]+)-3[′’']", text)
    expanded = []
    for product in products:
        product = product.replace("cDNA_Insert", "X" * 40)
        product = re.sub(r"N(\d+)", lambda m: "N" * int(m[1]), product)
        expanded.append(product.replace("-", ""))
    c("GEX endpoint matches a complete original guide product", M.RNA.top().replace("B", "N") in expanded)
    c("secretion endpoint matches a complete original guide product", M.SEC.top().replace("B", "N") in expanded)
if have(c, article, label="published human ASC read lengths"):
    text = " ".join(article.path.read_text().split())
    match = re.search(r"(\d+) bp Read1, (\d+) bp i7 Index, (\d+) bp i5 Index, and (\d+) bp Read2", text)
    c("published complete four-read run located", match is not None)
    if match:
        run = tuple(map(int, match.groups()))
        for heading, *_ in M.FINAL_LIBRARIES:
            c(heading + " follows the source run", tuple(M.READ_LENGTHS[heading][role] for role in M.ROLES), run)
c.report()
