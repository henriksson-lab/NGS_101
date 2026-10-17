#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import indrop
OUT = Path(__file__).resolve().parents[1] / "indrop-v1.html"
SOURCE_OLIGO_EXCEPTIONS = {
    "pcr_forward": "Cell 2015 Supplemental Information Table S5 wraps PCR enrichment primer 1 "
                   "across its 'primers:' row label (…CCCTA / primers: / CACGA), so the text "
                   "twin never holds it contiguously; selftest.py rejoins the two printed "
                   "lines and compares the whole primer.",
}


def protocol_graph():
    return indrop.protocol_graph(indrop.V1)


if __name__ == "__main__": render_page(indrop.V1, OUT, graph=protocol_graph())

