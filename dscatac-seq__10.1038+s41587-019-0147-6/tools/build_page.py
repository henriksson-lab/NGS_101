#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import dscatac
HERE=Path(__file__).resolve().parent; OUT=HERE.parent/"dscatac-seq.html"


def protocol_graph():
    return dscatac.protocol_graph("dscATAC-seq")


if __name__=="__main__": OUT.write_text(render_page("dscATAC-seq", graph_factory=protocol_graph),encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
