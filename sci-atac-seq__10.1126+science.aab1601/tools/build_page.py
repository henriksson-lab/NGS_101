#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page, protocol_graph as sci_atac_graph
import sciatac

OUT = Path(__file__).resolve().parents[1] / "sci-atac-seq.html"


def protocol_graph():
    return sci_atac_graph(sciatac.SCI15)


if __name__ == "__main__":
    render_page(sciatac.SCI15, OUT, graph_factory=protocol_graph)

