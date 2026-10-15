#!/usr/bin/env python3
from pathlib import Path
from page_parts import render_page
import indrop
OUT = Path(__file__).resolve().parents[1] / "indrop-v1.html"


def protocol_graph():
    return indrop.protocol_graph(indrop.V1)


if __name__ == "__main__": render_page(indrop.V1, OUT, graph=protocol_graph())

