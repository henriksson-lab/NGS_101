#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "indrop-v1__10.1016+j.cell.2015.04.044" / "tools"))
import page_parts
import indrop
import seqprimers as sp
from page import table
from indrop_v2_graph import TABLE1_CYCLES, READ_LIMITATION, protocol_graph
OUT = HERE.parent / "indrop-v2.html"


def v2_sequencing(protocol: str) -> str:
    """v2 read layout from the main text (custom primers, Table 1 and Step 166)."""
    if protocol != indrop.V2:
        raise ValueError(protocol)
    lib = indrop.library(protocol)
    if any(sp.locate(lib, p) is None for p in indrop.V2_PRIMERS):
        raise ValueError("inDrop read site missing")
    specs = [
        sp.read_spec("Read 1", "instrument-dependent (Table 1)", "transcript (sense cDNA)", "include"),
        sp.read_spec("Index 1 (i7)", "6 when multiplexed, else 0 (Table 1)",
                     "6-nt sample index from the PE1 primer", "demultiplex"),
        sp.read_spec("Read 2", "≥50 (Step 166); 51 or 55 in Table 1",
                     "barcode 1 (8–11), W1 (22), barcode 2 (8), UMI (6), then poly(T)", "extract"),
    ]
    rows = [(name, kit, str(r1), str(ix), str(r2)) for name, kit, r1, ix, r2 in TABLE1_CYCLES]
    return ('<h2>Read layout</h2>' + sp.diagram(lib, indrop.V2_PRIMERS)
            + sp.read_spec_table(specs, READ_LIMITATION + " Custom Read 1, Index and Read 2 "
                                 "primers are printed in the main-text Reagent Setup and "
                                 "Supplementary Table 2.")
            + '<h3>Published sequencing examples (main text Table 1)</h3>'
            + table(("Configuration", "Kit", "Read 1", "Index", "Read 2"), rows)
            + '</div>')


if __name__ == "__main__":
    page_parts.sequencing = v2_sequencing
    page_parts.render_page(indrop.V2, OUT, protocol_graph())
