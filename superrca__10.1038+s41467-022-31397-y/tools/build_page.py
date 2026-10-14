#!/usr/bin/env python3
"""Build the superRCA assay schematic."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import superrca as M
import seqprimers as sp
from chemdraw import workflow_panel
from page import head, info

OUT = HERE.parent / "superrca.html"


def render() -> str:
    return "\n".join([
        head(M.TITLE), '<div class="wrap">',
        f'<h1>{M.TITLE}</h1>',
        '<p class="research-notes"><a href="01_superrca.html">Research notes</a></p>',
        info(M.SOURCE),
        info("This is a role-level assay schematic. The article and available supplement "
             "establish the reaction order but do not supply the exact oligo sequences. "
             "The peer-review response mentions a missing Supplementary Table S2. "
             "Bracketed regions have no asserted sequence or length."),
        '<h2>Molecular workflow</h2>', workflow_panel(M.workflow(), cls="long"),
        '<h2>Detection output</h2>',
        info("Fluorescent probes distinguish mutant and wild-type superRCA objects. "
             "A flow cytometer counts those objects. The NGS and ddPCR methods in "
             "the paper are independent comparisons."),
        sp.read_spec_table([
            sp.read_spec("Sequencing output", None,
                         "No sequencing library or nucleotide reads", "not_applicable")
        ], "superRCA is quantified by flow cytometry or fluorescence microscopy."),
        '</div>',
    ])


if __name__ == "__main__":
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}")
