#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import r2c2 as R
from chemdraw import Scene, annotation_rows, panel, workflow_panel
from page import head, info

OUT = H.parent / "r2c2.html"

def render():
    lib = R.final_library()
    return "\n".join([
        head("R2C2 library chemistry"), '<div class="wrap">',
        '<h1>R2C2 &mdash; rolling circle to concatemeric consensus</h1>',
        '<p class="research-notes"><a href="01_r2c2.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1073/pnas.1806447115">Volden et al., <i>PNAS</i> (2018)</a>.'),
        info('A terminally matched DNA splint closes each full-length cDNA. Phi29 copies that circle repeatedly so one nanopore read contains several observations of the same molecule.'),
        '<h2>Reaction workflow</h2>', workflow_panel(R.workflow(), cls="long"),
        '<h2>Final nanopore library</h2>',
        panel(R.final_rows(), cls="long",
              caption="Debranched tandem repeats receive the paper's ONT 1D ligation adapters. Three repeat units are drawn; the molecular copy number varies. ** marks adapter ligation."),
        '<h2>Nanopore entry</h2>',
        info('There is no sequencing primer. The motor-bearing ONT 1D adapter engages the pore and feeds one intact concatemer strand through it; C3POa detects splint copies and computes the repeated-molecule consensus.'),
        '</div>'])

if __name__ == '__main__':
    OUT.write_text(render(), encoding='utf-8'); print(f'wrote {OUT}')
