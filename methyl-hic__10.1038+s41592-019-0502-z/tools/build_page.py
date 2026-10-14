#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import methyl_hic as M
from adapter_ligation import AdapterLigationGeometry
from chemdraw import panel
from single_cell_hic_page import render
OUT=H.parent/"methyl-hic.html"
if __name__ == "__main__":
    page = render(M)
    limitation = panel(AdapterLigationGeometry.unavailable(
        "The supplementary protocol names indexed Illumina adapters and TruSeq PCR primers, but does not print the ligated adapter strands or their ligation-ready end geometry.").rows(),
        cls="small", caption="Source limitation before on-bead adapter ligation.")
    page = page.replace('<h2>Final sequencing molecule</h2>',
                        '<h2>Sequencing-adapter end</h2>' + limitation +
                        '<h2>Final sequencing molecule</h2>')
    OUT.write_text(page, encoding="utf-8"); print(f"wrote {OUT}")
