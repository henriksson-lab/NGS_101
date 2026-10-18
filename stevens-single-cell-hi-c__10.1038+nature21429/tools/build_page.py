#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import stevens as M
from batch_page import render
from chemdraw import protocol_graph_panels
OUT = HERE.parent / "stevens-single-cell-hi-c.html"
page=render(M)
page=page.rsplit('</div>',1)[0]+'<h2>Published alternative: Nextera XT</h2><p>This is an alternative route from captured contacts, not a second physical aliquot of the CAA-adapter preparation. Native i5/i7 widths and acquisition cycles remain unassigned.</p>'+protocol_graph_panels(M.nextera_graph(),cls="long")+'</div>'
OUT.write_text(page, encoding="utf-8")
print(f"wrote {OUT}")
