#!/usr/bin/env python3
"""Build the sparse BD Rhapsody WTA chemistry page."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import rhapsody as R
import seqprimers as sp
from chemdraw import Construct,oligo,panel,protocol_graph_panels,strand_row
from page import head,info,table
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations
OUT=HERE.parent/"bd-rhapsody.html"

def render():
    final=R.final_library()
    graph=R.protocol_graph()
    r1=Construct([R._s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
                  R._s("UMI","[UMI]","umi",inferred=True,placeholder=True),
                  R._s("poly(T)","TTT",placeholder=True)])
    r2=Construct([R._s("cDNA","XXXXXXXX...",placeholder=True)])
    return "\n".join([head("BD Rhapsody WTA library chemistry"),'<div class="wrap" data-protocol-complete="true">',
      '<h1>BD Rhapsody WTA &mdash; microwell bead capture and random-primer second strand</h1>',
      info('Defining platform-use paper: <a href="https://doi.org/10.1038/s41592-018-0255-0">Yoon et al., <i>Nature Methods</i> (2019)</a>. Chemistry authority: BD Whole Transcriptome Analysis Alpha Protocol 23-21179-00 (December 2018).'),
      '<div class="caveat"><b>Workflow and sequence boundary.</b> This page follows BD\'s documented random-primer-extension workflow, not the older ligation route summarized in the paper. The checked BD alpha protocol prints the Randomer but does not resolve the bead, Universal Oligo or library-primer sequences; other relevant BD kit documents remain unchecked or unavailable. Dotted regions are architecture placeholders, not a claim that BD never published those bases.</div>',
      render_annotations(graph,AnnotationSlot.INTRO),
      '<h2>Key oligos</h2><seq>'+oligo("BD Randomer",R.randomer_segments())+
      oligo("bead capture oligo — published architecture",R.bead_capture_segments())+'</seq>',
      '<h2>Library construction</h2>',
      protocol_graph_panels(graph,cls="long"),
      '<h2>Read layout</h2>',
      sp.unavailable_diagram(final,'Primer-binding bases are unresolved in the checked BD sources'),
      panel([strand_row(r1,prefix="Read 1  ",suffix=""),strand_row(r2,prefix="Read 2  ",suffix="")],cls="small",caption="INFERRED — read-to-feature assignment follows the bead cartoon order and the Randomer TruSeq Read 2 handle; no checked source assigns reads to features. The platform paper used 101 × 2 sequencing (HiSeq 2500). Exact primer-binding bases remain unverified."),
      table(("Read","Reports"),[("Read 1","cell label, UMI, then poly(T)"),("Index 1","library index"),("Read 2","cDNA from the random-primed end")]),
      sp.read_spec_table([
          sp.read_spec("Read 1",101,"cell label, UMI, then poly(T)","extract"),
          sp.read_spec("Index 1 (i7)",None,"library index","demultiplex"),
          sp.read_spec("Read 2",101,"cDNA from the random-primed end","include")],
          "The run layout is source-declared because the checked sources do not establish the exact primer-binding bases."),'</div>'])

if __name__=="__main__":
    OUT.write_text(render(),encoding="utf-8");print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
