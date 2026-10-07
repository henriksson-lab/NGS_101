#!/usr/bin/env python3
"""Build the sparse BD Rhapsody WTA chemistry page."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import rhapsody as R
from chemdraw import Construct,annotation_rows,oligo,panel,strand_row
from page import head,info,table
OUT=HERE.parent/"bd-rhapsody.html"

def render():
    cdna=R.bead_cdna();rpe=R.rpe_product();final=R.final_library()
    r1=Construct([R._s("cell label","[CELL LABEL]","cbc",inferred=True,placeholder=True),
                  R._s("UMI","[UMI]","umi",inferred=True,placeholder=True),
                  R._s("poly(T)","TTT",placeholder=True)])
    r2=Construct([R._s("cDNA","XXXXXXXX...",placeholder=True)])
    return "\n".join([head("BD Rhapsody WTA library chemistry"),'<div class="wrap">',
      '<h1>BD Rhapsody WTA &mdash; microwell bead capture and random-primer second strand</h1>',
      info('Defining platform-use paper: <a href="https://doi.org/10.1038/s41592-018-0255-0">Yoon et al., <i>Nature Methods</i> (2019)</a>. Chemistry authority: BD Whole Transcriptome Analysis Alpha Protocol 23-21179-00 (December 2018).'),
      '<div class="caveat"><b>Workflow and sequence boundary.</b> This page follows BD\'s documented random-primer-extension workflow, not the older ligation route summarized in the paper. BD publishes the Randomer but not the bead, Universal Oligo or library-primer sequences; dotted regions are architecture placeholders.</div>',
      '<h2>Key oligos</h2><seq>'+oligo("BD Randomer",R.randomer_segments())+
      oligo("bead capture oligo — published architecture",R.bead_capture_segments())+'</seq>',
      '<h2>Library construction</h2>',
      '<h3>(1) Capture poly(A) RNA in microwells, retrieve beads, reverse transcribe</h3>',
      panel(R.capture_scene().rows(),cls="long",caption="The cell label and UMI are installed at the transcript's 3' end."),
      panel([strand_row(cdna),*annotation_rows(cdna)],cls="long"),
      '<h3>(2) Exonuclease I removes unused capture oligos; RNA is removed</h3>',
      '<h3>(3) Randomer anneals along the bead-bound first strand</h3>',
      panel(R.random_priming_scene().rows(),cls="long",caption="The published N9 randomer carries the fixed R2-side handle."),
      '<h3>(4) Klenow exo- extends toward the bead; heat releases the new strand</h3>',
      panel([strand_row(rpe),*annotation_rows(rpe)],cls="long",caption="Only products that reach the bead-proximal label contain cell identity and UMI."),
      '<h3>(5) RPE PCR and indexed library PCR</h3>',
      panel([strand_row(final),*annotation_rows(final)],cls="long",caption="INFERRED — supported inner order with proprietary bead/PCR arms dotted."),
      '<h2>Read layout</h2>',
      panel([strand_row(r1,prefix="Read 1  ",suffix=""),strand_row(r2,prefix="Read 2  ",suffix="")],cls="small",caption="The platform paper used 101 × 2 sequencing. Exact primer bases are proprietary."),
      table(("Read","Reports"),[("Read 1","cell label, UMI, then poly(T)"),("Index 1","library index"),("Read 2","cDNA from the random-primed end")]),'</div>'])

if __name__=="__main__":
    OUT.write_text(render(),encoding="utf-8");print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")
