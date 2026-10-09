#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scnmt_seq as M
import seqprimers as sp
from chemdraw import Scene, annotation_rows, panel
from page import head, info

OUT=H.parent/"scnmt-seq.html"

def main():
    p=[head(M.TITLE),'<div class="wrap">',f'<h1>{M.TITLE}</h1>',
       f'<p class="research-notes"><a href="{M.NOTES}">Research notes</a></p>',
       info(M.SOURCE),info(M.SUMMARY)]
    for h,rows,cap in M.sections(): p += [f'<h2>{h}</h2>',panel(rows,cls="long",caption=cap)]
    for name,lib,primers,caption,intro in (
        ("RNA branch",M.RNA_LIBRARY,M.RNA_PRIMERS,"Nextera XT library made from Smart-seq2 cDNA.","Nextera sequencing reads the RNA-derived insert."),
        ("DNA branch",M.DNA_LIBRARY,M.DNA_PRIMERS,"Historical scBS-seq random-primed library.","The dedicated historical iPCRTag primer reads i7; there is no i5.")):
        p += [f'<h2>{name} final library</h2>',panel([*Scene.duplex(list(lib),label=name).rows(),*annotation_rows(lib)],cls="long",caption=caption),
              sp.section(lib,primers,heading=name+" sequencing primers",intro=intro,required_roles=tuple(x.role for x in primers))]
    p.append('</div>')
    OUT.write_text('\n'.join(p),encoding='utf-8'); print(f'wrote {OUT}')
if __name__ == '__main__': main()
