"""Maintainer-only molecular inventory for the local C18 construct follow-up."""
from pathlib import Path
from irclip import (rt_primer_backbone, linear_cdna_backbone, conditional_first_extension,
                    protocol_graph, BOUNDARY)
from chemdraw import panel
from page import head, caveat

OUT = Path(__file__).resolve().parents[1] / 'irclip_partial.html'

def build():
    assert not protocol_graph().publishable
    primer=rt_primer_backbone();circle=linear_cdna_backbone().close()
    first=conditional_first_extension(circle)
    text=(head('irCLIP — C18 construct follow-up')+'<div class="wrap"><h1>irCLIP C18 backbone</h1>'
          +caveat('INFERRED — cDNA landing and first-extension mechanism. The inspected procedural source is an accepted manuscript; applicability to the final published procedure remains unverified.')
          +'<p>Molecular inventory states. The complete reaction workflow and sequencing output remain unresolved.</p>'
          +panel(primer.rows(),caption='Source-written barcode1 primer: two non-nucleotide C18 links')
          +panel(circle.rows(),caption='INFERRED — cDNA circle, opened only for display; both C18 links retained')
          +panel(first.product.rows(),caption='INFERRED — conditional P6 first extension; copy stops before second C18')
          +'<p>'+BOUNDARY+'</p></div>')
    OUT.write_text(text);print(OUT)

if __name__=='__main__': build()
