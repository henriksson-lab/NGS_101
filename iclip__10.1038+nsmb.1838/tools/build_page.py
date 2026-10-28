"""Maintainer-only partial molecular inventory; no finished-protocol marker."""
from pathlib import Path
from iclip import l3_scene, linear_cdna, BOUNDARY, protocol_graph
from circular import circularize_ssdna
from chemdraw import panel
from page import head, caveat

OUT = Path(__file__).resolve().parents[1] / 'iclip_partial.html'

def build():
    # Gate is exercised, not waived: no complete marker and no sequencing claim are emitted.
    assert not protocol_graph().publishable
    text = (head('iCLIP — partial source pass') + '<div class="wrap"><h1>iCLIP: partial source pass</h1>'
            + caveat(BOUNDARY)
            + '<p>Two molecular inventory states; the complete reaction workflow is unfinished.</p>'
            + panel(l3_scene().rows(), caption='Source-written L3 with terminal modifications')
            + panel(circularize_ssdna(linear_cdna(), five_prime_phosphate=True).rows(),
                    caption='Derived replicate-1 cDNA circle; variable RNA copy is a placeholder')
            + '</div>')
    target = OUT
    target.write_text(text)
    print(target)

if __name__ == '__main__':
    build()
