#!/usr/bin/env python3
"""Build the diagram-centric HyDrop-ATAC chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import hydrop_atac as H
from chemdraw import Scene, annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "hydrop-atac.html"


def oligos() -> str:
    bead = H.bead_oligo()
    rows = [
        oligo("Final split-pool-barcoded bead oligo", list(bead),
              mods="/5Acryd//iThioMC6-D/"),
        oligo("HYi7 bulk-PCR primer", H.hyi7_primer(), mods="last two bonds are phosphorothioate"),
        oligo("HYi5 bulk-PCR primer", H.hyi5_primer()),
    ]
    return ("<h2>Bead oligo and PCR primers</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + panel([strand_row(bead), *annotation_rows(bead)], cls="small",
                    caption="Three 96-way split-pool extensions add barcode 1, barcode 2 and barcode 3. The completed 117-nt oligo ends in the s7 capture sequence."))


def construction() -> str:
    template = H.gap_filled_template()
    product = H.bead_primed_product()
    lib = H.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Tagment nuclei in bulk with s5 and s7 Tn5 adaptors</h3>
{panel(H.tagmented_scene().rows(), cls="small",
       caption="The amplifiable s5/s7 fragment is shown before the two 9-nt gaps are filled.")}

<h3>(2) Co-encapsulate with one bead; DTT releases its primers</h3>
{panel([strand_row(template), *annotation_rows(template)], cls="small",
       caption="The initial 72 °C step fills the gaps, creating s7' at the end captured by the bead primer.")}
{panel(H.capture_scene().rows(), cls="small",
       caption="The bead oligo's 3' s7 pairs with s7'. Thirteen droplet cycles linearly copy fragments; the cell barcode remains attached to every new strand.")}
{panel([strand_row(product), *annotation_rows(product)], cls="small",
       caption="One cell-barcoded linear-extension product.")}

<h3>(3) HYi5/HYi7 bulk PCR adds flow-cell adapters and sample indices</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small",
       caption="Final HyDrop-ATAC library, P5 to P7'. HYi7 primes within the bead backbone, so T8 and most of the T7 promoter are absent.")}
'''


def sequencing() -> str:
    lib = H.final_library()
    landing_rows = []
    for primer, hit in H.primer_landings(lib):
        landing_rows.append((html.escape(primer.role), html.escape(primer.name),
                             html.escape(", ".join(hit.covers)),
                             f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 50 cycles; cell-barcode Index 1: 52; sample Index 2: 10; Read 2: 50. All four use standard Nextera sequencing primers.')}
{table(("Read", "Primer", "Primer site", "First bases"), landing_rows)}
<h3>Index 1 cell-barcode read</h3>
{table(("Cycles", "Content"), H.index1_layout(), scroll=False)}
'''


def render() -> str:
    return "\n".join([
        head("HyDrop-ATAC library chemistry"), '<div class="wrap">',
        '<h1>HyDrop-ATAC &mdash; droplet single-cell chromatin accessibility</h1>',
        info('Defining source: <a href="https://doi.org/10.7554/eLife.73971">De Rop et al., <i>eLife</i> (2022)</a>.'),
        oligos(), construction(), sequencing(), '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
