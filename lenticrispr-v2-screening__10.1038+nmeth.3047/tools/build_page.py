#!/usr/bin/env python3
"""Build the lentiCRISPR v2 screening page."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P

OUT = HERE.parent / "lenticrispr-v2-screening.html"


def render() -> str:
    body = f'''<h2>Transfer vector</h2>
{P.map_panel("v2")}
<p>v2 moves cPPT/CTS upstream of the hU6 guide cassette. A cPPT reverse primer that gives a short v1 product therefore must not be assumed to work on v2.</p>
{P.cloning()}
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The Moffat LCV2::TKOv3 two-step indexed readout is shown on its own page; it is not an extra step in the Joung one-PCR protocol above.</p>'''
    return P.page("lentiCRISPR v2 screening", "01_lenticrispr-v2.html",
        'Defining source: <a href="https://doi.org/10.1038/nmeth.3047">Sanjana, Shalem &amp; Zhang, <i>Nature Methods</i> (2014)</a>. The final-library example uses the separately cited Joung one-PCR readout.', body)


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
