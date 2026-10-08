#!/usr/bin/env python3
"""Build the lentiCRISPR v1 screening page."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P

OUT = HERE.parent / "lenticrispr-v1-screening.html"


def render() -> str:
    body = f'''<h2>Transfer vector</h2>
{P.map_panel("v1")}
<p>The guide cassette is upstream of cPPT/CTS in v1. That order matters when a reverse readout primer is placed in cPPT.</p>
{P.cloning()}
{P.one_pcr_readout("Compatible one-PCR genomic readout")}
<p>The original Shalem screen used a nested two-PCR readout; its details remain in the research notes and are not blended into the compatible Joung one-PCR example above.</p>'''
    return P.page("lentiCRISPR v1 screening", "01_lenticrispr-v1.html",
        'Defining source: <a href="https://doi.org/10.1126/science.1247005">Shalem et al., <i>Science</i> (2014)</a>. The final-library example uses the separately cited Joung one-PCR readout.', body)


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
