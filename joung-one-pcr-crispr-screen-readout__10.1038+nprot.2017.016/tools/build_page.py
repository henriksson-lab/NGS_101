#!/usr/bin/env python3
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import pooled_crispr_page as P

OUT = HERE.parent / "joung-one-pcr-crispr-screen-readout.html"


def render() -> str:
    body = P.one_pcr_readout("Genomic DNA to sequencing library — one PCR")
    return P.page("Joung one-PCR CRISPR screen readout", "01_joung-one-pcr.html",
        'Source: <a href="https://doi.org/10.1038/nprot.2017.016">Joung et al., <i>Nature Protocols</i> (2017)</a>. This is one PCR, not the PCR2 half of a nested workflow.', body)


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
