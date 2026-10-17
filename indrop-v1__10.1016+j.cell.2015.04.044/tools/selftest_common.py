"""Focused source/model checks for inDrop v1 and v2."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from checks import Check, Source, have
from storage import sources_root
import indrop as I

CELL_SI = Source(sources_root() / "indrop-v1__10.1016+j.cell.2015.04.044"
                 / "j.cell.2015.04.044_1-s2.0-S0092867415005000-mmc1.pdf.txt",
                 "Cell 2015 Supplemental Information mmc1 (ars.els-cdn.com, see MANIFEST.tsv); "
                 "python3 tools/doctext.py makes the .txt")


def _table_s5_pcr1(text: str) -> str:
    """PCR enrichment primer 1 rejoined across the wrapped 'primers:' row label."""
    lines = text.splitlines()
    i = next(n for n, line in enumerate(lines) if line.startswith("PCR enrichment "))
    if lines[i + 1].strip() != "primers:":
        raise ValueError("Table S5 PCR primer layout changed")
    return lines[i].split()[-1] + lines[i + 2].strip()


def run(protocol: str) -> None:
    check = Check(); check.section("published sequence assembled from model components")
    if protocol == I.V1:
        assembled = I.LEADER + I.T7_PROMOTER + I.MID + I.PE1[:24]
        check("acrydite primer bases (Cell 2015 Table S5; same in the v2 supplement)", assembled, I.SOURCE_ACRYDITE_BASES)
        if have(check, CELL_SI):
            check("Table S5 PCR enrichment primer 1, rejoined across its line wrap",
                  _table_s5_pcr1(CELL_SI.path.read_text(encoding="utf-8")), I.V1_PCR1)
    elif protocol == I.V2:
        check("PE2-N6", I.PE2 + "N" * 6, I.SOURCE_PE2_N6)
    else: raise ValueError(protocol)
    check.report()

