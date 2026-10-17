#!/usr/bin/env python3
"""Source transcription of the wellDA-seq oligo families (skips when sources are absent)."""
from pathlib import Path
import re
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
import illumina as il
import nextera as nx
from checks import Check, Source, have
from storage import data_root

check = Check()
chem.protocol_graph()
root = data_root() / "sources" / HERE.parent.name
article = Source(root / "2024.09.06.611519_v1.full.pdf.txt", "bioRxiv v1 full text; tools/get_sources.py")
supp = Source(root / "2024.09.06.611519_media-1.pdf.txt", "bioRxiv v1 supplement media-1; tools/get_sources.py")
if have(check, article, label="Methods oligo transcription"):
    text = re.sub(r"\s+", "", article.path.read_text())
    for name, seq in (("DEpi_ME_t1", chem.DEPI_ME_T1), ("DEpi_ME_t2", chem.DEPI_ME_T2),
                      ("ME reverse", "/phos/" + chem.ME_REVERSE), ("Idx7", chem.ATAC_IDX7),
                      ("Idx5", chem.ATAC_IDX5)):
        check(name + " in Methods", seq in text, True)
if have(check, supp, label="Supplementary Tables 5-6 primer families"):
    lines = supp.path.read_text().splitlines()
    def family(prefix, suffix):
        pat = re.compile("^" + prefix + "[ACGT]{8}" + suffix + "$")
        return sum(bool(pat.match(l.strip())) for l in lines)
    check("ATAC-CB1 S5XX primers", family(il.P5 + chem.ATAC_S5_SPACER, chem.ATAC_HANDLE1), 72)
    check("ATAC-CB2 N7XX primers", family(il.P7, chem.ATAC_HANDLE2), 72)
    check("DNA-CB1 S5XX primers", family(il.P5, nx.S5), 72)
    check("DNA-CB2 IN_N7XX primers", family(chem.DNA_IN_LINK, nx.S7), 72)
    check("DNA-CB3 OUT_N7XX primers", family(il.P7, chem.DNA_IN_LINK), 11)
check.report()
