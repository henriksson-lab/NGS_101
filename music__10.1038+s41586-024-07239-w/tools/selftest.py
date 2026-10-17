#!/usr/bin/env python3
"""MUSIC: build the graph and check oligo transcriptions against Supplementary Table 2.

The authors' public sheet is on disk as export.html (an XLSX) with its text twin
export.xlsx.txt; tools/source_oligo_check.py cannot map that twin to a manifest URL,
so the transcriptions are checked here against the sheet text.
"""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import music as M
from checks import Check, Source, have
from storage import sources_root

c=Check()
M.protocol_graph()
SHEET=Source(sources_root()/H.parent.name/"export.xlsx.txt",
             "python3 tools/get_sources.py --slugs "+H.parent.name+" (catalogue/extra_sources.tsv Google-sheet XLSX export)")
c.section("Supplementary Table 2 (public oligo sheet)")
if have(c,SHEET):
    rows={}
    sheet=None
    for line in SHEET.path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## sheet "):
            sheet=line[9:].strip(); continue
        cells=line.split("\t")
        if len(cells)>1 and cells[0]:
            rows.setdefault((sheet,cells[0]),cells[1])
    def bases(text):
        """Drop IDT modification tokens, phosphorothioate stars and RNA r-prefixes."""
        text=re.sub(r"/[^/]*/","",text).replace("*","")
        return re.sub(r"r([ACGU])",lambda m:"T" if m.group(1)=="U" else m.group(1),text)
    first=lambda name:next((v for (s,k),v in rows.items() if s==name and k!="Name"),"")
    c("RNA linker (RNA portion as DNA letters)",bases(rows.get(("Linkers and Adaptors","RNA_Linker_UMI_F"),"")),M.RNA_LINKER)
    c("DNA linker top",bases(rows.get(("Linkers and Adaptors","DNA_Linker_Top_PolA"),"")),M.DNA_TOP)
    c("DNA linker bottom",bases(rows.get(("Linkers and Adaptors","DNA_Linker_Bot_5OH"),"")),M.DNA_BOTTOM)
    c("Universal Adaptor",bases(rows.get(("Linkers and Adaptors","Universal Adaptor"),"")),M.UNIVERSAL)
    c("index adaptor A1",bases(rows.get(("Linkers and Adaptors","Truseq_R2_index_A1"),"")),M.INDEX_A1)
    for sheet,seq in (("BC1_Top",M.BC1_TOP),("BC1_Bottom",M.BC1_BOTTOM),("BC2_TOP",M.BC2_TOP),
                      ("BC2_Bottom",M.BC2_BOTTOM),("BC3_TOP",M.BC3_TOP),("BC3_Bottom",M.BC3_BOTTOM)):
        c(f"{sheet} member 1",bases(first(sheet)),seq)
c.report()
