#!/usr/bin/env python3
from pathlib import Path
import re,sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import native_barcoding as M
c=Check()
s=Source(sources_root()/HERE.parent.name/"chemistry-technical-document.html.txt","ONT CHTD_500 v1 revAV; ref/MANIFEST.md")
if have(c,s):
 text="".join(s.path.read_text().split())
 c("source identity","CHTD_500_v1_revAV_24July2026" in text,True)
 c("source transcription LA_TOP", M.LA_TOP in text, True)
 c("source transcription NATIVE_ADAPTER_BOTTOM", M.NATIVE_ADAPTER_BOTTOM in text, True)
 for name,seq in M.NATIVE_BARCODES.items():
  c("source barcode row "+name, name+seq+M.revcomp(seq) in text, True)
c.report()
