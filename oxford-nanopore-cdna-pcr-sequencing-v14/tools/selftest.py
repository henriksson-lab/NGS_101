#!/usr/bin/env python3
from pathlib import Path
import re,sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
import cdna_pcr as M
c=Check()
s=Source(sources_root()/HERE.parent.name/"chemistry-technical-document.html.txt","ONT CHTD_500 v1 revAV; ref/MANIFEST.md")
if have(c,s):
 text="".join(s.path.read_text().split())
 c("source identity","CHTD_500_v1_revAV_24July2026" in text,True)
 c("source transcription CRTA", M.CRTA in text, True)
 c("source transcription RTP", M.RTP in text, True)
 c("source transcription SSPII_AS_WRITTEN", M.SSPII_AS_WRITTEN in text, True)
 c("source transcription CDNA_PRIMER_FORWARD", M.CDNA_PRIMER_FORWARD in text, True)
 c("source transcription CDNA_PRIMER_REVERSE", M.CDNA_PRIMER_REVERSE in text, True)
 c("source transcription RAPID_ADAPTER", M.RAPID_ADAPTER in text, True)
c.report()
