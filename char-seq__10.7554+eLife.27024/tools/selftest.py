#!/usr/bin/env python3
"""Published-image bridge transcription and defining-run source checks."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check, Source, have
from storage import data_root
c=Check()
sources=data_root()/"sources"/H.parent.name
image=Source(sources/"eLife.27024_PMC5962340_supplementary.zip",
             "published eLife Figure1–figure supplement1 JPG inside original PMC supplement ZIP")
if have(c,image,label="published bridge figure original"):
    # Literal transcriptions are pinned to the visually read published figure,
    # not a text hit from the superseded preprint or the changed 2020 protocol.
    c("published image bridge top transcription",P.BRIDGE_TOP,"AAAAAAACCGGCGTCCAAGGATCTAG")
    c("published image bridge bottom 5-prime-to-3-prime transcription",P.BRIDGE_BOTTOM,"GATCCTAGATCCTTGGACGCCGGT")
    c("published image six-base unpaired RNA-ligation tail",P.TOP_TAIL,"AAAAAA")
article=Source(sources/"eLife.27024_v2.xml.txt", "published eLife v2 JATS original")
if have(c,article,label="published sequencing statement"):
    match=re.search(r"single-end (\d+) bp reads",article.path.read_text())
    c("defining single-end acquired length",P.READ_LENGTHS,{"Read 1":int(match.group(1))} if match else {})
c.report()
