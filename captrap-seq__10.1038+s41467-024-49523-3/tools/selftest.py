#!/usr/bin/env python3
"""Source-transcription checks for CapTrap-seq (published Methods, PMC11211341)."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
import storage
from checks import Check, Source, have

c=Check()
jats=Source(storage.sources_root()/"captrap-seq__10.1038+s41467-024-49523-3"/"PMC11211341_PMC11211341.xml.txt",
            "PMC Cloud https://pmc-oa-opendata.s3.amazonaws.com/PMC11211341.1/PMC11211341.1.xml (doctext twin)")
if have(c,jats,label="published Methods oligo transcription"):
    text=jats.path.read_text(encoding="utf-8")
    for name,seq in [("OligodT 16VN",P.OLIGO_DT_16VN),("5' linker upper",P.LINKER5_UPPER+"-P"),
                     ("5' linker lower",f"P-{P.LINKER5_LOWER_PRINTED}-P"),
                     ("3' linker upper",P.LINKER3_UPPER+"-P"),("3' linker lower",P.LINKER3_LOWER_PRINTED),
                     ("second primer UMI8",f"5′ {P.UMI8_PRIMER} 3′"),
                     ("forward LA-PCR",f"forward 5′ {P.LA_FORWARD} 3′"),("reverse LA-PCR",f"reverse 5′ {P.LA_REVERSE} 3′")]:
        c(f"{name} printed verbatim in published Methods",seq in text)
c.report()
