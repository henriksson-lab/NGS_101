#!/usr/bin/env python3
"""Original primer transcription and PCR1-to-readout preservation regressions."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import illumina_16s as M
import seqprimers as sp
from checks import Check, Source, have
c=Check()
root=data_root()/"sources"/H.parent.name
source=Source(root/"illumina_16s-metagenomic-library-prep-guide-15044223-b.pdf.txt",
              "Defining Illumina15044223B original; tools/get_sources.py --slugs "+H.parent.name)
c.section("complete source PCR1 primer strings")
if have(c,source,label="Illumina15044223B original"):
    text=source.path.read_text()
    start=text.index("16S Amplicon PCR Forward Primer =")
    end=text.index("This method can also",start)
    compact=re.sub(r"\s+","",text[start:end])
    c("full modified/degenerate forward ordering sequence",M.FWD in compact)
    c("full modified/degenerate reverse ordering sequence",M.REV in compact)
    # Source-positive endpoints must survive the second PCR rather than being replaced
    # by a generic insert shell; reverse orientation preserves the full IUPAC alphabet.
    c("final library retains the entire source forward PCR1 primer",M.FWD in M.FINAL_LIBRARY.top())
    c("final library retains the entire source reverse PCR1 primer",
      M.iupac_complement(M.REV)[::-1] in M.FINAL_LIBRARY.top())
c.section("read entry and variable-length biological interval")
for primer,locus in ((M.SEQ_PRIMERS[0],M.FWD_LOCUS),(M.SEQ_PRIMERS[-1],M.REV_LOCUS)):
    hit=sp.locate(M.FINAL_LIBRARY,primer,n=len(locus))
    c(primer.role+" starts with the complete source locus-primer region",hit.reads,locus)
    preview=sp.locate(M.FINAL_LIBRARY,primer,n=300).reads
    c(primer.role+" stops its preview at the unassigned biological interval",preview,locus+"…")
contents={x.read:x.content for x in sp.read_contents(M.FINAL_LIBRARY,M.SEQ_PRIMERS,M.READ_LENGTHS,
                                                   acquired_roles=M.ACQUIRED_ROLES)}
c("Read1 does not acquire the reverse-primer region at a fabricated fixed coordinate",
  "reverse V4" in contents["Read 1"],False)
c("Read2 does not acquire the forward-primer region at a fabricated fixed coordinate",
  "forward V3" in contents["Read 2"],False)
c.report()
