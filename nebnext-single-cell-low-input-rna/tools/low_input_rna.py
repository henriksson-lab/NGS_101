"""NEBNext Single Cell/Low Input RNA Library Prep, E6420."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from nebnext import STANDARD_PRIMERS, indexed_library
TITLE="NEBNext Single Cell / Low Input RNA"
NOTES="01_nebnext-single-cell-low-input-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina">NEB E6420 product manual</a>.'
SUMMARY="A proprietary RT primer and template-switching oligo generate full-length amplifiable cDNA; Ultra II FS then fragments and converts that cDNA into a dual-index Illumina library."
RT=[Row(chunks=[("mRNA  5′-[transcript]-poly(A)-3′",None,False)]), Row(chunks=[("                         ||||||  proprietary RT primer →", "r1",True)])]
TS=[Row(chunks=[("TSO (proprietary) → template switch at the mRNA 5′ end", "r2",True)]), Row(chunks=[("known handle — full-length first-strand cDNA — RT-primer handle",None,False)])]
WTA=Construct([Segment("TSO-derived handle", "X"*12, "r2", placeholder=True, inferred=True), Segment("full-length cDNA", "X"*48, placeholder=True), Segment("RT-primer-derived handle", "X"*12, "r1", placeholder=True, inferred=True)], name="amplified cDNA")
STEPS=(("Prime reverse transcription", RT, "The kit names but does not disclose its RT primer sequence; the poly(A)-directed structural role is shown without inventing bases."),
       ("Template-switch at the transcript 5′ end", TS, "Template switching installs a second PCR handle, enabling full-length cDNA amplification."),
       ("Amplify full-length cDNA", WTA, "The cDNA PCR primer recognizes kit-defined terminal handles; their sequences remain proprietary."),
       ("Ultra II FS fragmentation and library conversion", Construct([Segment("FS-fragmented cDNA", "X"*42, placeholder=True)], name="FS cDNA fragment"), "Enzymatic fragmentation, end repair and dA-tailing precede NEBNext adaptor ligation and indexed PCR."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Paired reads cover fragments distributed across the amplified full-length cDNA; the two index reads identify the library, not the cell."
FINAL_CAPTION="Dual-index library after Ultra II FS conversion. Proprietary RT/TSO handles are generally internal to the pre-amplified cDNA and are not asserted at every final fragment end."
def final_library(): return indexed_library(Construct([Segment("cDNA insert", "X"*42, placeholder=True)]), "NEBNext low-input RNA library")
