"""EMP ITS1f–ITS2 fusion-primer library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, revcomp

TITLE = "Earth Microbiome Project ITS amplicon sequencing"
NOTES = "01_emp-its.html"
SOURCE = '<a href="https://earthmicrobiome.org/protocols-and-standards/its/">Earth Microbiome Project ITS Illumina Amplicon Protocol (EMP.ITSkabir)</a> and its November 2016 oligo workbook.'
SUMMARY = "One fusion-primer PCR amplifies fungal ITS1, installs the flow-cell arms and a twelve-base Golay sample barcode, and creates landing sites for three locus-specific custom sequencing primers."
CAVEAT = "The published Read 2 primer is a degenerate pool (V and R). The landing diagram uses one valid member, V=G and R=A, so its exact binding can be computed; the ordered degenerate sequence is retained below. The fungal insert length and Golay identity vary."

FWD_ARM = "AATGATACGGCGACCACCGAGATCTACAC"
FWD_LINK = "GG"
FWD_LOCUS = "CTTGGTCATTTAGAGGAAGTAA"
REV_ARM = "CAAGCAGAAGACGGCATACGAGAT"
REV_LINK = "CG"
REV_LOCUS = "GCTGCGTTCTTCATCGATGC"

# November 2016 EMP workbook, sheets Read1, Read2 and Index.
READ1_PRIMER = "TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC"
READ2_PRIMER_ORDERED = "CGTTCTTCATCGATGCVAGARCCAAGAGATC"
READ2_PRIMER_REPRESENTATIVE = "CGTTCTTCATCGATGCGAGAACCAAGAGATC"  # V=G, R=A
INDEX1_PRIMER = "TCTCGCATCGATGAAGAACGCAGCCG"

# The workbook primers extend beyond the PCR primer's locus-specific 3' ends.
R1_EXTENSION = READ1_PRIMER[len(FWD_LOCUS) - 1:]
R2_EXTENSION = READ2_PRIMER_REPRESENTATIVE[len(REV_LOCUS) - 4:]

FINAL_LIBRARY = Construct([
    seg("P5-side fusion arm", FWD_ARM, "p5"),
    seg("forward linker", FWD_LINK),
    seg("ITS1f site", FWD_LOCUS, "r1"),
    seg("Read 1 locus extension", R1_EXTENSION),
    seg("fungal ITS1 insert", "N" * 34, placeholder=True),
    seg("Read 2 locus extension reverse complement", revcomp(R2_EXTENSION)),
    seg("ITS2 site reverse complement", revcomp(REV_LOCUS), "r2"),
    seg("reverse linker", revcomp(REV_LINK)),
    seg("twelve-base Golay barcode", "B" * 12, "cbc", placeholder=True),
    seg("P7-side fusion arm reverse complement", revcomp(REV_ARM), "p7"),
], name="EMP ITS fusion-primer library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "EMP ITS Read 1", READ1_PRIMER,
              "EMP November 2016 workbook, Read1 sheet"),
    sp.custom("Index 1 (i7)", "EMP ITS barcode primer", INDEX1_PRIMER,
              "EMP November 2016 workbook, Index sheet"),
    sp.custom("Read 2", "EMP ITS Read 2 (V=G, R=A member)",
              READ2_PRIMER_REPRESENTATIVE,
              "EMP November 2016 workbook, Read2 sheet",
              "One explicit member of the published degenerate primer pool."),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2")):
    raise ValueError("EMP ITS final library: " + "; ".join(problems))

FINAL_CAPTION = "EMP.ITSkabir paired-end, single-index fusion-primer product. The twelve-base Golay barcode is read in its own custom index read."
SEQUENCING_INTRO = "Three published custom primers sequence the fungal insert from both ends and read the twelve-base Golay barcode. There is no i5/Index 2 read. The Read 2 oligo shown for landing is one member of the ordered V/R-degenerate pool."

def sections():
    return [
        ("One PCR installs the complete fusion-primer construct", [
            Row(chunks=[("P5—GG—ITS1f → fungal ITS1 ← ITS2—CG—Golay12—P7", "cbc", False)]),
        ], "Both locus selection and platform-arm installation occur in the same amplification."),
        ("Custom sequencing primers cross into the locus", [
            Row(chunks=[("Read 1: ITS1f suffix + 19 locus bases →", "r1", False)]),
            Row(chunks=[("← Read 2: ITS2 suffix + 15 degenerate locus bases", "r2", False)]),
            Row(chunks=[("Index 1: four locus-extension bases + ITS2/CG → Golay12", "cbc", False)]),
        ], "The longer locus-matching primers raise melting temperature in the low-complexity amplicon. The index primer approaches the barcode through the reverse-primer side."),
    ]
