#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import sci_car as M
import nextera as nx
import illumina as il
from chemdraw import revcomp
from checks import Check
c=Check(); c("two distinct modality libraries",[x[0] for x in M.FINAL_LIBRARIES],["RNA library","ATAC library"])
# Amini 2014 Supplementary Table 4 (publisher MOESM2, corrected 20 May 2015), the
# architecture the inferred ATAC transposons follow.
c.section("inferred ATAC architecture against Amini 2014 Table 4")
c("P5_i5_1 transposon",nx.S5+M.CONNECTOR_A+"TATAGCCT"+M.EXT_A+nx.ME,"TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG")
c("P7_i7_1 transposon",nx.S7+M.CONNECTOR_B+"CGAGTAAT"+M.EXT_B+nx.ME,"GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG")
c("Read 1 primer",M.ATAC_R1,"GCGATCGAGGACGGCAGATGTGTATAAGAGACAG")
c("Read 2 primer",M.ATAC_R2,"CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG")
c("Index 1 primer",revcomp(M.ATAC_R2),"CTGTCTCTTATACACATCTGAGGCGGAGACGGTG")
# Published Supplementary Materials (NIHMS1033869-supplement-1.pdf.txt lines 88, 168,
# 178, 218; pp. 16-18) and GEO GSE117089 sample extraction protocol (line 344488),
# transcribed primer parts versus the library segments the model uses.
c.section("GEO GSE117089 primer transcription")
c("RNA P5 primer tail is the drawn Read 1 arm",il.TRUSEQ_READ1,"ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("RT primer arm is the 3' end of that arm",il.TRUSEQ_READ1.endswith("ACGACGCTCTTCCGATCT"),True)
c("P5 grafting sequence",il.P5,"AATGATACGGCGACCACCGAGATCTACAC")
c("P7 primer 5' part",revcomp(il.P7_RC),"CAAGCAGAAGACGGCATACGAGAT")
c("P7 primer 3' S7 part",nx.S7,"GTCTCGTGGGCTCGG")
c("ATAC P5 primer 3' S5 part",nx.S5,"TCGTCGGCAGCGTC")
c("RNA cycles",M.READ_LENGTHS["RNA library"],{"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52})
c.report()
