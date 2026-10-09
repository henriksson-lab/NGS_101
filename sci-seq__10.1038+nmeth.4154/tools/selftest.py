#!/usr/bin/env python3
import build_page
from checks import Check
M=build_page.M; c=Check(); c.section("source transcription")
c("Amini Read 1 primer", M.R1_SEQ, "GCGATCGAGGACGGCAGATGTGTATAAGAGACAG")
c("Amini Read 2 primer", M.R2_SEQ, "CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG")
c("Amini Index 1 primer", M.I1_SEQ, "CTGTCTCTTATACACATCTGAGGCGGAGACGGTG")
build_page.render(M); c.report()
