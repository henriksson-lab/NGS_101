#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import direct_seq as P
from checks import Check
c=Check();c.section("Direct-seq source transcription")
c("8A8G tract",P.CAPTURE,"AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA")
c("human Gln tRNA", P.HUMAN_TRNA_GLN,
  "GGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCT")
c("nested-PCR reverse primer", P.TRNA_READ2,
  "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCTCTGGACTCTGAATCCAGCGAT")
c("guide library preserves one guide identifier",
  [s.feature.role for s in P.GUIDE_LIB if s.feature],
  ["cell_barcode", "umi", "guide_barcode", "sample_index"])
# A matching library and sequencing-primer list can share a wrong handle.
# The independent Table S3d PCR pair must still reject that inconsistent endpoint.
from chemdraw import Construct, Segment
import nextera
wrong_read1 = Construct([
    Segment(s.name, nextera.READ1_PRIMER, s.tag) if s.name == "TruSeq Read 1 arm" else s
    for s in P.GUIDE_LIB], "wrong Nextera endpoint")
c.raises("published guide PCR pair rejects a Nextera Read1 substitution",
         lambda: P.require_guide_primer_sites(wrong_read1))
c.report()
