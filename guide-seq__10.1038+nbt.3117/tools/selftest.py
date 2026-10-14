#!/usr/bin/env python3
"""Source-boundary checks for GUIDE-seq.

The sequence derivations are enforced at import time by ``guide_seq._require``, so
nothing here restates a constant. What remains are claims about the published run
that the construction path cannot make true on its own: which arm each read enters,
that the P7 side is deliberately NOT the stock TruSeq arm, and where the 16-cycle
index-2 read puts the two P5-side identifiers.
"""
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import guide_seq as M
import seqprimers as sp
from checks import Check

c = Check()

c.section("the P7 side is an AMP-style arm, not the stock TruSeq arm")
# Supplementary Table 4 gives the P7 arm as the GSP2 tail completed by P701-P708, which
# is why cartridge positions 13 and 14 need custom primers. A model that had silently
# installed illumina.INDEX1_PRIMER would land these stock primers instead.
for key in ("I1", "R2"):
    c(f"stock TruSeq {key} primer has no site on the GUIDE-seq library",
      sp.locate(M.FINAL_LIBRARY, sp.TRUSEQ[key]), None)
c("the declared Index 1 primer is the reverse complement of the declared Read 2 primer",
  M.INDEX1_SEQ_PRIMER, sp.revcomp(M.READ2_SEQ_PRIMER))

c.section("read entry points, 151 | 8 | 16 | 151")
entry = {}
for primer in M.SEQ_PRIMERS:
    hit = sp.locate(M.FINAL_LIBRARY, primer)
    entry[primer.role] = None if hit is None else hit.reads_from
c("Read 1 enters the sheared end and reports genomic DNA",
  entry["Read 1"], "adjacent genomic DNA")
c("Read 2 enters the integrated dsODN", entry["Read 2"], "integrated dsODN end")
c("Index 1 reports the P7 barcode", entry["Index 1 (i7)"], "P7 barcode reverse complement")
c("Index 2 reports the Y-adapter barcode first", entry["Index 2 (i5)"], "8-nt sample barcode")

c.section("identifier cycles derived from the declared run lengths")
spans = {(s.read, s.feature.id): s.cycles
         for s in sp.feature_spans(M.FINAL_LIBRARY, M.SEQ_PRIMERS, M.READ_LENGTHS)}
c("the molecular index is read by index 2, cycles 9-16",
  spans.get(("Index 2 (i5)", "guide_seq_molecular_index")), "9–16")
c("the Y-adapter sample barcode is read by index 2, cycles 1-8",
  spans.get(("Index 2 (i5)", "guide_seq_y_adapter_barcode")), "1–8")
c("the P7 barcode is read by index 1, cycles 1-8",
  spans.get(("Index 1 (i7)", "guide_seq_p7_barcode")), "1–8")
content = {r.read: r.content for r in sp.read_contents(
    M.FINAL_LIBRARY, M.SEQ_PRIMERS, M.READ_LENGTHS, acquired_roles=M.ACQUIRED_ROLES)}
c("Read 1 reports no identifier: the molecular index is behind its primer",
  [label for label in ("UMI", "sample index") if label in content["Read 1"]], [])
c("Read 2 reads 22 bp of dsODN before crossing into genomic DNA",
  content["Read 2"].startswith("integrated dsODN end (cycles 1–22); "
                               "adjacent genomic DNA (from cycle 23"), True)

c.section("one ligated end, one PCR-installed end")
geometry = M.adapter_geometry()
c("the adapter-end geometry is published, not unavailable", geometry.availability, "known")
c("the ligating end is cohesive", geometry.kind, "cohesive")
c("the library carries exactly one dA/dT junction, beside the Read 1 arm",
  [s.name for s in M.FINAL_LIBRARY if "dT" in s.name or "dA" in s.name],
  ["adapter 3' dT : fragment dA"])

c.report()
