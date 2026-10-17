#!/usr/bin/env python3
"""Source-transcription checks for MARS-seq (2014)."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import marsseq as M  # noqa: E402
from mars_graph import protocol_graph  # noqa: E402
from page_parts import render_page  # noqa: E402
from checks import Check  # noqa: E402
check = Check()

# Jaitin et al. 2014, Science 343:776, Supplementary Material, Table S7 "Primers used
# during MARS-seq library construction" -- the publisher's own PDF, on disk as
# scg_jaitin-sm.pdf.  These three primers are assembled from illumina constants, so a
# change there would silently alter what this page presents as verbatim.  The scg_lib_structs
# page is not evidence for any of them.  MARS2_RT2 (CTACACGACGCTCTTCCGATCT) belongs to
# the separate MARS-seq2.0 author protocol and is not used by the 2014 graph.
check("Table S7 barcoded RT primer (XXXXXX/NNNN modelled as N)", M.MARS1_RT1,
      "CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCT"
      + "N" * 10 + "T" * 20 + "N")
check("Table S7 Second RT primer", M.MARS1_RT2_PRINTED, "TCTAGCCTTCTCGCAGCACATC")
check("Table S7 P5_Rd1 PCR forward", M.MARS1_P5_PCR,
      "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
check("Table S7 P7_Rd2 PCR reverse", M.MARS1_P7_PCR,
      "CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
# Table S8, row lig_NNNX4_ix1: /5Phos/GACTNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
check("Table S8 ligation-adapter constant", M.LIG_CONSTANT, "AGATCGGAAGAGCGTCGTGTAG")
# Jaitin 2014 full Supplementary Methods, steps 2–4 (visible pp. 4–5): the
# first-strand pool and its Exo I cleanup precede final pooling and second strand.
graph = protocol_graph()
required = ("reverse_transcription", "pool", "remove_free_primers",
            "pool_after_exonuclease", "second_strand", "ivt")
order = tuple(graph.steps)
check("original Methods reaction order", tuple(order.index(s) for s in required),
      tuple(sorted(order.index(s) for s in required)))
check("Exo I fragments have a discard terminal",
      graph.terminals["exo_primer_fragments"].disposition.value, "discard")
check("Exo I fragments appear on the public page",
      "digested unincorporated RT-primer fragments → discard" in
      render_page("MARS-seq", graph=graph), True)
check("unresolved RT2-dependent Read 1 junction remains inferred",
      next(s.inferred for s in M.final_library("MARS-seq")
           if s.name == "TruSeq Read 1, non-overlap"), True)
check.report()
