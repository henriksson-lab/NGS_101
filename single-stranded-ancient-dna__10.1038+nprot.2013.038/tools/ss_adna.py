"""Gansauge & Meyer single-stranded ancient-DNA library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene

TITLE = "Single-stranded library preparation for ancient or damaged DNA"
NOTES = "01_single-stranded-library.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nprot.2013.038">Gansauge & Meyer, <i>Nature Protocols</i> (2013)</a>.'
SUMMARY = "Heat-denatured individual DNA strands receive a biotinylated 3′ adapter, are captured on streptavidin beads, copied into duplexes, blunt-ligated to a second adapter and completed by indexed PCR."
CAVEAT = "Modified oligos are shown with their published modifications. The damaged insert remains a placeholder because every recovered strand can have a different length and lesion pattern."

CL72 = "ACACTCTTTCCCTACACGACGCTCTTCC"
GES_INDEX2 = "GGAAGAGCGTCGTGTAGGGAAAGAGTGT"
FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("i5", "J" * 8, "cbc", placeholder=True),
    seg("CL72 custom Read 1 site", CL72, "r1"),
    seg("ancient DNA insert", "X" * 34, placeholder=True), seg("dA junction", "A"),
    seg("CL9 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="single-stranded ancient-DNA library")
SEQ_PRIMERS = (sp.custom("Read 1", "CL72", CL72, "Gansauge & Meyer 2013"),
               sp.TRUSEQ["I1"],
               sp.custom("Index 2 (i5)", "Gesaffelstein index 2", GES_INDEX2,
                         "Gansauge single-strand library sequencing setup"),
               sp.TRUSEQ["R2"])
problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                     required_roles=tuple(p.role for p in SEQ_PRIMERS))
if problems: raise ValueError("ancient-DNA library: " + "; ".join(problems))
FINAL_CAPTION = "Indexed PCR supplies full flow-cell arms around the recovered ancient strand. The P5-side primer site is deliberately truncated and uses the published CL72 custom Read 1 primer."
SEQUENCING_INTRO = "CL72 replaces the standard Read 1 primer because the P5-side library arm is truncated. The read begins directly in the captured ancient DNA strand."

def sections():
    ss = Scene(); ss.strand("damaged strand", [seg("damaged DNA", "X" * 30, placeholder=True)], label="heat-denatured input")
    return [
        ("Denature and ligate the first adapter to individual strands", ss.rows() +
         [Row(chunks=[("3′-OH ** 5′-p-CL78-[C3 spacer]10-[TEG-biotin]", "cbc", False)])],
         "CircLigase joins the 3′ end of each input strand to phosphorylated CL78. ** marks the single-strand ligation; the terminal biotin enables capture."),
        ("Immobilize and copy the captured strand",
         [Row(chunks=[("streptavidin bead—biotin—CL78—ancient DNA", "cbc", False)]),
          Row(chunks=[("                         ||||||||||||||  ← polymerase copy", None, False)])],
         "The bead retains the adapter-bearing template through washes. An extension primer copies the complete captured strand into a duplex."),
        ("Blunt-ligate the second adapter and index by PCR",
         [Row(chunks=[("CL53/CL73 duplex ** copied ancient-DNA duplex ** first-adapter end", "me", False)]),
          Row(chunks=[("                    indexed P5/P7 primers complete the library", None, False)])],
         "A second adapter is added by blunt-end ligation. PCR releases and completes the sequenceable molecules; ** marks ligation boundaries."),
    ]
