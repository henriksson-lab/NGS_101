"""Schwalb et al. 2016 TT-seq molecular workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row

TITLE = "TT-seq — transient transcriptome sequencing"
NOTES = "01_tt-seq.html"
SOURCE = ('Defining paper: <a href="https://doi.org/10.1126/science.aad9841">'
          'Schwalb et al., <i>Science</i> (2016)</a>.')
SUMMARY = ("A five-minute 4-thiouridine pulse marks newly made RNA. Fragmenting RNA before "
           "thiol capture makes only the newly synthesized pieces recoverable for sequencing.")
CAVEAT = ("The paper names the NuGEN Ovation Human Blood RNA-seq kit but does not print its "
          "adapter oligos. The Illumina-compatible library below is therefore an inferred "
          "run geometry; the TT-seq label, fragmentation and capture steps are source-defined.")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5", inferred=True),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
    seg("strand-specific cDNA from captured RNA", "X" * 38, placeholder=True),
    seg("dA junction", "A", inferred=True),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc",
        placeholder=True, inferred=True),
    seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="inferred single-index TT-seq library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("TT-seq final library: " + "; ".join(problems))

FINAL_CAPTION = ("INFERRED — Illumina-compatible paired-end, six-base single-index geometry; "
                 "the defining supplement reports paired 50-base reads plus a six-base barcode.")
SEQUENCING_INTRO = ("The inferred standard primer sites show how the named HiSeq run enters the "
                    "captured cDNA. Adapter bases are not asserted as NuGEN kit disclosures.")


def sections():
    return [
        ("Pulse-label nascent RNA", [
            Row(chunks=[("pre-existing RNA     --------------------------", None, False)]),
            Row(chunks=[("new RNA (5 min)      ----4sU---------4sU-------", "umi", False)]),
        ], "Cells receive 500 µM 4sU for five minutes; labelled spike-ins are added during extraction."),
        ("Fragment before selection", [
            Row(chunks=[("long labelled RNA -- BioRuptor 30 s ON / 30 s OFF --> shorter labelled RNA pieces", None, False)]),
        ], "Fragmentation precedes labelled-RNA purification; that ordering is the distinguishing TT-seq step."),
        ("Biotinylate and capture labelled pieces", [
            Row(chunks=[("RNA—4sU  + thiol-reactive biotin  ->  RNA—4sU—biotin", "umi", False)]),
            Row(chunks=[("RNA—4sU—biotin  -> streptavidin capture -> DTT elution", None, False)]),
        ], "Only fragments carrying a pulse-labelled thiol are enriched; unlabeled RNA is washed away."),
        ("INFERRED — convert captured RNA to the sequencing library", [
            Row(chunks=[("DNase treatment -> strand-specific NuGEN Ovation library preparation -> paired-end library", None, True)]),
        ], "The defining supplement delegates this construct-changing step to the commercial kit."),
    ]
