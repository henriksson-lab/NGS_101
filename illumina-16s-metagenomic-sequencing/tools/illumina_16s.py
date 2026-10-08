"""Illumina two-PCR V3-V4 16S library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg
from chemdraw import Row
from targeted_ngs import two_stage_amplicon

TITLE="Illumina 16S metagenomic sequencing"
NOTES="01_illumina-16s.html"
SOURCE='Illumina <a href="https://assets.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf">16S Metagenomic Sequencing Library Preparation, 15044223 B</a>.'
SUMMARY="A V3–V4 locus PCR installs Nextera overhangs; a limited second PCR adds dual sample indexes and complete Illumina flow-cell arms."
CAVEAT="The locus contains the guide's degenerate IUPAC positions; the schematic uses an explicit target placeholder rather than pretending one bacterial sequence represents the community."
FINAL_LIBRARY,SEQ_PRIMERS=two_stage_amplicon("16S V3–V4 amplicon",40)
FINAL_CAPTION="Dual-indexed V3–V4 amplicon after the second PCR."
SEQUENCING_INTRO="All four landing sites come from the Nextera overhang/index-PCR architecture."
FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"

def sections():
    return [
      ("PCR1 — select V3–V4 and install universal tails",[
       Row(chunks=[("5'- Nextera overhang — CCTACGGGNGGCWGCAG -3'", "s5", False)]),
       Row(chunks=[("              bacterial 16S V3–V4 target", None, False)]),
       Row(chunks=[("3'- Nextera overhang — GACTACHVGGGTATCTAATCC -5'", "s7", False)]),
      ],"The locus-specific 3′ ends determine the amplicon; their 5′ overhangs become PCR2 priming sites."),
      ("PCR2 — add dual indexes and flow-cell arms",[
       Row(chunks=[("P5—i5—S5/Read1 — V3–V4 insert — Read2/S7—i7—P7", "cbc", False)]),
      ],"Limited-cycle index PCR turns the cleaned PCR1 product into a sequencing-ready library."),
    ]
