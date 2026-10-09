"""Illumina two-PCR V3-V4 16S library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg
from chemdraw import Construct, Scene, duplex_rows, revcomp
from targeted_ngs import two_stage_amplicon

TITLE="Illumina 16S metagenomic sequencing"
NOTES="01_illumina-16s.html"
SOURCE='Illumina <a href="https://assets.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf">16S Metagenomic Sequencing Library Preparation, 15044223 B</a>.'
SUMMARY="A V3–V4 locus PCR installs Nextera overhangs; a limited second PCR adds dual sample indexes and complete Illumina flow-cell arms."
CAVEAT="The locus contains the guide's degenerate IUPAC positions; the schematic uses an explicit target placeholder rather than pretending one bacterial sequence represents the community."
FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
FWD_LOCUS="CCTACGGGNGGCWGCAG"
REV_LOCUS="GACTACHVGGGTATCTAATCC"
FWD_OVERHANG=FWD[:-len(FWD_LOCUS)]
REV_OVERHANG=REV[:-len(REV_LOCUS)]

PCR1_PRODUCT=Construct([
 seg("forward Nextera overhang",FWD_OVERHANG,"s5"),
 seg("forward V3 locus primer",FWD_LOCUS,"r1",placeholder=True),
 seg("16S V3–V4 target","N"*40,placeholder=True),
 seg("reverse V4 locus primer complement",revcomp(REV_LOCUS),"r2",placeholder=True),
 seg("reverse Nextera overhang complement",revcomp(REV_OVERHANG),"s7"),
],name="tailed V3–V4 PCR1 product")

FINAL_LIBRARY,SEQ_PRIMERS=two_stage_amplicon("16S V3–V4 amplicon",40)
FINAL_CAPTION="Dual-indexed V3–V4 amplicon after the second PCR."
SEQUENCING_INTRO="All four landing sites come from the Nextera overhang/index-PCR architecture."


def pcr1_scene():
    scene=Scene.duplex(list(PCR1_PRODUCT),label="PCR1 product")
    scene.labels("top")
    return scene

def sections():
    return [
      ("PCR1 — select V3–V4 and install universal tails",pcr1_scene().rows(),
       "The locus-specific 3′ ends determine the amplicon; their 5′ overhangs become PCR2 priming sites."),
      ("PCR2 — add dual indexes and flow-cell arms",
       duplex_rows(FINAL_LIBRARY,label="PCR2 product"),
       "Limited-cycle index PCR turns the cleaned PCR1 product into a sequencing-ready library."),
    ]
