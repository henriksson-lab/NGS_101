"""Illumina DNA Prep shotgun-metagenomics model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from targeted_ngs import shotgun_library

TITLE="Illumina DNA Prep shotgun metagenomics"
NOTES="01_shotgun-metagenomics.html"
SOURCE='Illumina <a href="https://support.illumina.com/downloads/illumina-dna-prep-reference-guide-1000000025416.html">DNA Prep Reference Guide 1000000025416 v12</a> and <a href="https://www.illumina.com/content/dam/illumina/gcs/assembled-assets/marketing-literature/novaseq-shotgun-metagenomics-app-note-m-na-00010/novaseq-shotgun-metagenomics-app-note-m-na-00010.pdf">shotgun-metagenomics application note M-NA-00010 v2.0</a>.'
SUMMARY="Total community DNA is tagmented by bead-linked transposomes and index-PCR-completed into a dual-index Nextera library without locus selection."
CAVEAT="Bead-linked transposome reagent composition is proprietary; its public molecular consequence—Nextera adapter ends on fragmented DNA—is modeled exactly."
FINAL_LIBRARY,SEQ_PRIMERS=shotgun_library()
FINAL_CAPTION="A representative community-DNA fragment after bead-linked tagmentation and dual-index PCR."
SEQUENCING_INTRO="Standard Nextera primers read a random community fragment and both sample indexes."

def sections():
 return [
  ("Extract total community DNA",[
   Row(chunks=[("bacteria + archaea + fungi + viruses → pooled genomic DNA",None,False)]),
   Row(chunks=[("no marker-gene primer; no locus selection",None,False)]),
  ],"Extraction determines which community molecules reach library preparation."),
  ("Tagment with bead-linked transposomes",[
   Row(chunks=[("BLT bead ●—Tn5 + community dsDNA → bead-bound tagged fragments", "me", False)]),
   Row(chunks=[("wash bead → remove unbound DNA and reaction components",None,False)]),
  ],"Tn5 fragmentation and adapter transfer occur together; bead immobilization permits washing and saturation-based normalization."),
  ("Index PCR and pool",[
   Row(chunks=[("tagged fragment → reduced-cycle PCR → P5—i5 … insert … i7—P7", "cbc", False)]),
  ],"PCR completes the platform arms and labels the sample, not the organism or molecule."),
 ]
