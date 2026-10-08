"""Stevens et al. 2017 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import ALUI
from single_cell_hic import ContactWorkflow, inferred_truseq_library
WORKFLOW = ContactWorkflow("Stevens single-cell Hi-C", ALUI,
    "before digestion and ligation", "biotin-blunt", "PCR", capture=True)
TITLE = "Stevens single-cell Hi-C — single nuclei before AluI ligation"
NOTES = "01_stevens-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21429">Stevens et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "Individual nuclei are isolated first; AluI creates blunt fragments that are biotin-marked and proximity-ligated."
STEPS = (("Isolate nuclei", "Lyse fixed cells and distribute individual nuclei before library construction."),
         ("Blunt digest", "AluI cuts AG|CT without a cohesive overhang."),
         ("Mark and ligate", "Biotin-mark blunt ends, then ligate spatial neighbours inside each nucleus."),
         ("Capture and PCR", "Reverse crosslinks, purify, enrich marked junctions and make an indexed Illumina library."))
READOUT = "Paired reads report the two AluI fragments joined within one nucleus."
JUNCTION_CAPTION = "AluI supplies a blunt boundary rather than the GATC cohesive-end fill used by classical Hi-C. ** marks ligation."
LIBRARY_CAPTION = "The paper-level schematic resolves the genomic insert; canonical TruSeq arms are shown as inferred kit structure."
def contact_product(): return WORKFLOW.junction()
def final_library(): return inferred_truseq_library(contact_product(), "Stevens scHi-C library")
