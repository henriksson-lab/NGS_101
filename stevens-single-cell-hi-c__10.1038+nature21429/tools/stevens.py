"""Stevens et al. 2017 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from restriction import MBOI
from single_cell_hic import (ContactWorkflow, HISTORICAL_PE_READ2,
                             historical_inline_pe_library)
WORKFLOW = ContactWorkflow("Stevens single-cell Hi-C", MBOI,
    "before digestion and ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Stevens single-cell Hi-C — single nuclei before AluI ligation"
NOTES = "01_stevens-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21429">Stevens et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "Individual nuclei are isolated first; MboI contact ends are biotin-filled and proximity-ligated. AluI is used later to fragment captured DNA for library preparation."
STEPS = (("Isolate nuclei", "Lyse fixed cells and distribute individual nuclei before library construction."),
         ("Contact digest", "MboI cuts GATC in each isolated nucleus."),
         ("Mark and ligate", "Fill MboI ends with biotin-14-dATP, then ligate spatial neighbours."),
         ("Capture and library prep", "Capture contacts, digest captured DNA with AluI, A-tail, ligate a 3-bp-tagged Illumina adapter and PCR. Nextera XT was a published alternative."))
READOUT = "Paired reads report the two MboI fragments joined within one nucleus; the first three bases identify the library adapter."
JUNCTION_CAPTION = "Biotin-filled MboI ends form the contact boundary. AluI acts only after contact capture. ** marks ligation."
LIBRARY_CAPTION = "The source-printed adapter-ligation route, with the CAA 3-bp identification tag shown. The paper also reports a dual-index Nextera XT alternative."
def contact_product(): return WORKFLOW.junction()
SEQ_PRIMERS = (
    sp.custom("Read 1", "historical Illumina Read 1", il.TRUSEQ_READ1,
              "Stevens et al. 2017 Supplementary Methods"),
    sp.custom("Read 2", "historical paired-end Read 2", HISTORICAL_PE_READ2,
              "Stevens et al. 2017 Supplementary Methods"),
)
REQUIRED_ROLES = ("Read 1", "Read 2")
def final_library(): return historical_inline_pe_library(contact_product(), "Stevens scHi-C library")
