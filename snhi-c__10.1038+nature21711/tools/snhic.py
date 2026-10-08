"""Flyamer et al. 2017 single-nucleus Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import DPNII
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
WORKFLOW = ContactWorkflow("snHi-C", DPNII, "before fixation", "none", "MDA")
TITLE = "snHi-C — DpnII contact ligation followed by MDA"
NOTES = "01_snhi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature21711">Flyamer et al., <i>Nature</i> (2017)</a>.'
SUMMARY = "Nuclei are isolated before fixation; DpnII contacts are ligated without biotin enrichment and whole molecules are amplified by phi29 MDA."
STEPS = (("Isolate and fix nuclei", "Prepare intact single nuclei before the restriction reaction."),
         ("Digest and ligate", "DpnII cuts GATC; compatible chromatin ends are proximity-ligated without biotin fill-in."),
         ("Reverse crosslinks", "Purify the complete single-nucleus ligation mixture rather than selecting streptavidin-bound junctions."),
         ("MDA then library prep", "Amplify with phi29 multiple-displacement amplification and convert amplified DNA to an indexed sequencing library."))
READOUT = "Paired reads identify restriction fragments; computational filtering distinguishes contact ligations from other MDA products."
JUNCTION_CAPTION = "A DpnII-compatible contact product. Unlike biotin Hi-C, this protocol does not select a filled, marked junction. ** marks ligation."
LIBRARY_CAPTION = "The contact product passes through phi29 MDA before paired-end library conversion. The supplementary methods do not identify the library kit or adapter sequences."
LIBRARY_CAVEAT = "INFERRED — only the presence of an unspecified Illumina library is established; adapter sequence and index architecture are not reported."
SEQUENCING_UNAVAILABLE = "The source specifies 2 × 125-bp HiSeq 2500 sequencing but not the library kit or sequencing-primer sequences, so primer binding cannot be placed base by base."
def contact_product(): return WORKFLOW.junction()
def final_library(): return unresolved_illumina_library(contact_product(), "snHi-C MDA-derived library", indexed=None)
