"""Bulk Methyl-HiC molecular model (Li et al. 2019)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import DPNII
from single_cell_hic import ContactWorkflow, unresolved_illumina_library

WORKFLOW = ContactWorkflow("Methyl-HiC", DPNII, "bulk nuclei", "biotin-fill", "bisulfite-PCR", capture=True)
TITLE = "Methyl-HiC — chromatin contacts plus DNA methylation"
NOTES = "01_methyl-hic.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0502-z">Li et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY = "In situ DpnII Hi-C preserves proximity junctions; on-bead adapter ligation followed by bisulfite conversion turns each paired-end contact molecule into a joint contact-and-methylation readout."
STEPS = (
    ("Digest nuclei and fill ends", "DpnII digestion is followed by fill-in with biotinylated nucleotides."),
    ("Proximity ligation", "T4 DNA ligase joins spatially adjacent filled restriction fragments in intact nuclei."),
    ("Shear and enrich", "Reverse crosslinks, sonicate to about 400 bp and capture biotin-marked junctions on streptavidin beads."),
    ("Ligate adapters, then convert", "Ligate sequencing adapters on beads before sodium-bisulfite conversion."),
    ("Uracil-tolerant PCR", "KAPA HiFi HotStart Uracil+ mix amplifies the converted library."),
)
READOUT = "Paired-end reads span the two DpnII-derived contact arms while C-to-T conversion reports cytosine methylation."
JUNCTION_CAPTION = "DpnII GATC ends are filled with biotinylated nucleotide and proximity-ligated; ** marks the contact junction."
LIBRARY_CAPTION = "Bisulfite-converted contact insert between adapter regions whose exact oligonucleotide sequences are not printed in the defining paper."
LIBRARY_CAVEAT = "INFERRED — the paper states adapter ligation but does not identify or print the exact adapter oligonucleotides; no TruSeq generation is assumed."
SEQUENCING_UNAVAILABLE = "The paper specifies paired-end Illumina sequencing on HiSeq 2500/4000, but not a complete adapter or custom sequencing-primer sequence set, so base-level primer placement cannot be verified."

def contact_product(): return WORKFLOW.junction()
def final_library(): return unresolved_illumina_library(contact_product(), "bisulfite-converted Methyl-HiC library", indexed=None)
