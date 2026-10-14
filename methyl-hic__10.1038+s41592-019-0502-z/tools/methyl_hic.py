"""Bulk Methyl-HiC molecular model (Li et al. 2019)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import DPNII
from single_cell_hic import ContactWorkflow, SEQ_PRIMERS, inferred_truseq_library

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
LIBRARY_CAPTION = "Bisulfite-converted contact insert in the TruSeq-compatible shell named by the supplementary protocol. The particular indexed adapter member is not reported."
LIBRARY_CAVEAT = ("INFERRED ADAPTER BASES — the supplement names TruSeq PCR primers and "
                   "an indexed Illumina adapter, but does not print the ligated adapter "
                   "oligos or selected index. Canonical TruSeq-compatible arms therefore "
                   "show primer placement without claiming a source transcription.")
REQUIRED_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")
ADAPTER_BEFORE_CONVERSION = True

def contact_product(): return WORKFLOW.junction()
def final_library():
    return inferred_truseq_library(contact_product(), "bisulfite-converted Methyl-HiC library")
