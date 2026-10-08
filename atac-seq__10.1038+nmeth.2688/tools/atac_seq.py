"""Buenrostro et al. 2013 bulk ATAC-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Row
from chromatin_epigenetics import tagmented_library, tagmentation_scene

TITLE = "ATAC-seq — transposition into accessible chromatin"
NOTES = "01_atac-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.2688">Buenrostro et al., <i>Nature Methods</i> (2013)</a>.'
SUMMARY = "Tn5 simultaneously cleaves accessible native chromatin and transfers sequencing-adapter ends; extension and indexed PCR complete the paired-end library."
CAVEAT = "The paper uses Nextera transposase and prints its PCR primers. The final construct is assembled from the repository’s canonical Nextera components."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("accessible chromatin fragment")
FINAL_CAPTION = "Productive heterologous Tn5 product after 72 °C extension and limited indexed PCR."
SEQUENCING_INTRO = "Paired reads enter the accessible chromatin fragment from its two transposition sites; the index reads identify the PCR-added sample indexes."

def sections():
    return [
        ("Expose native chromatin", [Row(chunks=[("nucleosome -- linker -- nucleosome     [open regulatory DNA]", None, False)])],
         "Mild lysis releases nuclei while retaining chromatin architecture."),
        ("Tagment accessible DNA", tagmentation_scene("accessible chromatin fragment").rows(),
         "Tn5 cleavage and mosaic-end transfer are coupled; ** marks both transferred adapter–insert boundaries."),
        ("Extend and amplify", [Row(chunks=[("72 °C fill/extension → indexed PCR → size-distributed library", None, False)])],
         "The initial extension makes both transposase-tagged ends amplifiable; PCR cycle number is limited by a qPCR side reaction."),
    ]
