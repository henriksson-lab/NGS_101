"""Mumbach et al. 2016 HiChIP molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import proximity_rows, tagmented_library, tagmentation_scene

TITLE = "HiChIP — protein-directed chromatin contacts"
NOTES = "01_hichip.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3999">Mumbach et al., <i>Nature Methods</i> (2016)</a>.'
SUMMARY = "Build biotin-marked in situ Hi-C contacts, immunoprecipitate contacts associated with a chosen protein, and tagment the bead-bound DNA into an Illumina library."
CAVEAT = "The contact panel derives the original MboI fill and GATCGATC ligation scar. The opposite-strand biotin is present in the duplex but only the top-strand biotin is separately coloured; the final adapter structure is canonical Nextera."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("protein-enriched contact fragment")
FINAL_CAPTION = "Tn5-built paired-end library from immunoprecipitated, biotin-bearing proximity-ligation products."
SEQUENCING_INTRO = "Paired reads recover the two genomic ends of each protein-enriched contact molecule; an internal ligation scar identifies informative contacts."

def sections():
    return [
        ("Create in situ proximity ligations", proximity_rows(),
         "MboI ends are filled with biotin-dATP and ligated in intact nuclei; ** marks the exact boundary in the derived GATCGATC contact scar."),
        ("Immunoprecipitate protein-associated contacts", [Row(chunks=[("contact chromatin + target antibody → protein-directed bead capture", "cbc", False)])],
         "Chromatin is sheared after proximity ligation and enriched with an antibody to the factor or histone mark of interest."),
        ("Tagment bead-bound DNA", tagmentation_scene("protein-enriched contact fragment").rows(),
         "Tn5 installs sequencing ends directly on immunoprecipitated chromatin; ** marks the two adapter-transfer boundaries."),
    ]
