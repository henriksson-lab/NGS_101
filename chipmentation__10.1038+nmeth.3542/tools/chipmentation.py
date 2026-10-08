"""Schmidl et al. ChIPmentation molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows, tagmented_library, tagmentation_scene

TITLE = "ChIPmentation — tagmentation of bead-bound ChIP DNA"
NOTES = "01_chipmentation.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3542">Schmidl et al., <i>Nature Methods</i> (2015)</a>.'
SUMMARY = "Immunoprecipitate crosslinked chromatin first, then let Tn5 fragment and adapter-tag the bead-bound ChIP DNA in a single reaction before indexed PCR."
CAVEAT = "The defining chemistry uses Nextera reagents; the final molecule is assembled from the repository’s canonical Nextera components."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("ChIP-selected DNA")
FINAL_CAPTION = "Indexed Nextera library whose insert was selected by chromatin immunoprecipitation before tagmentation."
SEQUENCING_INTRO = "Paired reads enter the immunoprecipitated fragment at the two Tn5 insertion sites; sample indexes are installed during PCR."

def sections():
    return [
        ("Fragment and immunoprecipitate chromatin", antibody_enzyme_rows("bead capture", "retain target-associated chromatin"),
         "Crosslinked, sonicated chromatin is enriched with a target-specific antibody."),
        ("Tagment on beads", tagmentation_scene("ChIP-selected DNA").rows(),
         "Tn5 replaces conventional end repair, dA addition and adapter ligation; ** marks both adapter-transfer boundaries."),
        ("Reverse crosslinks and amplify", [Row(chunks=[("bead-bound tagged DNA → proteinase/reverse crosslink → indexed PCR", None, False)])],
         "PCR adds complete flow-cell adapters and sample indexes to productive two-ended fragments."),
    ]
