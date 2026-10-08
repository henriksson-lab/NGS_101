"""Henikoff et al. low-salt CUTAC molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows, tagmented_library, tagmentation_scene

TITLE = "CUTAC — tethered tagmentation of nearby accessible DNA"
NOTES = "01_cutac.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1186/s13059-020-02266-4">Henikoff et al., <i>Genome Biology</i> (2020)</a>.'
SUMMARY = "An antibody tethers pA/G–Tn5 to a transcription-associated chromatin mark; low-salt activation shifts tagmentation into the adjacent nucleosome-depleted accessible region."
CAVEAT = "The link from the antibody-bound nucleosome to the neighboring accessible site is spatial, not covalent. The final adapter structure is canonical Nextera."
FINAL_LIBRARY, SEQ_PRIMERS = tagmented_library("adjacent accessible DNA")
FINAL_CAPTION = "Paired-end library from low-salt, antibody-tethered tagmentation of an adjacent nucleosome-depleted region."
SEQUENCING_INTRO = "Paired reads enter the short accessible-DNA fragment from the two targeted Tn5 insertion sites."

def sections():
    return [
        ("Target a transcription-associated nucleosome", antibody_enzyme_rows("pA/G–Tn5", "loaded Tn5 tethered beside an NDR"),
         "An H3K4me2/3 or RNA-polymerase antibody positions adapter-loaded transposase on chromatin flanking accessible DNA."),
        ("Activate in low salt", [Row(chunks=[("marked nucleosome | low-salt reach → [nucleosome-depleted region] ← low-salt reach | marked nucleosome", None, False)])],
         "Reduced ionic strength promotes tagmentation of the nearby accessible gap rather than ordinary nucleosomal CUT&Tag fragments."),
        ("Tag accessible DNA and amplify", tagmentation_scene("adjacent accessible DNA").rows(),
         "Tn5 cleavage and adapter transfer are coupled; ** marks both adapter–insert boundaries."),
    ]
