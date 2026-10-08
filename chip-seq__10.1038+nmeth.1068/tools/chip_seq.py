"""Robertson et al. conventional ChIP-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import antibody_enzyme_rows, illumina_ligation_library, ligated_insert_scene

TITLE = "ChIP-seq — conventional immunoprecipitation library"
NOTES = "01_chip-seq.html"
SOURCE = 'Defining protocol: <a href="https://doi.org/10.1038/nmeth.1068">Robertson et al., <i>Nature Methods</i> (2007)</a>.'
SUMMARY = "Crosslink and fragment chromatin, immunoprecipitate DNA associated with a chosen protein or histone mark, reverse crosslinks, and make a conventional ligation-based Illumina library."
CAVEAT = "The historical protocol predates current dual-index TruSeq kits. Its vendor adapter bases were not printed, so the page marks the canonical modern paired-end adapter layout as inferred."
FINAL_LIBRARY, SEQ_PRIMERS = illumina_ligation_library("ChIP-enriched DNA", inferred=True)
FINAL_CAPTION = "INFERRED — PCR-amplified sequencing library from ChIP-enriched DNA; dotted canonical arms are a modern readout for the historical vendor adapter."
SEQUENCING_INTRO = "Reads begin at the ends of the immunoprecipitated DNA fragment; paired-end sequencing additionally preserves fragment-span information."

def sections():
    return [
        ("Crosslink and fragment chromatin", [Row(chunks=[("protein—DNA chromatin → formaldehyde crosslink → sonicated fragments", None, False)])],
         "Mechanical shearing defines fragment ends independently of the target protein."),
        ("Immunoprecipitate the target", antibody_enzyme_rows("bead capture", "retain target-associated fragments"),
         "A target-specific antibody enriches the covalently crosslinked protein–DNA complexes; crosslinks are then reversed."),
        ("Repair, dA-tail and ligate adapters", ligated_insert_scene("ChIP-enriched DNA").rows(),
         "INFERRED — canonical TruSeq placement represents the historical vendor Y-adapter; ** marks both adapter–insert ligations."),
    ]
