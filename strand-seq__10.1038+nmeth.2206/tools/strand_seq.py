"""Falconer et al. Strand-seq parental-template model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import illumina_ligation_library, ligated_insert_scene, strand_selection_rows

TITLE = "Strand-seq — parental template-strand sequencing"
NOTES = "01_strand-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.2206">Falconer et al., <i>Nature Methods</i> (2012)</a>.'
SUMMARY = "Cells replicate once in BrdU; after single-cell isolation, Hoechst/UV selectively nicks the BrdU-containing nascent strands so only the two parental template strands enter the sequencing library."
CAVEAT = "The defining paper describes the strand-selection principle but uses historical Illumina library reagents. The final canonical paired-end arms are therefore shown as inferred."
FINAL_LIBRARY, SEQ_PRIMERS = illumina_ligation_library("retained parental template fragment", inferred=True)
FINAL_CAPTION = "INFERRED — library derived from one retained parental template strand. Dotted adapter regions are the canonical paired-end implementation."
SEQUENCING_INTRO = "Mapped read direction records which parental template strand each daughter cell inherited; this orientation signal is the assay readout."

def sections():
    return [
        ("Label nascent DNA during one S phase", strand_selection_rows()[:2],
         "Semiconservative replication leaves each chromatid with one parental template and one BrdU-containing newly synthesized strand."),
        ("Destroy the nascent strand", strand_selection_rows()[2:],
         "Hoechst sensitization and UV nick BrdU DNA; exonucleolytic/fragment-processing steps prevent the nascent strand from contributing a library."),
        ("Build the single-cell library", ligated_insert_scene("retained parental template fragment").rows(),
         "INFERRED — canonical TruSeq placement represents the historical adapter kit; ** marks the two adapter–insert ligations."),
    ]
