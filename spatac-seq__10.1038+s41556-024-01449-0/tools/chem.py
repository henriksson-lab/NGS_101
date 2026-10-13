"""SPATAC-seq protocol page model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import split_pool_atac as spa


TITLE = "SPATAC-seq — split-pool ligation-based single-cell ATAC-seq"
NOTES = "01_spatac_seq.html"
SOURCE = ('<a href="https://doi.org/10.1038/s41556-024-01449-0">Sun et al., '
          '<i>Nature Cell Biology</i> (2024)</a>, Figure 1 and Supplementary Table 1; '
          '<a href="https://github.com/Lan-lab/SPATAC-seq">authors\' processing code</a>.')
SUMMARY = ("A two-ended barcode is installed during custom-Tn5 transposition, then "
           "successive split-pool ligations append a second barcode to the Read-1 end "
           "and a third barcode to the Read-2 end. Four 8-nt pieces identify one nucleus.")
CAVEAT = ("SPATAC-seq is a split-pool single-cell method, not Spatial-ATAC-seq. The "
          "representative sequences below are one member of each published oligo panel; "
          "the linked research note identifies the complete tables and whitelist.")
FINAL_LIBRARY, SEQ_PRIMERS = spa.library("accessible genomic DNA")
READ_LENGTHS = {
    "Read 1": 115,
    "Index 1 (i7)": 8,
    "Index 2 (i5)": 8,
    "Read 2": 100,
}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("The hybrid final library uses a TruSeq Read 1 site on the left and a "
                 "Nextera S7/Read 2 site on the right. ** marks the two split-pool "
                 "ligation bonds. Barcode pieces are structured parts of one 32-nt "
                 "cell identifier; i5 and i7 identify the PCR sublibrary.")
SEQUENCING_INTRO = ("The authors' demonstrated preprocessing retains Read 1 cycles "
                    "66–115 and Read 2 cycles 51–100 as the 50-nt genomic mates. Before "
                    "those positions, Read 1 reports round 2 and the forward half of "
                    "round 1; Read 2 reports round 3 and the reverse half of round 1. "
                    "The two eight-cycle index reads report the PCR sublibrary.")


def oligos():
    return spa.oligos()


def workflow():
    return spa.workflow("accessible genomic DNA")


def final_rows():
    return spa.final_rows(FINAL_LIBRARY)


def protocol_graph():
    return spa.protocol_graph(
        "accessible genomic DNA",
        authority="Defining source: https://doi.org/10.1038/s41556-024-01449-0")
