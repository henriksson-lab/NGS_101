"""Spatial-ATAC-seq orthogonal in-tissue barcoding."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import spatial_epigenome as se


TITLE = "Spatial-ATAC-seq — spatial chromatin accessibility"
NOTES = "01_spatial-atac-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41586-022-05094-1">'
          'Deng et al., <i>Nature</i> (2022)</a>; '
          '<a href="https://doi.org/10.1038/s41596-024-01013-y">Farzad et al., '
          '<i>Nature Protocols</i> (2024)</a> supplies the exact oligo tables.')
SUMMARY = ("Tn5 tags accessible chromatin in a fixed tissue section. Two perpendicular "
           "microfluidic ligations append 8-nt A and B barcodes whose ordered pair "
           "identifies a tissue pixel.")
ATAC_TN5_ME_A = "/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG"
FINAL_LIBRARY, SEQ_PRIMERS = se.library("accessible genomic fragment")
READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("The custom-primer end reads barcode B, linker 2, barcode A, linker 1 "
                 "and then accessible genomic DNA. ** marks the two in-tissue ligations.")
SEQUENCING_INTRO = ("Paired-end 150 sequencing was used. Raw Read 1 carries barcode B "
                    "at cycles 1–8 and barcode A at cycles 39–46 before genomic DNA; "
                    "Read 2 reports the opposite genomic end.")


def oligos():
    return (f"<b>Spatial-ATAC Tn5ME-A</b> — 5′-{ATAC_TN5_ME_A}-3′ "
            "(Deng et al. 2022 Methods; distinct from the CUT&amp;Tag adapter in Farzad et al. Table 1)",
            *se.oligos())


def workflow():
    return se.workflow(
        "accessible genomic fragment",
        first_action="Expose the tissue section to linker-loaded Tn5 and activate tagmentation",
        first_note="Tn5 transfers ligation-linker 1 and the S7-side adapter into accessible chromatin.",
    )


def final_rows():
    return se.final_rows(FINAL_LIBRARY)


def protocol_graph():
    return se.protocol_graph(
        "accessible genomic fragment", title="Spatial-ATAC-seq",
        authority="Defining source: https://doi.org/10.1038/s41586-022-05094-1; exact oligos: https://doi.org/10.1038/s41596-024-01013-y",
        antibody_targeted=False)
