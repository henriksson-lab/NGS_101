"""Spatial-CUT&Tag / hsrChST-seq molecular model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import spatial_epigenome as se


TITLE = "Spatial-CUT&Tag (hsrChST-seq)"
NOTES = "01_spatial-cut-tag.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1126/science.abg7216">'
          'Deng et al., <i>Science</i> (2022)</a>; first named hsrChST-seq in the '
          '<a href="https://doi.org/10.1101/2021.03.11.434985">2021 preprint</a>. '
          'Exact oligos: <a href="https://doi.org/10.1038/s41596-024-01013-y">'
          'Farzad et al., <i>Nature Protocols</i> (2024)</a>.')
SUMMARY = ("Antibody-tethered pA–Tn5 transfers linker-bearing adapters beside a chosen "
           "histone mark. Two perpendicular in-tissue ligations append 8-nt A and B "
           "barcodes whose ordered pair identifies a tissue pixel.")
FINAL_LIBRARY, SEQ_PRIMERS = se.library("histone-mark-associated genomic fragment")
READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("The custom-primer end reads barcode B, linker 2, barcode A, linker 1 "
                 "and then genomic DNA. ** marks the two in-tissue ligations.")
SEQUENCING_INTRO = ("HiSeq 4000 paired-end 150 sequencing was used. Raw Read 1 carries "
                    "barcode B at cycles 1–8 and barcode A at cycles 39–46 before "
                    "entering genomic DNA; Read 2 reports the opposite genomic end.")


def oligos():
    return se.oligos()


def workflow():
    return se.workflow(
        "histone-mark-associated genomic fragment",
        first_action="Bind primary/secondary antibody, recruit pA–Tn5 and activate tagmentation",
        first_note="Antibody-tethered pA–Tn5 transfers ligation-linker 1 at chromatin carrying the selected histone modification.",
    )


def final_rows():
    return se.final_rows(FINAL_LIBRARY)


def protocol_graph():
    return se.protocol_graph(
        "histone-mark-associated genomic fragment", title="Spatial-CUT&Tag",
        authority="Defining source: https://doi.org/10.1126/science.abg7216; exact oligos: https://doi.org/10.1038/s41596-024-01013-y",
        antibody_targeted=True)
