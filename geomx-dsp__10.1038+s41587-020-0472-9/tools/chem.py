"""GeoMx DSP NGS readout page model."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import digital_spatial_tags as dsp


TITLE = "GeoMx DSP — photocleaved spatial tags with NGS readout"
NOTES = "01_geomx_dsp.html"
SOURCE = ('Defining method: <a href="https://doi.org/10.1038/s41587-020-0472-9">'
          'Merritt et al., <i>Nature Biotechnology</i> (2020)</a>. Library preparation: '
          '<a href="https://www.brukerspatialbiology.com/wp-content/uploads/2022/06/'
          'MAN-10153-01-GeoMx-DSP-NGS-Readout-User-Manual.pdf">GeoMx DSP NGS Readout '
          'User Manual MAN-10153-01</a>.')
SUMMARY = ("Patterned UV releases target-identifying synthetic tags from a selected "
           "tissue area. A well-specific unique-dual-index PCR preserves that AOI, while "
           "a 12-nt RTS ID identifies the RNA or protein target and a 14-nt UMI counts molecules.")
CAVEAT = ("The vendor publishes the construct order and exact lengths but not the tag or "
          "Seq Code bases. Canonical TruSeq bases are shown with inferred styling solely "
          "to explain standard-primer binding; UMI, RTS-ID and UDI values remain placeholders.")
FINAL_LIBRARY, SEQ_PRIMERS = dsp.library()
READ_LENGTHS = {"Read 1": 27, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 27}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("Vendor-published final order and lengths. Inferred canonical adapter "
                 "bases are visibly marked; the central 14-nt UMI and 12-nt target ID, "
                 "and the two AOI indexes, are structured features with undisclosed values.")
SEQUENCING_INTRO = ("MAN-10153-01 specifies paired-end 2×27 plus two 8-cycle index "
                    "reads. Read 1 sees UMI then RTS ID; Read 2 sees RTS ID then UMI. "
                    "Index 2 orientation depends on the Illumina instrument workflow.")


def oligos():
    return ()


def workflow():
    return dsp.workflow()


def final_rows():
    return dsp.final_rows(FINAL_LIBRARY)
