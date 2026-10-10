"""Stergachis et al. Fiber-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from longread_footprinting import conventional_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows

TITLE = "Fiber-seq — long-read chromatin fibers"
NOTES = "01_fiber-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.aaz1646">Stergachis et al., <i>Science</i> (2020)</a>.'
SUMMARY = "Hia5 methylates exposed adenines in intact nuclei; PCR-free PacBio sequencing preserves that m6dA stencil together with the underlying long DNA molecule."
CAVEAT = "The defining paper establishes Hia5 marking and PacBio single-molecule readout. The page uses the current standard PCR-free SMRTbell geometry; hairpin bases are proprietary and are not invented."
SMRTBELL = conventional_smrtbell()
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL)

def sections():
    return [
        ("Stencil accessible DNA with Hia5", marked_fiber("Hia5"),
         "SAM-dependent m6dA is deposited where chromatin does not protect adenines."),
        ("Extract high-molecular-weight marked DNA", [Row(chunks=[("one native, m6dA-bearing DNA fiber", "w1", False)])],
         "No amplification is introduced: the modification pattern must remain on the sequenced template."),
        ("Shear, repair, dA-tail and ligate PacBio hairpins", smrtbell_rows(SMRTBELL),
         "Both hairpins close the duplex into a nuclease-resistant SMRTbell; the drawing marks both sealed insert-adapter junctions."),
    ]
