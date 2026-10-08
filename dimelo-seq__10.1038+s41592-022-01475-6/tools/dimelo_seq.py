"""Altemose et al. DiMeLo-seq directed-methylation model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Row
from chromatin_epigenetics import dimelo_mark_rows

TITLE = "DiMeLo-seq — antibody-directed methylation on long DNA"
NOTES = "01_dimelo-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41592-022-01475-6">Altemose et al., <i>Nature Methods</i> (2022)</a>.'
SUMMARY = "An antibody recruits Protein A–Hia5 to a chromatin target; SAM activation deposits nearby m6A marks that are read directly on unamplified long DNA molecules."
CAVEAT = "The scientific construct is the marked native DNA molecule. Oxford Nanopore library-adapter bases are kit-specific and are intentionally represented by role rather than invented sequence."
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "Extract high-molecular-weight native DNA, attach the chosen Oxford Nanopore ligation-sequencing motor adapter, and feed one strand through the pore. Modified-base calling reports exogenous m6A together with endogenous CpG methylation on the same molecule; no sequencing primer binds."

def sections():
    return [
        ("Bind antibody and pA–Hia5", dimelo_mark_rows()[:4],
         "The antibody supplies target specificity while the Protein A fusion positions Hia5 near the chromatin feature."),
        ("Activate directed methylation", dimelo_mark_rows()[4:],
         "SAM activates Hia5; deposited m6A is a covalent mark on nearby adenines, not an adapter or cleavage site."),
        ("Extract native long DNA", [Row(chunks=[("intact marked chromatin DNA → deproteinize gently → long native molecules", None, False)])],
         "The DNA is not PCR-amplified, preserving co-occurring marks and haplotype information along each molecule."),
    ]
