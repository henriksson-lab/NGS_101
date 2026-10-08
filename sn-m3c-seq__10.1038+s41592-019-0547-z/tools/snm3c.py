"""Lee et al. 2019 sn-m3C-seq chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from restriction import DPNII
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
WORKFLOW = ContactWorkflow("sn-m3C-seq", DPNII, "after ligation", "none", "bisulfite-PCR")
TITLE = "sn-m3C-seq — contacts and DNA methylation from one nucleus"
NOTES = "01_sn-m3c-seq.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/s41592-019-0547-z">Lee et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY = "Single-nucleus chromosome-conformation capture is followed by bisulfite conversion and snmC-seq2-style indexed random-priming library construction."
STEPS = (("Digest and ligate nuclei", "Generate proximity-ligation products in fixed nuclei without a biotin pull-down."),
         ("Isolate one nucleus", "Sort nuclei so the contact and methylation observations retain a common cellular identity."),
         ("Bisulfite conversion", "Convert unmethylated cytosine to uracil while methylated cytosine remains protected."),
         ("Indexed random priming and PCR", "Copy converted fragments with an indexed P5-bearing random primer, add the opposite adapter and amplify."))
READOUT = "Paired bisulfite reads jointly encode the two contact fragments and cytosine-conversion state."
JUNCTION_CAPTION = "The unconverted contact architecture is shown here; after ligation, bisulfite treatment changes sequence content but not the contact boundary. ** marks ligation."
LIBRARY_CAPTION = "Bisulfite-compatible snmC-seq2 product. The 25-base random-primer-derived prefix and 3-base Adaptase tail are experimentally trimmed; the sn-m3C paper does not print the complete adapter oligos."
LIBRARY_CAVEAT = "INFERRED — the trim-defined protocol-specific regions are established, but the complete outer adapter bases are not reported in the defining sn-m3C source."
SEQUENCING_UNAVAILABLE = "The source specifies paired-end reads and trimming 25 bases from the start and 3 bases from the end of each read, but does not print a complete sequencing-primer set on which base-level placement can be verified."
def contact_product(): return WORKFLOW.junction()
RANDOM_PRIMER_TRIM_NT = 25
ADAPTASE_TAIL_TRIM_NT = 3
def final_library():
    converted = Construct([
        Segment("indexed random-primer-derived sequence", "N"*RANDOM_PRIMER_TRIM_NT,
                "cbc", placeholder=True),
        *list(contact_product()),
        Segment("Adaptase low-complexity tail", "N"*ADAPTASE_TAIL_TRIM_NT,
                "me", placeholder=True),
    ], name="converted and tagged contact")
    return unresolved_illumina_library(converted, "sn-m3C-seq bisulfite library", indexed=True)
