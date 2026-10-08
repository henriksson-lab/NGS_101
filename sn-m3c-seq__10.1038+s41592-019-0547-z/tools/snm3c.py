"""Lee et al. 2019 sn-m3C-seq chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from restriction import DPNII
from single_cell_hic import ContactWorkflow, inferred_truseq_library
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
LIBRARY_CAPTION = "Representative bisulfite-compatible paired-end product. The canonical outer adapter shell is inferred; the protocol-specific random-priming transition is called out above."
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
    return inferred_truseq_library(converted, "sn-m3C-seq bisulfite library")
