"""Tan et al. 2018 Dip-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import MBOI
from chemdraw import Construct, Segment
from single_cell_hic import ContactWorkflow, unresolved_illumina_library
WORKFLOW = ContactWorkflow("Dip-C", MBOI, "after ligation", "none", "META")
TITLE = "Dip-C — contact ligation followed by META amplification"
NOTES = "01_dip-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1126/science.aat5641">Tan et al., <i>Science</i> (2018)</a>; detailed protocol: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8225968/">Tan et al. (2021)</a>.'
SUMMARY = "Dip-C omits biotin capture and uses multiplex end-tagging amplification (META), a transposon-based single-cell WGA route."
STEPS = (("Digest and ligate", "Process fixed chromatin with a compatible restriction workflow (MboI is drawn) and proximity-ligate without biotin fill-in."),
         ("Isolate one nucleus", "Sort a single ligated nucleus, lyse it and reverse crosslinks."),
         ("META", "Tag many molecule ends by transposition and extend/amplify them as whole-genome material."),
         ("Sequencing library", "Convert the amplified material to a paired-end indexed library for contact calling and haplotype imputation."))
READOUT = "Paired reads spanning chimeric products report contacts; dense single-cell coverage supports diploid haplotype reconstruction."
JUNCTION_CAPTION = "Dip-C retains ligation products without a biotin-selected fill junction. ** marks ligation."
LIBRARY_CAPTION = "A representative contact molecule after META. The defining paper identifies 39 transposon-derived bases at each read start but the accessible text does not disclose their sequence or a complete outer adapter structure."
LIBRARY_CAVEAT = "INFERRED — the 39-base META read prefixes are source-defined by length and role only; their bases and the complete sequencing shell remain unresolved."
SEQUENCING_UNAVAILABLE = "Dip-C reports paired-end sequencing and removal of the first 39 transposon-derived bases from both reads, but the defining source available here does not print the primer or transposon oligos needed for base-level placement."
def contact_product(): return WORKFLOW.junction()
META_READ_PREFIX_NT = 39
def final_library():
    tagged = Construct([
        Segment("META transposon prefix", "T"*META_READ_PREFIX_NT, "me", placeholder=True),
        *list(contact_product()),
        Segment("opposite META transposon prefix", "T"*META_READ_PREFIX_NT, "me", placeholder=True),
    ], name="META-tagged contact")
    return unresolved_illumina_library(tagged, "Dip-C META-derived library", indexed=None)
