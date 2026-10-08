"""Nagano et al. 2013 single-cell Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from restriction import BGLII
from single_cell_hic import ContactWorkflow, inferred_truseq_library
from chemdraw import Construct, Segment

WORKFLOW = ContactWorkflow("Nagano single-cell Hi-C", BGLII,
    "after bulk proximity ligation", "biotin-fill", "PCR", capture=True)
TITLE = "Nagano single-cell Hi-C — ligate first, isolate nuclei second"
NOTES = "01_nagano-single-cell-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nature12593">Nagano et al., <i>Nature</i> (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2015.127">Nagano et al. (2015)</a>.'
SUMMARY = "BglII contact ends are biotin-filled and ligated in fixed nuclei; individual nuclei are isolated only after ligation."
STEPS = (("Crosslink and digest", "Digest fixed chromatin in intact nuclei with BglII."),
         ("Fill and ligate", "Fill the 5′ overhang with biotin-14-dATP, then proximity-ligate."),
         ("Isolate one nucleus", "Dilute and sort individual ligated nuclei before reversal of crosslinks."),
         ("Capture and amplify", "Fragment DNA, enrich biotin-bearing junctions on streptavidin and PCR-amplify an Illumina library."))
READOUT = "Paired reads start in the two genomic fragments on either side of a captured contact junction."
JUNCTION_CAPTION = "BglII-compatible ends form the contact boundary; the experimental junction is selected through the biotin-bearing fill-in. ** marks ligation."
LIBRARY_CAPTION = "PCR-completed paired-end library. Dotted adapter arms are canonical TruSeq structure because the paper names the Illumina workflow without printing those bases."
JUNCTIONS = (("left filled end", "right filled end", "proximity ligation"),)
def contact_product():
    j = WORKFLOW.digest().fill_in(biotin_base="A").junction
    return Construct([
        Segment("locus A", "X"*22, placeholder=True),
        Segment("junction before top biotin", j[:2], "me"),
        Segment("top-strand biotin-dA", j[2], "w1"),
        Segment("left filled end", j[3:5], "me"),
        Segment("right filled end", j[5:7], "me"),
        Segment("opposite-strand biotin-dA", j[7], "w1"),
        Segment("junction after opposite biotin", j[8:], "me"),
        Segment("locus B", "X"*22, placeholder=True),
    ], name="BglII biotin-filled contact")
def final_library(): return inferred_truseq_library(contact_product(), "Nagano scHi-C library")
