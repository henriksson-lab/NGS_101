"""Ramani et al. 2017 sci-Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment
from restriction import DPNII
from single_cell_hic import ContactWorkflow, inferred_truseq_library
WORKFLOW = ContactWorkflow("sci-Hi-C", DPNII, "combinatorial indexing", "bridge-adaptor", "PCR", barcoding_rounds=2)
TITLE = "sci-Hi-C — two-level combinatorial indexing of contacts"
NOTES = "01_sci-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4155">Ramani et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "A first well barcode enters through indexed bridge adaptors before proximity ligation; a second well barcode enters through indexed library adaptors after redistribution."
STEPS = (("Digest in 96 wells", "Distribute fixed nuclei and digest chromatin with DpnII."),
         ("Barcode round 1", "Ligate biotinylated indexed bridge adaptors to restriction ends in each first-round well."),
         ("Pool, ligate, redistribute", "Pool nuclei, proximity-ligate bridged ends, then sort a small number of nuclei into second-round wells."),
         ("Barcode round 2 and PCR", "Add indexed Illumina Y adaptors, enrich biotin-bearing molecules and PCR-amplify; the barcode pair identifies one nucleus."))
READOUT = "The two combinatorial indices assign each read pair to a nucleus; paired genomic ends report a chromatin contact."
JUNCTION_CAPTION = "The contact junction contains the first-round indexed bridge (BC1); ** marks each ligation boundary."
LIBRARY_CAPTION = "BC1 remains internal while the second-round indexed Illumina adapter supplies the library index and sequencing-primer sites."
JUNCTIONS = (("locus A", "DpnII end A", "restriction end"),
             ("DpnII end A", "round-1 barcode A", "adapter ligation"),
             ("round-1 barcode B", "DpnII end B", "adapter ligation"))
def contact_product():
    return Construct([Segment("locus A", "X"*22, placeholder=True), Segment("DpnII end A", "GATC", "me"),
        Segment("round-1 barcode A", "B"*8, "cbc", placeholder=True), Segment("bridge", "N"*12, "umi", placeholder=True),
        Segment("round-1 barcode B", "B"*8, "cbc", placeholder=True), Segment("DpnII end B", "GATC", "me"),
        Segment("locus B", "X"*22, placeholder=True)], name="sci-Hi-C bridged contact")
def final_library(): return inferred_truseq_library(contact_product(), "sci-Hi-C two-index library")
