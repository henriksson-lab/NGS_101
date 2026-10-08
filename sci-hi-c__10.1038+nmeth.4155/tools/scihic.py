"""Ramani et al. 2017 sci-Hi-C chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp
from restriction import DPNII
from single_cell_hic import ContactWorkflow
WORKFLOW = ContactWorkflow("sci-Hi-C", DPNII, "combinatorial indexing", "bridge-adaptor", "PCR", barcoding_rounds=2)
TITLE = "sci-Hi-C — two-level combinatorial indexing of contacts"
NOTES = "01_sci-hi-c.html"
CITATION = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4155">Ramani et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "A first well barcode enters through indexed bridge adaptors before proximity ligation; a second well barcode enters through indexed library adaptors after redistribution."
STEPS = (("Digest in 96 wells", "Distribute fixed nuclei and digest chromatin with DpnII."),
         ("Barcode round 1", "Ligate biotinylated indexed bridge adaptors to restriction ends in each first-round well."),
         ("Pool, ligate, redistribute", "Pool nuclei, proximity-ligate bridged ends, then sort a small number of nuclei into second-round wells."),
         ("Barcode round 2 and PCR", "Add indexed Illumina Y adaptors, enrich biotin-bearing molecules and PCR-amplify; the barcode pair identifies one nucleus."))
READOUT = "The first eight bases of each genomic read carry BC2; the two BC1 copies flank the internal bridge. The barcode pair assigns the contact to a nucleus."
JUNCTION_CAPTION = "The contact junction contains the first-round indexed bridge (BC1); ** marks each ligation boundary."
LIBRARY_CAPTION = "BC1 remains internal. BC2 is an 8-bp inline barcode at both read starts, followed by the source-printed CGT fixed bases; sci-Hi-C does not use an index read for BC2."
JUNCTIONS = (("locus A", "DpnII end A", "restriction end"),
             ("DpnII end A", "round-1 barcode A", "adapter ligation"),
             ("round-1 barcode B", "DpnII end B", "adapter ligation"))
def contact_product():
    return Construct([Segment("locus A", "X"*22, placeholder=True), Segment("DpnII end A", "GATC", "me"),
        Segment("round-1 barcode A", "B"*8, "cbc", placeholder=True), Segment("bridge EcoRI site", "GAATTC", "me"),
        Segment("round-1 barcode B", "B"*8, "cbc", placeholder=True), Segment("DpnII end B", "GATC", "me"),
        Segment("locus B", "X"*22, placeholder=True)], name="sci-Hi-C bridged contact")
BC2_EXAMPLE = "TGACCTTG"  # A1 adapter, Supplementary Data oligo table
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["R2"])
REQUIRED_ROLES = ("Read 1", "Read 2")
def final_library():
    read_prefix = BC2_EXAMPLE + "CGT"
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        Segment("BC2 A1", BC2_EXAMPLE, "cbc"), Segment("fixed read prefix", "CGT", "me"),
        *list(contact_product()),
        Segment("opposite fixed prefix and BC2", revcomp(read_prefix), "cbc"),
        Segment("Read 2 site", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-Hi-C two-barcode library")
    errors = sp.verify(lib, SEQ_PRIMERS, required_roles=REQUIRED_ROLES)
    if errors: raise ValueError("invalid sci-Hi-C library: " + "; ".join(errors))
    return lib
