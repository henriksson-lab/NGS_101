"""Robertson et al. 2007 ChIP-seq and original single-read Solexa library."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment, revcomp
from chromatin_epigenetics import antibody_enzyme_rows


def seg(name, top, tag=None, **kw):
    return Segment(name, top, tag, **kw)


TITLE = "ChIP-seq — conventional immunoprecipitation library"
NOTES = "01_chip-seq.html"
SOURCE = 'Defining protocol: <a href="https://doi.org/10.1038/nmeth.1068">Robertson et al., <i>Nature Methods</i> (2007)</a>.'
SUMMARY = "Crosslink and sonicate chromatin, immunoprecipitate target-associated DNA, reverse crosslinks, and prepare the recovered fragments for the original single-read Solexa Genome Analyzer."
CAVEAT = "This 2007 protocol predates TruSeq and sample indexing. It used Illumina’s Genomic DNA Sample Prep adapter and PCR primers 1.1/2.1, followed by one genomic-DNA sequencing read—not a modern paired-end dual-index endpoint."

GENOMIC_ADAPTER_BOTTOM = il.TRUSEQ_READ1
GENOMIC_PCR2 = "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"
# PCR primer 2.1 has one 3'-terminal T paired to the insert's dA junction; the rest
# determines the phosphorylated adapter oligo rather than duplicating it independently.
GENOMIC_ADAPTER_TOP = revcomp(GENOMIC_PCR2)[1:]


def adapter_scene() -> Scene:
    sc = Scene()
    sc.strand("bottom", [
        seg("Read 1 fork arm", GENOMIC_ADAPTER_BOTTOM[:-13], "r1"),
        seg("genomic stem complement", GENOMIC_ADAPTER_BOTTOM[-13:-1], "r1"),
        seg("3′ T", GENOMIC_ADAPTER_BOTTOM[-1], "r1"),
    ], label="Genomic Adapter Oligo")
    sc.anneal("top", [
        seg("genomic stem", GENOMIC_ADAPTER_TOP[:12], "p7"),
        seg("genomic fork arm", GENOMIC_ADAPTER_TOP[12:], "p7"),
    ], to="bottom", pair=("genomic stem", "genomic stem complement"),
              label="Genomic Adapter Oligo", mod5="p",
              unpaired=("genomic fork arm",))
    sc.mark("bottom", "3′ T", "ligation overhang")
    return sc


SEQ_PRIMERS = (sp.custom(
    "Read 1", "Genomic DNA sequencing primer", il.TRUSEQ_READ1,
    'Illumina "Genomic DNA Sample Prep Kit" / adapter-sequence guide'),)


def final_library() -> Construct:
    lib = Construct([
        seg("P5 / PCR primer 1.1", il.P5, "p5"),
        seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("ChIP-enriched DNA", "X" * 36, placeholder=True),
        seg("dA junction", "A"),
        seg("Genomic adapter / PCR primer 2.1 arm", GENOMIC_ADAPTER_TOP, "p7"),
    ], name="2007 ChIP-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=("Read 1",))
    if problems:
        raise ValueError("invalid original ChIP-seq library: " + "; ".join(problems))
    return lib


FINAL_LIBRARY = final_library()
FINAL_CAPTION = "Original unindexed single-read Solexa library made with Genomic Adapter Oligo Mix and PCR primers 1.1/2.1."
SEQUENCING_INTRO = "The defining experiment generated single-end approximately 27-nt reads. It had no sample index and no Read 2."


def sections():
    return [
        ("Crosslink and fragment chromatin",
         [Row(chunks=[("protein—DNA chromatin → formaldehyde crosslink → sonicated fragments", None, False)])],
         "Mechanical shearing defines fragment ends independently of the target protein."),
        ("Immunoprecipitate the target",
         antibody_enzyme_rows("bead capture", "retain target-associated fragments"),
         "A target-specific antibody enriches crosslinked protein–DNA complexes; crosslinks are then reversed."),
        ("Repair, dA-tail and ligate the Genomic Adapter", adapter_scene().rows(),
         "The original Illumina single-read fork has a 3′-T ligation overhang. PCR primers 1.1 and 2.1 complete the flow-cell arms."),
    ]
