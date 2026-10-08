"""Miura et al. amplification-free PBAT model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Segment
from chromatin_epigenetics import pbat_priming_scene

TITLE = "PBAT — post-bisulfite adaptor tagging"
NOTES = "01_pbat.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gks454">Miura et al., <i>Nucleic Acids Research</i> (2012)</a>.'
SUMMARY = "Bisulfite-convert genomic DNA first, then install the two sequencing handles by successive random-primer extensions, avoiding conversion-induced loss of completed library molecules."
BIOPEA2 = il.TRUSEQ_READ1
PE_REVERSE = il.P7
PRIMER3 = il.P5 + il.TRUSEQ_READ1[4:]
FINAL_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"), Segment("Read 1 arm remainder", BIOPEA2[4:], "r1"),
    Segment("first random N4", "N" * 4, placeholder=True),
    Segment("directional PBAT insert", "X" * 36, placeholder=True),
    Segment("second random N4 reverse complement", "N" * 4, placeholder=True),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
], name="amplification-free PBAT sequencing template")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("invalid PBAT read layout: " + "; ".join(problems))
FINAL_CAPTION = "Amplification-free PBAT template after Primer-3 extension. The two N4 tracts record the random priming sites; the original protocol has no sample index."
SEQUENCING_INTRO = "The original GAIIx/HiSeq protocol produces a single Read 1 from the strand complementary to bisulfite-converted DNA; reads are consequently generally G-poor."

def sections():
    return [
        ("Bisulfite-convert untagged genomic DNA", [Row(chunks=[("genomic duplex → denature + bisulfite → fragmented C→U single strands", "w1", False)])],
         "Conversion precedes adapter tagging, the defining reversal relative to conventional WGBS."),
        ("First random-primer extension and capture", pbat_priming_scene("first-strand", "BioPEA2", BIOPEA2, mod5="biotin").rows(),
         "A 5′-biotinylated adapter primer ending in random N4 copies the converted strand; streptavidin immobilizes the product."),
        ("Second random-primer extension", pbat_priming_scene("second-strand", "PE-reverse", PE_REVERSE).rows(),
         "After alkaline denaturation, PE-reverse-N4 copies the bead-bound first strand; Primer-3 then adds P5 and completes the sequencing template."),
    ]
