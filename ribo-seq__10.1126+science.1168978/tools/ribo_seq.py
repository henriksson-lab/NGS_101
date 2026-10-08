"""Ingolia ribosome-profiling library with the 2012 indexed run geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row
from rna_special import adapter_ligation_scene, circular_cdna_rows

TITLE = "Ribo-seq — ribosome profiling"
NOTES = "01_ribo-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1126/science.1168978">Ingolia et al. (2009)</a>; library details: <a href="https://doi.org/10.1038/nprot.2012.086">Ingolia et al. (2012)</a>.'
SUMMARY = "Nuclease leaves ribosome-protected RNA footprints. Size-selected footprints receive a 3′ adapter, are reverse-transcribed and circularized, then PCR adds the Illumina run structure."
CAVEAT = "The six-base sample index is selected from the published indexed reverse-primer table; the model leaves only that position variable."

LINKER = "CTGTAGGCACCATCAAT"
RT_PRIMER = ("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGC"
             "CACTCA"
             "TTCAGACGTGTGCTCTTCCGATCTATTGATGGTGCCTACAG")
FORWARD_PCR = il.P5

# The P5 PCR primer and Read 1 site share ACAC physically.  Store the overlap once.
FINAL_LIBRARY = Construct([
    seg("P5-exclusive", il.P5[:-4], "p5"),
    seg("P5 / Read 1 shared ACAC", il.P5[-4:], "r1"),
    seg("Read 1 site after shared ACAC", il.TRUSEQ_READ1[4:], "r1"),
    seg("ribosome footprint cDNA", "X" * 30, placeholder=True),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="single-index ribosome-profiling library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)")):
    raise ValueError("Ribo-seq final library: " + "; ".join(problems))

FINAL_CAPTION = "PCR-linearized product from the circular cDNA: single-end, six-base single-index, with the physical P5/Read 1 ACAC overlap represented once."
SEQUENCING_INTRO = "The protocol explicitly uses the standard Illumina genomic Read 1 and indexing primers. Read 1 enters the footprint; one six-base i7 read identifies the sample."

def sections():
    return [
        ("Digest unprotected RNA", [Row(chunks=[("mRNA ——— [80S ribosome protects ~28–30 nt] ——— mRNA", None, False)]), Row(chunks=[("RNase I ↓          retained ribosome footprint", "umi", False)])], "RNase digestion and monosome isolation select ribosome-protected fragments."),
        ("Ligate the footprint 3′ adapter", adapter_ligation_scene(fragment_name="ribosome footprint").rows(), "A preadenylated adapter is joined to the footprint 3′ end before reverse transcription."),
        ("Reverse-transcribe and circularize cDNA", circular_cdna_rows([seg("footprint cDNA", "X"*30, placeholder=True), seg("published RT-primer handles", "X"*24, "r2", placeholder=True)]), "The phosphorylated RT primer contains the sequencing/PCR handles separated by two C18 spacers. CircLigase closes first-strand cDNA, providing the second PCR boundary without a second RNA ligation."),
    ]
