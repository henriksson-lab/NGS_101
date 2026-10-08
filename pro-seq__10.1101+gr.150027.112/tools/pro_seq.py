"""Mahat et al. PRO-seq library, including its TruSeq Small RNA run geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row
from rna_special import adapter_ligation_scene

TITLE = "PRO-seq — precision nuclear run-on sequencing"
NOTES = "01_pro-seq.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1101/gr.150027.112">Kwak et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1038/nprot.2016.086">Mahat et al. (2016)</a>.'
SUMMARY = "Engaged polymerases extend nascent RNA by one or a few biotin-NTPs. Biotin purification and sequential RNA-adapter ligations preserve the nucleotide-resolution 3′ end for sequencing."
CAVEAT = "The six-base sample index varies between RPI-n primers; the model leaves only those six bases variable."

# Mahat et al. 2016, Table 1. RNA U is represented as DNA T in the final cDNA.
VRA3 = "GATCGTCGGACTGTAGAACTCTGAAC"
VRA5 = "CCTTGGCACCCGAGAATTCCA"
RP1 = "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RPI_BIND = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
SMALL_RNA_R1 = "GTTCAGAGTTCTACAGTCCGACGATC"
SMALL_RNA_I1 = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("TruSeq Small RNA Read 1 site", SMALL_RNA_R1, "r1"),
    seg("reverse-complement nascent RNA", "X" * 34, placeholder=True),
    seg("VRA5-derived Read 2 side", "TGGAATTCTCGGGTGCCAAGG", "r2"),
    seg("Index 1 continuation", "AACTCCAGTCAC", "cbc"),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc", placeholder=True),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="single-index PRO-seq library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Small RNA Read 1", SMALL_RNA_R1,
              "Mahat et al. 2016 Table 1; TruSeq Small RNA geometry"),
    sp.custom("Index 1 (i7)", "TruSeq Small RNA Index 1", SMALL_RNA_I1,
              "Mahat et al. 2016 Table 1; RPI-n"),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)")):
    raise ValueError("PRO-seq final library: " + "; ".join(problems))

FINAL_CAPTION = "Mahat et al. single-end, six-base single-index library. Read 1 sequences the reverse complement of nascent RNA from its informative 3′ end."
SEQUENCING_INTRO = "The TruSeq Small RNA Read 1 primer enters the reverse-complement nascent-RNA insert; the custom small-RNA index primer reads one six-base i7 index. There is no i5 read."

def sections():
    return [
        ("Biotin nuclear run-on", [Row(chunks=[("DNA template ————— RNAP — nascent RNA—3′", None, False)]), Row(chunks=[("                            + biotin-NTP → RNA—●—3′ (stalled)", "umi", False)])], "Permeabilized cells undergo a short run-on with biotinylated NTPs; the bulky incorporated nucleotide both marks and limits extension."),
        ("Ligate the 3′ RNA adapter", adapter_ligation_scene(fragment_name="biotinylated nascent RNA").rows(), "The adapter is joined at the informative nascent-RNA 3′ end."),
        ("Capture, decap and ligate the 5′ adapter", [Row(chunks=[("biotin-RNA—3′ adapter  → streptavidin capture → 5′ end repair", None, False)]), Row(chunks=[("5′ adapter ** nascent RNA ** 3′ adapter", "r1", False)])], "Repeated streptavidin enrichment retains biotin-run-on RNA; 5′ processing permits the second RNA-adapter ligation."),
        ("Reverse-transcribe and PCR-index", [Row(chunks=[("RP1/P5 — reverse-complement nascent RNA — VRA5-derived arm — i7(6) — P7", None, True)])], "RP1 primes reverse transcription and is reused for PCR with one RPI-n primer, which installs a six-base i7 index."),
    ]
