"""RAMPAGE promoter-profiling workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, feature, revcomp

TITLE = "RAMPAGE — paired-end promoter profiling"
NOTES = "01_rampage.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1101/gr.139618.112">Batut et al. (2013)</a>; detailed protocol: <a href="https://doi.org/10.1002/0471142727.mb25b11s104">Batut &amp; Gingeras (2013)</a>.'
SUMMARY = "Template switching tags 5′-complete cDNA, cap trapping enriches genuine capped molecules, and paired-end sequencing connects a base-resolution transcription start site to downstream transcript sequence."
CAVEAT = "RAMPAGE was selected over generic CAGE because its defining paper and complete public protocol specify the paired-end molecular path and exact oligos. The six-base inline barcode is shown as a placeholder because its identity varies by sample."
TSO_HANDLE = "TAGTCGAACTGAAGGTCTCCAGCA"
RT_HANDLE = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
R1_PRIMER = TSO_HANDLE
R2_PRIMER = RT_HANDLE
FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("custom Read 1 site", R1_PRIMER, "r1"),
    seg("inline sample barcode", "B"*6, "cbc", placeholder=True,
        feature=feature("sample_inline", "sample_index", "whitelist", whitelist="published RAMPAGE TSO set")),
    seg("5-prime-complete cDNA", "X"*42, placeholder=True),
    seg("custom Read 2 site reverse complement", revcomp(R2_PRIMER), "r2"),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="RAMPAGE library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "RAMPAGE custom Read 1", R1_PRIMER,
              "Batut & Gingeras 2013, Synthetic oligonucleotides"),
    sp.custom("Read 2", "RAMPAGE SBS8 Read 2", R2_PRIMER,
              "Batut & Gingeras 2013, Synthetic oligonucleotides"),
)
if sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1", "Read 2")):
    raise ValueError("RAMPAGE custom sequencing primers do not land on the final library")
FINAL_CAPTION = "Paired-end RAMPAGE library: one insert end marks the capped transcription start and the mate samples downstream transcript sequence."
SEQUENCING_INTRO = "Read 1 identifies the transcription start boundary; the mate supplies transcript connectivity."

def sections():
    rna = [seg("RNA body", "X"*36, placeholder=True)]
    rt = [seg("RT / Read 2 handle", RT_HANDLE, "r2"),
          seg("random 15-mer", "N"*15, placeholder=True)]
    sc = Scene(); sc.strand("RNA", rna, label="ribosome-depleted total RNA")
    sc.anneal("random RT primer", rt, to="RNA", pair=("random 15-mer", "RNA body"),
              shift=12, label="rampage_RT", unpaired=("RT / Read 2 handle",))
    sc.arrow("random RT primer", "reverse transcription toward the RNA 5-prime end")
    sc.strand("TSO", [seg("Read 1 handle", TSO_HANDLE, "r1"),
                       seg("inline barcode", "B"*6, "cbc", placeholder=True),
                       seg("rGrGrG", "GGG")], label="rampage_TS")
    return [
        ("Synthesize 5′-complete cDNA", sc.rows(), "A random 15-mer carrying the Read 2 arm primes RT. At a capped RNA 5′ end, the barcoded rGrGrG-ended TSO installs the custom Read 1 arm."),
        ("Cap-trap full-length RNA–cDNA hybrids", [Row(chunks=[("m7G cap—RNA:cDNA—poly(A)  → cap oxidation/biotinylation", None, False)]), Row(chunks=[("       streptavidin retains capped, 5′-complete hybrids", "umi", False)])], "Cap trapping adds a second, independent selection for 5′ completeness."),
        ("Remove incomplete molecules and release cDNA", [Row(chunks=[("uncapped / incomplete RNA → phosphate-dependent exonuclease", None, False)]), Row(chunks=[("captured 5′-complete cDNA → released for library PCR", "r1", False)])], "Enzymatic depletion and cap capture create the high-specificity promoter boundary."),
        ("Add paired-end run structure", [Row(chunks=[("P5 / custom Read 1 — inline barcode — 5′-complete cDNA — custom Read 2 / P7", None, False)])], "Final PCR adds the flow-cell ends. RAMPAGE uses the TSO barcode inline at the start of Read 1 rather than a separate index read."),
    ]
