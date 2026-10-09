"""Original paired-end eCLIP molecular workflow and exact ENCODE SOP adapters."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
from rna_special import adapter_ligation_scene

TITLE = "eCLIP — enhanced crosslinking and immunoprecipitation"
NOTES = "01_eclip.html"
SOURCE = 'Defining paper and protocol: <a href="https://doi.org/10.1038/nmeth.3810">Van Nostrand et al. (2016)</a>; <a href="https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf">ENCODE eCLIP SOP v1.P</a>.'
SUMMARY = "UV-crosslinked RBP–RNA complexes are immunoprecipitated and RNase-trimmed. A 3′ RNA adapter precedes reverse transcription; a second adapter is ligated to cDNA so both read-through and crosslink-truncated cDNAs can amplify."
CAVEAT = "The SOP provides several color-balanced RNA barcode and D5/D7 index choices. The final panel uses the published X1A, D501 and D701 examples; their positions, widths and primer family do not vary."
I5_FEATURE=feature("sample_index_i5","sample_index","fixed")
I7_FEATURE=feature("sample_index_i7","sample_index","fixed")
INLINE_FEATURE=feature("inline_rna_barcode","inline_barcode","fixed")
RANDOMER_5=feature("molecular_barcode_rna_adapter","umi","random",group="molecular_barcode",part="RNA adapter")
RANDOMER_10=feature("molecular_barcode_cdna_adapter","umi","random",group="molecular_barcode",part="cDNA adapter")

RIL19 = "AGATCGGAAGAGCGTCGTGTG"
RNA_X1A_BARCODE = "ATATAGG"
AR17 = "ACACGACGCTCTTCCGA"
RAND103TR3_FIXED = "AGATCGGAAGAGCACACGTCTG"
D501_I5 = "TATAGCCT"
D701_I7_OLIGO = "CGAGTAAT"

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("D501 i5", D501_I5, "cbc", feature=I5_FEATURE),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("X1A inline RNA barcode", RNA_X1A_BARCODE, "cbc", feature=INLINE_FEATURE),
    seg("X1A five-base randomer", "U" * 5, "umi", placeholder=True, feature=RANDOMER_5),
    seg("RBP-bound RNA cDNA", "X" * 34, placeholder=True),
    seg("rand103Tr3 ten-base randomer", "V" * 10, "umi", placeholder=True, feature=RANDOMER_10),
    seg("Read 2 arm first base", "A", "r2"),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("D701 i7 reverse complement", "ATTACTCG", "cbc", feature=I7_FEATURE),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="paired-end dual-index eCLIP library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS):
    raise ValueError("eCLIP final library: " + "; ".join(problems))

FINAL_CAPTION = "Original paired-end eCLIP amplicon using the SOP's X1A/D501/D701 examples. Read 1 starts with the seven-base inline barcode and five-base randomer; Read 2 starts with the ten-base cDNA-adapter randomer."
SEQUENCING_INTRO = "Original eCLIP uses paired-end reads and both D5/i5 and D7/i7 index reads. Read 1 begins at the RNA-adapter barcode; Read 2 begins at rand103Tr3's ten-base randomer and approaches the RNA 5′ end/crosslink-stop boundary."

def sections():
    return [
        ("Crosslink, trim and immunoprecipitate", [Row(chunks=[("RNA ————— ×RBP× ————— RNA", None, False)]), Row(chunks=[("UV crosslink   + partial RNase   + antibody beads", "umi", False)]), Row(chunks=[("                 retained RBP–RNA fragment", None, False)])], "UV fixes direct contacts; limited RNase and RBP immunoprecipitation define the captured RNA fragment."),
        ("Ligate the 3′ RNA adapter on beads", adapter_ligation_scene(fragment_name="RBP-bound RNA fragment").rows(), "The indexed RNA adapter is ligated before gel purification and reverse transcription."),
        ("Reverse transcription can stop at the crosslink", [Row(chunks=[("RNA  5′ ————— × peptide ————— adapter — 3′", None, False)]), Row(chunks=[("cDNA 3′ ————— ← RT stop", "r1", False)]), Row(chunks=[("or   3′ ——————————————————— ← read-through cDNA", "r2", False)])], "Both truncated and read-through cDNAs are retained."),
        ("Ligate rand103Tr3 to cDNA", [Row(chunks=[("5′-p NNNNNNNNNN—AGATCGGAAGAGCACACGTCTG ** cDNA ending at crosslink", "umi", False)]), Row(chunks=[("5′-p NNNNNNNNNN—AGATCGGAAGAGCACACGTCTG ** full-length cDNA — RNA-adapter complement", "umi", False)])], "The cDNA 3′ ligation adds a ten-base randomer and bypasses iCLIP circularization, making both cDNA classes amplifiable."),
    ]
