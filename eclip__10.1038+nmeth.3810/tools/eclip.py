"""Original paired-end eCLIP molecular workflow and exact ENCODE SOP adapters."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Scene, complement_segments, feature
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


def crosslinked_rna_rows():
    rna = [
        seg("RNA 5-prime flank", "X" * 14, placeholder=True),
        seg("crosslinked nucleotide", "X", "umi", placeholder=True),
        seg("RNA 3-prime flank", "X" * 14, placeholder=True),
    ]
    sc = Scene()
    sc.strand("RBP-bound RNA", rna)
    sc.mark("RBP-bound RNA", "crosslinked nucleotide", "UV-crosslinked RBP")
    sc.note("RBP-bound RNA", "partial RNase + antibody-bead enrichment")
    return sc.rows()


def reverse_transcription_rows():
    rna = [
        seg("RNA 5-prime flank", "X" * 14, placeholder=True),
        seg("peptide-remnant crosslink", "X", "umi", placeholder=True),
        seg("RNA 3-prime fragment", "X" * 20, placeholder=True),
        seg("3-prime RNA adapter", "X" * 18, "r2", placeholder=True),
    ]
    sc = Scene()
    sc.strand("RNA", rna)
    truncated = complement_segments(rna[2:])
    read_through = complement_segments(rna)
    sc.anneal("crosslink-truncated cDNA", truncated, to="RNA",
              pair=("3-prime RNA adapter'", "3-prime RNA adapter"))
    sc.arrow("crosslink-truncated cDNA", "RT stops at peptide remnant")
    sc.anneal("read-through cDNA", read_through, to="RNA",
              pair=("3-prime RNA adapter'", "3-prime RNA adapter"))
    sc.arrow("read-through cDNA", "read-through RT")
    sc.mark("RNA", "peptide-remnant crosslink", "crosslink site")
    return sc.rows()


def cdna_adapter_ligation_rows():
    def product(name, cdna_name, *, read_through=False):
        parts = [
            seg("rand103Tr3 ten-base randomer", "N" * 10, "umi", placeholder=True,
                feature=RANDOMER_10),
            seg("rand103Tr3 fixed region", RAND103TR3_FIXED, "r2"),
            seg(cdna_name, "X" * 28, placeholder=True),
        ]
        if read_through:
            parts.append(seg("RNA-adapter complement", "X" * 18, "r2", placeholder=True))
        sc = Scene()
        sc.strand(name, parts, mod5="p")
        sc.junction(name, "rand103Tr3 fixed region", cdna_name, "cDNA ligation")
        return sc.rows()

    return [
        *product("crosslink-truncated product", "cDNA ending at crosslink"),
        *product("read-through product", "full-length cDNA", read_through=True),
    ]


def sections():
    return [
        ("Crosslink, trim and immunoprecipitate", crosslinked_rna_rows(), "UV fixes direct contacts; limited RNase and RBP immunoprecipitation define the captured RNA fragment."),
        ("Ligate the 3′ RNA adapter on beads", adapter_ligation_scene(fragment_name="RBP-bound RNA fragment").rows(), "The indexed RNA adapter is ligated before gel purification and reverse transcription."),
        ("Reverse transcription can stop at the crosslink", reverse_transcription_rows(), "Both truncated and read-through cDNAs are retained."),
        ("Ligate rand103Tr3 to cDNA", cdna_adapter_ligation_rows(), "The cDNA 3′ ligation adds a ten-base randomer and bypasses iCLIP circularization, making both cDNA classes amplifiable."),
    ]
