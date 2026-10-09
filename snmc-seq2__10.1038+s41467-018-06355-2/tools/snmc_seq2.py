"""Luo et al. 2018 snmC-seq2 chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from base_conversion import bisulfite_path
from batch_ngs import seg, truseq_library
from chemdraw import Row, Scene, feature

TITLE = "snmC-seq2 — indexed single-nucleus methylomes"
NOTES = "01_snmc-seq2.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41467-018-06355-2">Luo et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY = "Bisulfite-converted nuclear DNA is copied with one of eight indexed random primers, completed by Adaptase, and amplified with unique dual indexes."
CAVEAT = ("The defining supplement prints the random-primer and dual-index families, but not the bases added by the proprietary Adaptase reaction. "
          "That short region is explicitly inferred below.")
INLINE = "CGATGT"  # representative P5L_AD002 random-primer barcode


def library():
    inserts = [seg("6-nt inline cell barcode", INLINE, "cbc",
                   feature=feature("cell_inline", "cell_barcode", "whitelist", whitelist="published random-primer barcode set")),
               seg("random N9", "N" * 9, placeholder=True),
               seg("bisulfite-converted genomic insert", "X" * 36, placeholder=True),
               seg("Adaptase-added bases", "[Adaptase bases]", placeholder=True, inferred=True,
                   note="proprietary Adaptase sequence is not printed")]
    return truseq_library(inserts, "snmC-seq2 library", dual_index=True)


FINAL_LIBRARY, SEQ_PRIMERS = library()
FINAL_CAPTION = "INFERRED — dual-index endpoint with the source-printed 6-nt inline barcode and N9; the Adaptase-added bases are unavailable."
SEQUENCING_INTRO = "Read 1 begins with the inline cell barcode and random-primer sequence; i5 and i7 form the unique dual index."


def sections():
    return [
        ("Bisulfite-convert each nucleus", [
            Row(chunks=[(bisulfite_path(protected=False).text(), None, False)]),
            Row(chunks=[(bisulfite_path(protected=True).text(), "w1", False)]),
        ], "Unmethylated C is copied as T; 5mC remains C."),
        ("Indexed random priming", [
            Row(chunks=[("5′ spacer — Read 1 tail — 6-nt cell barcode — N9 → extension", "cbc", False)])
        ], "Eight source-printed random primers barcode plate quadrants before pooling."),
        ("INFERRED — Adaptase completion and indexed PCR", [
            Row(chunks=[("random-primer product → Adaptase tail → P5/i5 + P7/i7 PCR", None, False)])
        ], "The PCR primer families are printed; only the proprietary Adaptase-added bases are unknown."),
    ]
