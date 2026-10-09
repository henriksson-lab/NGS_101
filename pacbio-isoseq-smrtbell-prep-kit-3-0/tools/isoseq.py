"""PacBio Iso-Seq v2 / SMRTbell prep kit 3.0 workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, feature, feature_rows, strand_row
from rna_special import isoseq_smrtbell, template_switch_scene

TITLE = "PacBio Iso-Seq v2 — SMRTbell prep kit 3.0"
NOTES = "01_isoseq.html"
SOURCE = 'Commercial protocol: PacBio <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-Iso-Seq-libraries-using-SMRTbell-prep-kit-3.0.pdf">Preparing Iso-Seq v2 libraries using SMRTbell prep kit 3.0, 102-396-000 Rev07</a>.'
SUMMARY = "Template-switch reverse transcription captures full-length poly(A) transcripts, cDNA is amplified, repaired and A-tailed, then hairpin adapters close both ends into a SMRTbell for circular-consensus sequencing."
CAVEAT = "Current cDNA-synthesis oligos and SMRTbell adapter bases are proprietary. Bracketed/unknown regions preserve molecular roles and topology without inventing sequences."
INSERT = Construct([
    seg("Iso-Seq primer barcode", "[Iso-Seq barcode]", "cbc", placeholder=True,
        inferred=True,
        feature=feature("isoseq_sample_index", "sample_index", "whitelist",
                        whitelist="PacBio Iso-Seq primer barcodes bc01–12")),
    seg("5-prime cDNA handle", "X"*14, "r1", placeholder=True, inferred=True),
    seg("full-length transcript cDNA", "X"*44, placeholder=True),
    seg("poly(A/T) boundary", "X"*12, placeholder=True),
    seg("3-prime cDNA handle", "X"*14, "r2", placeholder=True, inferred=True),
], name="barcoded full-length Iso-Seq cDNA")
SMRTBELL, SMRTBELL_ROWS = isoseq_smrtbell(INSERT)
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "The SMRTbell sequencing primer anneals within the proprietary hairpin adapter. A bound polymerase repeatedly traverses forward cDNA, the opposite hairpin, reverse-complement cDNA and the starting hairpin; those passes form one circular-consensus (HiFi) read."

def sections():
    return [
        ("Make full-length cDNA by template switching", template_switch_scene(inferred=True).rows(), "INFERRED — the vendor discloses oligo roles but not current bases. Poly(A)-primed RT and template switching place PCR handles at both cDNA ends."),
        ("Amplify, barcode and repair full-length cDNA",
         [strand_row(INSERT), *feature_rows(INSERT, prefix_width=4),
          Row(chunks=[("PCR → damage repair → end repair → dA tailing", None, False)])],
         "The required bc01–12 primer installs the sample index during cDNA PCR. Repair and A-tailing make both duplex ends adapter-compatible."),
        ("Ligate hairpins to make the SMRTbell", SMRTBELL_ROWS, "Each ** boundary is a hairpin-adapter ligation. The Dumbbell constructor rejects missing ends or incompatible A:T overhangs."),
        ("Polymerase repeatedly traverses the insert", [Row(chunks=[("primer → insert forward → hairpin → insert reverse complement → hairpin ↻", "r1", False)])], "Multiple passes over the same full-length cDNA molecule yield the circular-consensus sequence."),
    ]
