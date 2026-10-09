"""Peterson et al. ddRAD-seq library chemistry."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Row, Scene, Segment
from restriction import ECORI, MSPI

TITLE = "ddRAD-seq — double-digest RAD sequencing"
NOTES = "01_ddrad-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1371/journal.pone.0037135">'
          'Peterson et al., <i>PLOS ONE</i> (2012)</a>.')
SUMMARY = ("Digest genomic DNA with two enzymes, ligate end-selective adapters carrying "
           "an inline sample barcode, pool and size-select, then add the Illumina outer "
           "adapter and a second sample index by PCR.")
CAVEAT = "The displayed GCATG inline barcode and ATCACG i7 are examples from the published oligo table."

INLINE = "GCATG"
I7 = "ATCACG"
P1_TOP = il.TRUSEQ_READ1 + INLINE
P1_BOTTOM = ECORI.overhang + "CATGC" + il.INDEX2_PRIMER_RC
P2_TOP = il.TRUSEQ_READ2
P2_BOTTOM = MSPI.overhang + "AGATCGGAAGAGCGAGAACAA"

# A source transcription is accepted only through the enzyme-derived end model.
P1_END = ECORI.adapter_end(P1_BOTTOM[:len(ECORI.overhang)])
P2_END = MSPI.adapter_end(P2_BOTTOM[:len(MSPI.overhang)])


def final_library() -> tuple[Construct, tuple]:
    lib = Construct([
        Segment("P5 before shared ACAC", il.P5[:-4], "p5"),
        Segment("Read 1 site (shares ACAC with P5)", il.TRUSEQ_READ1, "r1"),
        Segment("inline sample barcode", INLINE, "cbc"),
        Segment("EcoRI remnant", ECORI.overhang),
        Segment("size-selected genomic insert", "X" * 34, placeholder=True),
        Segment("MspI remnant", MSPI.overhang),
        Segment("Read 2 site start", "A", "r2"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("6-nt i7", I7, "cbc"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="ddRAD-seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid ddRAD-seq endpoint: " + "; ".join(problems))
    return lib, primers


FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("A P1–insert–P2 molecule after PCR. Read 1 encounters the inline barcode "
                 "before the EcoRI-associated genomic end; Index 1 reports i7, and Read 2 "
                 "enters from the MspI-associated end.")
SEQUENCING_INTRO = ("Standard TruSeq Read 1, Index 1 and Read 2 primers bind. Demultiplexing "
                    "uses the inline P1 barcode together with the optional six-base i7.")


def sections():
    digest = [
        Row(chunks=[(f"EcoRI  G^{ECORI.overhang}C  →  5′-{ECORI.overhang} cohesive end", None, False)]),
        Row(chunks=[(f"MspI   C^{MSPI.overhang}G    →  5′-{MSPI.overhang} cohesive end", None, False)]),
    ]
    lib = FINAL_LIBRARY
    sc = Scene.duplex(list(lib), label="selected ddRAD molecule")
    sc.junction("top", "inline sample barcode", "EcoRI remnant", "P1 ligation")
    sc.junction("top", "MspI remnant", "Read 2 site start", "P2 ligation")
    sc.labels("top")
    return [
        ("Double digest", digest,
         "EcoRI-HF and MspI make distinguishable 5′ cohesive ends; their overhangs are derived from the cut sites."),
        ("Ligate end-selective adapters", sc.rows(),
         "P1 carries the inline barcode and fits the EcoRI end; P2 fits the MspI end. ** marks each ligation."),
        ("Pool, size-select and PCR", [Row(chunks=[
            ("barcoded ligation products → pooled size window → 8–12 PCR cycles → indexed library", None, False)])],
         "Pooling occurs after inline barcoding; a narrow size window defines the reduced genomic representation."),
    ]
