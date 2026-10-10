"""Molecular model for LAST-seq (Lyu and Chen, 2023)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, feature, oligo)
from rna_amplification import direct_ssrna_ivt, poly_a_capture


def seg(name, bases, tag=None, **kw):
    return Segment(name, bases, tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
RANDOM_RT_HANDLE = "TACACGACGCTCTTCCGATCT"
RANDOM_RT = RANDOM_RT_HANDLE + "NNNNNN"
DTRU_LINKER = "CCACCTTTCATTCACCCT"
DTRU_DNA = "TTTTTT"
DTRU_RNA = "U" * 19
DTRU_TAIL = DTRU_LINKER + DTRU_DNA + DTRU_RNA
CELL_BARCODE_LEN = 6
UMI_LEN = 8

CELL_FEATURE = feature("last_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="Supplementary Table S3")
UMI_FEATURE = feature("last_umi", "umi", "random")


def capture_primer():
    return [
        seg("hairpin module", "[hairpin: PCR handle + T7 promoter + UMI + cell barcode]",
            "r2", placeholder=True),
        seg("complementary linker", DTRU_LINKER, "tso"),
        seg("dT7", DTRU_DNA),
        seg("rU19", DTRU_RNA),
    ]


def capture_scene():
    sc = poly_a_capture(capture_primer(), annealing_segment="rU19",
                        primer_label="3′-ddC-blocked LAST-seq primer")
    sc.junction("primer", "hairpin module", "complementary linker", "primer assembly ligation")
    sc.note("primer", "the dT/rU tract is nicked by RNase H; 3′ ddC blocks extension from the primer end")
    return sc


def promoter_scene():
    promoter = seg("T7 promoter", T7_PROMOTER, "t7")
    template = [seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
                    feature=CELL_FEATURE),
                seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
                seg("mRNA template", "X" * 34, placeholder=True)]
    return direct_ssrna_ivt(template, promoter)


def arna_scene():
    sc = Scene(); sc.strand("aRNA", [
        seg("aRNA transcript", "X" * 34, placeholder=True),
        seg("copied UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
        seg("copied cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
            feature=CELL_FEATURE)], label="antisense RNA copies")
    sc.note("aRNA", "hundreds of linear-amplification products from one original ssRNA")
    return sc


def random_rt_scene():
    arna = [seg("aRNA body", "X" * 34, placeholder=True)]
    primer = [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
              seg("random N6", "N" * 6, placeholder=True)]
    sc = Scene(); sc.strand("aRNA", arna, label="amplified antisense RNA")
    sc.anneal("random RT primer", primer, to="aRNA", pair=("random N6", "aRNA body"),
              shift=len(arna[0]) - 6, label="random RT primer",
              unpaired=("Read 1 handle", "random N6", "aRNA body"))
    sc.arrow("random RT primer", "SuperScript IV")
    sc.labels("random RT primer")
    return sc


FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("random-primer N6", "N" * 6, placeholder=True),
    seg("aRNA-derived cDNA", "X" * 38, placeholder=True),
    seg("cell barcode", "B" * CELL_BARCODE_LEN, "cbc", placeholder=True,
        feature=CELL_FEATURE),
    seg("UMI", "U" * UMI_LEN, "umi", placeholder=True, feature=UMI_FEATURE),
    seg("Read 2 arm completion base", "A", "r2"),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
    seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
        feature=feature("last_sample_index", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="LAST-seq library")

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_errors = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                    required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _errors:
    raise ValueError("LAST-seq final library: " + "; ".join(_errors))


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="LAST-seq")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Poly(A) capture by the assembled hairpin/rU-dT primer",
                                 tuple(capture_scene().rows())))
    wf.react("RNase H nicking and Klenow exo− extension", promoter_scene().rows(),
             note="A short rA/dT hybrid is nicked and extended so the original RNA gains a duplex T7 promoter without reverse transcription.")
    wf.react("T7 in-vitro transcription directly from ssRNA", arna_scene().rows(),
             note="The original single-stranded RNA is the transcription template; amplification is linear.")
    wf.react("Random-prime the pooled aRNA and reverse-transcribe", random_rt_scene().rows(),
             note="Reverse transcription occurs only after direct RNA amplification.")
    wf.react("P5/P7 PCR and library enrichment", final_rows(),
             note="PCR completes the Illumina library; the inline cell barcode and UMI originate in the capture hairpin.")
    return wf


TITLE = "LAST-seq — direct linear amplification of single-stranded RNA"
NOTES = "01_last-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1186/s13059-023-03025-5">'
          'Lyu and Chen, <i>Genome Biology</i> (2023)</a>.')
SUMMARY = ("A ligated hairpin/rU-dT primer installs a short duplex T7 promoter at an mRNA’s "
           "3′ end. T7 polymerase then amplifies the original RNA directly, before any reverse transcription.")
FINAL_CAPTION = "Representative PCR-completed LAST-seq molecule. Read 2 reaches the hairpin-derived UMI and cell barcode; Read 1 enters the randomly primed aRNA-derived insert."
SEQUENCING_INTRO = "NextSeq 550: Read 1, 50 cycles; Read 2, 25 cycles. The paper reports 6-nt cellular barcodes and 8-nt UMIs."
READ_LENGTHS = {"Read 1": 50, "Read 2": 25}


def oligos():
    yield oligo("rU-dT module", [seg("linker", DTRU_LINKER, "tso"),
                                  seg("dT7", DTRU_DNA), seg("rU19", DTRU_RNA)],
                mods="5′ phosphate; 3′ ddC")
    yield oligo("RandomRT primer6", [seg("Read 1 handle", RANDOM_RT_HANDLE, "r1"),
                                      seg("N6", "N" * 6, placeholder=True)])
