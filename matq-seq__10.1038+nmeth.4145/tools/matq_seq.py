"""Molecular model for MATQ-seq (Sheng et al., 2017)."""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Workflow, annotation_rows, oligo
from rna_amplification import poly_a_capture, terminal_tail

GAT27 = "GTGAGTGATGGTTGAGGATGTGTGGAG"
GAT27_DT = GAT27 + "N" * 5 + "T" * 20
GAT27_5N3G = GAT27 + "N" * 5 + "GGG"
GAT27_5N3T = GAT27 + "N" * 5 + "TTT"
GAT21_6N3G = "GATGGTTGAGGATGTGTGGAG" + "N" * 6 + "GGG"
THREE_N_GAT24 = "N" * 3 + "AGTGATGGTTGAGGATGTGTGGAG"


def initial_scene():
    primer = [seg("GAT27 handle", GAT27, "tso"),
              seg("five random bases", "N" * 5, placeholder=True),
              seg("dT20", "T" * 20)]
    sc = poly_a_capture(primer, annealing_segment="dT20",
                        primer_label="GAT27dT primer")
    sc.note("primer", "the same RT also contains GAT27-5N3G and GAT27-5N3T random primers")
    return sc


FIRST_STRAND = Construct([
    seg("GAT27 handle", GAT27, "tso"), seg("primer N5", "N" * 5, placeholder=True),
    seg("first-strand cDNA", "X" * 38, placeholder=True)], name="MATQ first strand")


def first_strand_scene():
    sc = Scene(); sc.strand("cDNA", list(FIRST_STRAND), label="first-strand cDNA")
    sc.arrow("cDNA", "ten stepped-temperature annealing/extension cycles")
    sc.labels("cDNA")
    return sc


TAILED, _TAIL_SCENE = terminal_tail(FIRST_STRAND, "C", 12, enzyme="TdT")


def second_strand_scene():
    primer = [seg("GAT21 handle", "GATGGTTGAGGATGTGTGGAG", "tso"),
              seg("random N6", "N" * 6, placeholder=True), seg("3G", "GGG")]
    sc = Scene(); sc.strand("dC-tailed cDNA", list(TAILED), label="dC-tailed first strand")
    sc.anneal("GAT21-6N3G", primer, to="dC-tailed cDNA",
              pair=("3G", "poly(dC) tail"), label="second-strand primer",
              unpaired=("GAT21 handle", "random N6"))
    sc.arrow("GAT21-6N3G", "Deep Vent exo−; repeated annealing and extension")
    sc.labels("GAT21-6N3G")
    return sc


FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("sheared MATQ cDNA", "X" * 42, placeholder=True)],
    "MATQ-seq paired-end library", inferred_adapters=True)


def final_rows():
    sc = Scene.duplex(list(FINAL_LIBRARY), label="MATQ-seq")
    sc.junction("top", "Read 1 arm", "sheared MATQ cDNA", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return [*sc.rows(), *annotation_rows(FINAL_LIBRARY)]


def workflow():
    wf = Workflow(MolecularState("Mixed oligo-dT and random priming",
                                 tuple(initial_scene().rows())))
    wf.react("SuperScript III reverse transcription", first_strand_scene().rows(),
             note="Ten low-to-high-temperature cycles repeatedly anneal the mixed primers across total RNA.")
    wf.react("T4 DNA polymerase primer digestion, then RNase H/RNase I",
             first_strand_scene().rows(),
             note="Unused primers and the RNA template are removed before homopolymer tailing.")
    wf.react("TdT adds a poly(dC) tail", _TAIL_SCENE.rows(),
             note="The common GAT27-bearing first strand receives a terminal dC tract.")
    wf.react("GAT21-6N3G priming and second-strand synthesis",
             second_strand_scene().rows(),
             note="The terminal GGG anchors at the dC tail; random bases promote repeated annealing/looping amplification.")
    wf.react("GAT27 PCR, 3NGAT24 conversion, shearing and Illumina library construction",
             final_rows(),
             note="The supplementary protocol defines WTA and shearing; the exact sequencing-adapter kit sequence was not printed.")
    return wf


TITLE = "MATQ-seq — multiple annealing and dC-tailing single-cell RNA-seq"
NOTES = "01_matq-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nmeth.4145">'
          'Sheng et al., <i>Nature Methods</i> (2017)</a>.')
SUMMARY = ("A mixture of oligo-dT and random primers repeatedly samples total RNA. After RT, "
           "TdT adds dC, a G-anchored primer makes the second strand, and one common handle amplifies the transcriptome.")
CAVEAT = "INFERRED — the source specifies paired-end Illumina sequencing and shearing but does not print the final adapter kit’s molecular sequence; those terminal regions are shown as an inferred canonical TruSeq geometry."
FINAL_CAPTION = "INFERRED — representative sheared MATQ amplicon in a canonical paired-end Illumina library; the GAT27 amplification products themselves are source-defined."
SEQUENCING_INTRO = "The supplementary figures report paired-end 85-base reads. Standard TruSeq primer sites are inferred from the unreported final adapter kit."
READ_LENGTHS = {"Read 1": 85, "Read 2": 85}


def oligos():
    yield oligo("GAT27dT", [seg("GAT27", GAT27, "tso"),
                             seg("N5", "N" * 5, placeholder=True), seg("dT20", "T" * 20)])
    yield oligo("GAT27-5N3G", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27-5N3T", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3T", "TTT")])
    yield oligo("GAT21-6N3G", [seg("GAT21", "GATGGTTGAGGATGTGTGGAG", "tso"),
                                seg("N6", "N" * 6, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27 PCR", [seg("GAT27", GAT27, "tso")])
    yield oligo("3NGAT24", [seg("N3", "N" * 3, placeholder=True),
                             seg("GAT24", "AGTGATGGTTGAGGATGTGTGGAG", "tso")])
