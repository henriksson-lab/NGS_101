"""SaferSeqS copied-duplex-UID and strand-specific target library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import duplex_fragment, seg
from chemdraw import Construct, MolecularState, Segment, Workflow, feature, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from molecular_consensus import CopiedDuplexUID, saferseqs_graph

TITLE = "SaferSeqS — copied duplex UID with strand-specific enrichment"
NOTES = "01_saferseqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-021-00900-z">'
          'Cohen et al., <i>Nature Biotechnology</i> (2021)</a>.')
SUMMARY = ("A 3′-adapter UID is copied into both strands of each original duplex before "
           "amplification. Separate Watson- and Crick-specific nested PCRs preserve strand "
           "identity, allowing paired family consensuses from the same molecule.")
CAVEAT = ("Supplementary Table 7 publishes the adapters, PCR shells and custom sequencing "
          "primers. Target-specific GSPs and eight-base sample indices vary; the exact "
          "library below instantiates one published TP53 GSP2 and a representative index.")

UID = feature("saferseqs_duplex_uid", "umi", "random", group="duplex_uid", part="shared")
I5 = feature("saferseqs_i5", "sample_index", "unknown")
I7 = feature("saferseqs_i7", "sample_index", "unknown")
REPRESENTATIVE_INDEX = "ACGTACGT"
TP53_342_GSP2 = "CTTGAGTTCCAAGGCCTCAT"
CRICK_TARGET_HANDLE = "CACACAGGAAACAGCTATGACCATG"
CRICK_ANCHOR_PRIMER = il.P5 + REPRESENTATIVE_INDEX + il.TRUSEQ_READ1
CRICK_GSP2_PRIMER = il.P7 + REPRESENTATIVE_INDEX + CRICK_TARGET_HANDLE + TP53_342_GSP2
CUSTOM_INDEX_PRIMER = "CATGGTCATAGCTGTTTCCTGTGTG"
CUSTOM_READ2_PRIMER = "GGCATACGAGAT" + REPRESENTATIVE_INDEX + CRICK_TARGET_HANDLE
FINAL_LIBRARY = Construct([
    Segment("P5 graft", il.P5, "p5"),
    Segment("i5 sample index", REPRESENTATIVE_INDEX, "cbc", feature=I5),
    Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
    seg("14-nt copied duplex UID", "N" * 14, "umi", placeholder=True, feature=UID),
    seg("strand-specific target amplicon", "X" * 42, placeholder=True),
    Segment("TP53 p.R342X GSP2 site", revcomp(TP53_342_GSP2), "r2"),
    Segment("Crick target-primer handle", revcomp(CRICK_TARGET_HANDLE), "r2"),
    Segment("i7 sample index reverse complement", revcomp(REPRESENTATIVE_INDEX), "cbc",
            feature=I7),
    Segment("P7 graft reverse complement", il.P7_RC, "p7"),
], name="SaferSeqS Watson-read-orientation representative library")
SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "SaferSeqS custom index primer",
              CUSTOM_INDEX_PRIMER, "Cohen et al. Supplementary Table 7"),
    sp.TRUSEQ["I2"],
    sp.custom("Read 2", "SaferSeqS custom Read 2 primer",
              CUSTOM_READ2_PRIMER, "Cohen et al. Supplementary Table 7"),
)
ACQUIRED_ROLES = sp.ROLES
READ_USES = {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "include"}
FINAL_CAPTION = ("Representative Watson-read-orientation endpoint. Standard Read 1 sees the "
                 "copied 14-nt UID first; the published custom Read 2 primer enters from "
                 "the target-specific end. The opposite Watson orientation swaps which "
                 "paired read encounters the UID first.")
SEQUENCING_INTRO = ("Supplementary Table 7 says the custom read, index and Read 2 primers "
                    "are spiked into the standard Illumina primers. This representative "
                    "Watson-read-orientation molecule uses the standard Read 1 and i5 primers "
                    "and the published custom i7 and Read 2 primers.")
SEQUENCING_UNAVAILABLE = ""


def workflow():
    chemistry = CopiedDuplexUID(14)
    w = Workflow(MolecularState("End-repaired genomic fragment",
                                tuple(duplex_fragment("blunt dephosphorylated DNA").rows())))
    uid_end = AdapterLigationGeometry.unavailable(
        "Supplementary Table 7 publishes both adapter oligos, but their two-ended placement "
        "on the genomic duplex is not yet encoded as an exact ligation geometry.")
    w.ligate_adapter("Ligate the UID-bearing 3′ adapter", uid_end,
            chemistry.three_prime_adapter().rows(), name="3′-adapter-ligated duplex",
            note="The adapter is installed at the two 3′ fragment ends after dephosphorylation and blunting.")
    five_prime_end = AdapterLigationGeometry.unavailable(
        "Supplementary Table 7 publishes the 5′ adapter, but its two-ended placement on the "
        "genomic duplex is not yet encoded as an exact ligation geometry.")
    w.ligate_adapter("Extend the UID and ligate the 5′ adapter", five_prime_end,
            chemistry.copied_barcode().rows(), name="duplex-barcoded library",
            note="Extension copies one random UID into both strands before redundant amplification.")
    w.react("Partition into Watson- and Crick-specific hemi-nested PCRs",
            chemistry.strand_families().rows(), name="oriented target families",
            note="A second nested PCR adds sample indices and further target specificity.")
    return w


def oligos():
    return (
        "3′ N14 Adapter Oligo #1: 5′-/5Phos/ATAAAACGACGGCNNNNNNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAG*T*C/3SpC3/-3′",
        "3′ Adapter Oligo #2: 5′-GCCGUCGUUUUA/33dT/-3′",
        f"5′ Adapter / library amplification primer #1: 5′-{il.TRUSEQ_READ1}-3′",
        "Library amplification primer #2: 5′-GACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′",
        "Watson GSP1 anchor primer: 5′-GACTGGAGTTCAGACGT*G*T-3′",
        "Crick GSP1 anchor primer: 5′-ACACTCTTTCCCTACACG*A*C-3′",
        "Watson GSP2 primer shell: 5′-AATGATACGGCGACCACCGAGATCTACACXXXXXXXXCGACGTAAAACGACGGCCAGT[GSP2]-3′",
        "Crick GSP2 primer shell: 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXCACACAGGAAACAGCTATGACCATG[GSP2]-3′",
        "Watson GSP2 anchor/indexing primer: 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT*C*T-3′",
        "Crick GSP2 anchor/indexing primer: 5′-AATGATACGGCGACCACCGAGATCTACACXXXXXXXXACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T-3′",
        "Custom Read 1 sequencing primer: 5′-GAGATCTACACNNNNNNNNCGACGTAAAACGACGGCCA*G*T-3′",
        "Custom index sequencing primer: 5′-CATGGTCATAGCTGTTTCCTGTG*T*G-3′",
        "Custom Read 2 sequencing primer: 5′-GGCATACGAGATNNNNNNNNCACACAGGAAACAGCTATGACCA*T*G-3′",
    )


protocol_graph = lambda: saferseqs_graph(
    SOURCE, FINAL_LIBRARY,
    "Target-specific and indexed nested-PCR oligos vary by assay; the page models one "
    "published TP53 GSP2 with a representative eight-base index.",
    primers=SEQ_PRIMERS, read_uses=READ_USES)
