"""TRIP reporter-barcode counting and integration-site mapping chemistry."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      circle_rows, feature, oligo, revcomp)
from pcr import PCRGeometry, PCRPrimer

TITLE = "TRIP — Thousands of Reporters Integrated in Parallel"
NOTES = "01_trip.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nprot.2014.072">'
          'Akhtar et al., <i>Nature Protocols</i> (2014)</a>; exact indexed primers and '
          'library structures are in Supplementary Table 1 and Supplementary Figure 3.')
SUMMARY = ("A random 16-nt barcode identifies each piggyBac reporter. Barcode counts from "
           "reporter cDNA measure expression, counts from genomic DNA normalize reporter "
           "abundance, and DpnII circularization plus inverse PCR maps each barcode to its "
           "genomic integration site.")
CAVEAT = ("The reporter and primer sequences below are published. Genomic flanks are "
          "variable biological sequence, represented as placeholders rather than bases.")

# Nature Protocols 9, 1760–1769 (2014), reagent list and Supplementary Table 1.
BARCODE_LONG = "GTGACACCTGCAGGATCA" + "N" * 16 + "CTCGAGTTGTGGCCGGCCCTTGTGACTG"
BSSH2_A = "AAAAGCGCGCATACTAGATTAACCCTAGAAAGATAATCATATTG"
INDEX_EXAMPLE = "ACGAGTGCGT"
PB_CDNA_F = il.TRUSEQ_READ1 + INDEX_EXAMPLE + "GTCACAAGGGCCGGCCACAA"
PB_CDNA_R = "CAAGCAGAAGACGGCATACGAGATGATCACATGGTCCTGCTGGAGTT"
PB_OUTER_F = "TTTTACGCATGATTATCTTTAACGTACGTC"
CDNA_AMPL_R = "CGCCAGGGTTTTCCCAGTCACAAG"
INV_PCR_F_TRUSEQ_PE2 = ("GAT" + il.TRUSEQ_READ2
                         + "GTACGTCACAATATGATTATCTTTCTAG")
SOLEXA_R_TRUSEQ_PE2 = il.P7 + il.TRUSEQ_READ2[:23]
EGFP_INT3 = "TCAAGGTGAACTTCAAGATCCG"

LEFT_SITE = "GTCACAAGGGCCGGCCACAA"
REPORTER_LEFT = "CTCGAG"
REPORTER_RIGHT = "TGATCCTGCAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTT"
CDNA_R_SITE = "GATCACATGGTCCTGCTGGAGTT"
PB_3TR_TOP = revcomp("GTACGTCACAATATGATTATCTTTCTAG")


def seg(name: str, bases: str, tag: str | None = None, **kwargs) -> Segment:
    return Segment(name, bases, tag, **kwargs)


def reporter_barcode(name: str = "16-nt reporter barcode") -> Segment:
    return seg(name, "N" * 16, "cbc", placeholder=True,
               feature=feature("reporter_barcode", "feature_barcode", "random"))


def sample_index() -> Segment:
    return seg("10-nt inline sample index", INDEX_EXAMPLE, "cbc",
               feature=feature("sample_index", "sample_index", "fixed"))


def p5_read1_parts() -> list[Segment]:
    # P5 and Read 1 share four ACAC bases in the physical oligo; splitting at P5 avoids
    # duplicating that overlap while seqprimers still finds the complete Read 1 site.
    return [seg("P5", il.P5, "p5"),
            seg("remainder of Read 1 site", il.TRUSEQ_P5_FULL[len(il.P5):], "r1")]


def barcode_payload(*, mapping: bool = False) -> list[Segment]:
    parts = [seg("reporter left constant", LEFT_SITE + REPORTER_LEFT, "me"),
             reporter_barcode()]
    if mapping:
        parts.extend([seg("reporter sequence before DpnII site", "T", "me"),
                      seg("restored DpnII junction", "GATC", "me"),
                      seg("genomic integration flank", "X" * 42, placeholder=True),
                      seg("piggyBac 3-prime terminal repeat", PB_3TR_TOP, "w1")])
    else:
        parts.append(seg("reporter right constant", REPORTER_RIGHT, "me"))
    return parts


def count_template(name: str) -> Construct:
    return Construct([seg("PB-cDNA-forward annealing site", LEFT_SITE, "r1"),
                      seg("reporter left constant", REPORTER_LEFT, "me"),
                      reporter_barcode(),
                      seg("reporter right constant", REPORTER_RIGHT, "me"),
                      seg("PB-cDNA-reverse site", revcomp(CDNA_R_SITE), "r2")], name=name)


def count_pcr_geometry(name: str) -> PCRGeometry:
    return PCRGeometry.exact(
        count_template(name), PCRPrimer("PB-cDNA-forward-1-BC", PB_CDNA_F),
        PCRPrimer("PB-cDNA-reverse", PB_CDNA_R),
        statement="The protocol prints both primers; their 3-prime regions select the reporter barcode cassette while their 5-prime tails install Read 1/index/P5-side and P7-side sequence.")


def count_library(name: str) -> Construct:
    con = Construct([*p5_read1_parts(), sample_index(), *barcode_payload(),
                     seg("PB-cDNA-reverse product", revcomp(PB_CDNA_R), "p7")], name=name)
    primers = sp.standard_primers(con, ("Read 1",))
    problems = sp.verify(con, primers, required_roles=("Read 1",))
    if problems:
        raise ValueError(f"invalid {name}: " + "; ".join(problems))
    return con


def mapping_template() -> Construct:
    return Construct([seg("PB-cDNA-forward annealing site", LEFT_SITE, "r1"),
                      seg("reporter left constant", REPORTER_LEFT, "me"),
                      reporter_barcode(),
                      seg("reporter sequence before DpnII site", "T", "me"),
                      seg("DpnII ligation junction", "GATC", "me"),
                      seg("genomic integration flank", "X" * 42, placeholder=True),
                      seg("piggyBac 3-prime terminal repeat", PB_3TR_TOP, "w1")],
                     name="opened self-ligated DpnII mapping template")


def mapping_pcr1_geometry() -> PCRGeometry:
    return PCRGeometry.exact(
        mapping_template(), PCRPrimer("PB-cDNA-forward-1-BC", PB_CDNA_F),
        PCRPrimer("InvPCR-F-TruSeqPE2", INV_PCR_F_TRUSEQ_PE2),
        statement="The indexed reporter primer faces the published piggyBac-terminal-repeat inverse-PCR primer across the DpnII-linked genomic flank.")


def mapping_pcr1_product() -> Construct:
    return Construct([*p5_read1_parts(), sample_index(), *barcode_payload(mapping=True),
                      seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "r2")],
                     name="TRIP mapping PCR1 product")


def mapping_pcr2_geometry() -> PCRGeometry:
    return PCRGeometry.exact(
        mapping_pcr1_product(), PCRPrimer("Solexa-Adapter-SE1", il.TRUSEQ_P5_FULL),
        PCRPrimer("Solexa-R-Ad-TruSeqPE2", SOLEXA_R_TRUSEQ_PE2),
        statement="PCR2 uses the printed universal P5/Read 1 primer and P7/Read 2 primer to complete the flow-cell arms around the inverse-PCR product.")


def mapping_library() -> Construct:
    con = Construct([*p5_read1_parts(), sample_index(), *barcode_payload(mapping=True),
                     seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "r2"),
                     seg("P7", il.P7_RC, "p7")], name="TRIP integration-mapping library")
    primers = sp.standard_primers(con, ("Read 1", "Read 2"))
    problems = sp.verify(con, primers, required_roles=("Read 1", "Read 2"))
    if problems:
        raise ValueError("invalid mapping library: " + "; ".join(problems))
    return con


EXPRESSION = count_library("TRIP expression barcode library")
NORMALIZATION = count_library("TRIP genomic normalization barcode library")
MAPPING = mapping_library()
COUNT_PRIMERS = sp.standard_primers(EXPRESSION, ("Read 1",))
MAPPING_PRIMERS = sp.standard_primers(MAPPING, ("Read 1", "Read 2"))

FINAL_LIBRARIES = (
    ("Expression library", EXPRESSION, COUNT_PRIMERS,
     "Reporter cDNA barcode amplicon; Read 1 counts the inline sample index and reporter barcode.",
     "The expression and normalization products have the same sequenced construct but originate from cDNA and genomic DNA, respectively."),
    ("Genomic normalization library", NORMALIZATION, COUNT_PRIMERS,
     "Reporter genomic-DNA barcode amplicon used to normalize clone abundance.",
     "Only Read 1 was acquired for barcode counting; the ten-base sample index is inline rather than an index read."),
    ("Integration-site mapping library", MAPPING, MAPPING_PRIMERS,
     "Paired-end inverse-PCR product linking one reporter barcode to its genomic integration flank.",
     "Read 1 identifies the sample and reporter barcode; Read 2 crosses the piggyBac terminal repeat into variable genomic DNA."),
)

READ_LENGTHS = {
    "Expression library": {"Read 1": 100},
    "Genomic normalization library": {"Read 1": 100},
    "Integration-site mapping library": {"Read 1": 100, "Read 2": 100},
}
READ_USES = {
    "Expression library": {"Read 1": "extract"},
    "Genomic normalization library": {"Read 1": "extract"},
    "Integration-site mapping library": {"Read 1": "extract", "Read 2": "include"},
}


def oligos():
    return [
        oligo("PB-barcode-long-7", [seg("left cloning arm", "GTGACACCTGCAGGATCA"),
              reporter_barcode(), seg("right cloning arm", "CTCGAGTTGTGGCCGGCCCTTGTGACTG")]),
        oligo("BssH2-A", [seg("BssH2-A", BSSH2_A)]),
        oligo("PB-cDNA-forward-1-BC (index 1 example)",
              [seg("Read 1", il.TRUSEQ_READ1, "r1"), sample_index(),
               seg("reporter annealing site", LEFT_SITE)]),
        oligo("PB-cDNA-reverse", [seg("P7", il.P7, "p7"),
              seg("reporter annealing region", CDNA_R_SITE)]),
        oligo("Solexa-Adapter-SE1", [seg("P5 + Read 1", il.TRUSEQ_P5_FULL, "p5")]),
        oligo("InvPCR-F-TruSeqPE2", [seg("5-prime GAT", "GAT"),
              seg("Read 2", il.TRUSEQ_READ2, "r2"),
              seg("piggyBac 3-prime-repeat annealing region", INV_PCR_F_TRUSEQ_PE2[-28:])]),
        oligo("Solexa-R-Ad-TruSeqPE2", [seg("P7", il.P7, "p7"),
              seg("Read 2 prefix", il.TRUSEQ_READ2[:23], "r2")]),
        oligo("eGFP-int3 sequencing primer", [seg("eGFP-int3", EGFP_INT3)]),
    ]


def _reporter_scene() -> Scene:
    sc = Scene.duplex([seg("piggyBac left arm", "X" * 22, placeholder=True),
                       seg("reporter cassette", "X" * 28, placeholder=True),
                       reporter_barcode(),
                       seg("piggyBac right arm", "X" * 22, placeholder=True)],
                      label="barcoded piggyBac reporter")
    sc.mark("top", "16-nt reporter barcode", "one random barcode identifies one integrated reporter")
    sc.labels("top")
    return sc


def _mapping_circle_rows() -> list[Row]:
    circle = Construct([seg("barcoded reporter fragment", "X" * 34, placeholder=True),
                        seg("DpnII junction", "GATC", "me"),
                        seg("genomic integration flank", "X" * 42, placeholder=True),
                        seg("piggyBac 3-prime terminal repeat", PB_3TR_TOP, "w1")],
                       name="self-ligated DpnII fragment")
    return circle_rows(circle, "dilute T4-ligase closure")


def workflow() -> Workflow:
    w = Workflow(MolecularState("Randomly barcoded reporter pool", tuple(_reporter_scene().rows())))
    w.react("Integrate reporters and split DNA/RNA aliquots",
            [Row(chunks=[("expression: poly(A) RNA -> oligo(dT) cDNA", None, False)]),
             Row(chunks=[("normalization: genomic reporter DNA", None, False)]),
             Row(chunks=[("mapping: genomic reporter plus insertion flank", None, False)])],
            name="three matched molecular readouts",
            note="Barcode identity joins expression, reporter abundance and integration position.")
    count_geom = count_pcr_geometry("expression or normalization reporter template")
    w.pcr("PCR1: amplify expression and normalization barcodes", count_geom,
          count_geom.rows(), name="indexed reporter-barcode amplicon",
          note="The 10-nt sample index is introduced inline by the forward primer; PCR2 with Solexa-Adapter-SE1 completes P5.")
    w.react("DpnII-digest mapping DNA and self-ligate at low concentration",
            _mapping_circle_rows(), name="circular mapping template",
            note="DpnII supplies compatible GATC ends; intramolecular ligation places the unknown genomic flank between known reporter and piggyBac sequences.")
    geom1 = mapping_pcr1_geometry()
    w.pcr("Inverse PCR across the reporter–genome junction", geom1, geom1.rows(),
          name="mapping PCR1 product",
          note="The primer pair faces outward on the original integrated reporter but inward after circularization.")
    geom2 = mapping_pcr2_geometry()
    w.pcr("PCR2 completes P5 and P7 flow-cell arms", geom2, geom2.rows(),
          name="paired-end mapping library",
          note="The second printed primer pair completes the standard TruSeq Read 1 and Read 2 sites.")
    return w
