"""Molecular model for Oxford Nanopore SQK-NBD114.24 native barcoding."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from endprep import dA_tailed_scene, repair_and_dA_tail
from nanopore_library import motor_library, no_primer_ending, role, single_strand_rows

TITLE = "Oxford Nanopore Native Barcoding V14 — SQK-NBD114.24"
NOTES = "01_native-barcoding-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24">'
          'NBE_9169_v114</a>.')
SUMMARY = ("PCR-free sample multiplexing: repair and dA-tail each DNA sample, ligate a "
           "dT-tailed native barcode, pool, then ligate the motor-bearing Native Adapter "
           "through the barcode adapter's cohesive end.")
CAVEAT = ("Barcode-adapter flanks, the cohesive hook, Native Adapter bases and motor "
          "attachment are proprietary. Bracketed regions encode disclosed roles only.")

NATIVE_BARCODES = {
    "NB01": "CACAAAGACACCGACAACTTTCTT", "NB02": "ACAGACGACTACAAACGGAATCGA",
    "NB03": "CCTGGTAACTGGGACACAAGACTC", "NB04": "TAGGGAAACACGATAGAATCCGAA",
    "NB05": "AAGGTTACACAAACCCTGGACAAG", "NB06": "GACTACTTTCTGCCTTTGCGAGAA",
    "NB07": "AAGGATTCATTCCCACGGTAACAC", "NB08": "ACGTAACTTGGTTTGTTCCCTGAA",
    "NB09": "AACCAAGACTCGCTGTGCCTAGTT", "NB10": "GAGAGGACAAAGGTTTCAACGCTT",
    "NB11": "TCCATTCCCTCCGATAGATGAAAC", "NB12": "TCCGATTCTGCTTCTTTCTACCTG",
    "NB13": "AGAACGACTTCCATACTCGTGTGA", "NB14": "AACGAGTCTCTTGGGACCCATAGA",
    "NB15": "AGGTCTACCTCGCTAACACCACTG", "NB16": "CGTCAACTGACAGTGGTTCGTACT",
    "NB17": "ACCCTCCAGGAAAGTACCTCTGAT", "NB18": "CCAAACCCAACAACCTAGATAGGC",
    "NB19": "GTTCCTCGTGCAGTGTCAAGAGAT", "NB20": "TTGCGTCCTGTTACGAGAACTCAT",
    "NB21": "GAGCCTCTCATTGTCCGTTCTCTA", "NB22": "ACCACTGCCATGTATCAAAGTACG",
    "NB23": "CTTACTACCCAGTGAACCTCCTCG", "NB24": "GCATAGTTCTGCATGATGGGTTAG",
}
NB01 = NATIVE_BARCODES["NB01"]
INSERT = Construct([Segment("genomic DNA", "X" * 42, placeholder=True)], name="genomic DNA")
END_PREP = repair_and_dA_tail(INSERT)
BARCODE_FEATURE = feature("native_barcode", "sample_index", "whitelist",
                          whitelist="Oxford Nanopore NB01–NB24 published barcode set")
BARCODE_PAYLOAD = [
    role("Native Adapter cohesive hook", "[cohesive hook]", "me"),
    role("barcode outer flank", "[barcode-adapter flank]", "r1"),
    Segment("NB01 representative barcode", NB01, "cbc", feature=BARCODE_FEATURE),
    role("barcode inner flank", "[barcode-adapter flank]", "r1"),
    Segment("dT:dA ligation junction", "T", bottom="A"),
    Segment("genomic DNA", "X" * 42, placeholder=True),
]
FINAL_LIBRARY = motor_library(BARCODE_PAYLOAD, adapter="Native Adapter",
                              name="SQK-NBD114.24 pore-entry strand")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = no_primer_ending("Native Adapter")
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing is motor-driven rather than primed; the "
                          "current Native Adapter and motor-bearing strand are undisclosed.")
FINAL_CAPTION = ("One productive pore-entry orientation. ** labels the adapter-hook "
                 "attachment and the barcode-adapter-to-dA-tailed-insert ligation.")
END_PREP_PROGRAM = ((20, 5), (65, 5))
BARCODE_LIGATION_MIN = 20
NATIVE_ADAPTER_LIGATION_MIN = 20


def oligos():
    for name, forward in NATIVE_BARCODES.items():
        yield oligo(f"{name} forward", [Segment("native barcode", forward, "cbc",
                                                feature=BARCODE_FEATURE)])
        yield oligo(f"{name} reverse", [Segment("native barcode reverse", revcomp(forward),
                                                "cbc", feature=BARCODE_FEATURE)])


def final_rows():
    return single_strand_rows(FINAL_LIBRARY, (
        ("motor-loaded sequencing adapter", "Native Adapter cohesive hook", "attachment"),
        ("barcode inner flank", "dT:dA ligation junction", "ligation"),
    ))


def sections():
    barcode = Construct(BARCODE_PAYLOAD[1:], name="barcoded DNA end")
    return [
        ("Repair and dA-tail each genomic DNA sample", dA_tailed_scene(END_PREP).rows(),
         "FFPE repair and Ultra II end prep: 20 °C for 5 min, then 65 °C for 5 min."),
        ("Ligate one dT-tailed Native Barcode per sample", single_strand_rows(barcode,
         (("barcode inner flank", "dT:dA ligation junction", "ligation"),)),
         "Blunt/TA ligase acts for 20 min at room temperature and joins barcode adapters to both dA-tailed insert ends; clean up, quantify and pool samples equimolarly."),
        ("Ligate the pooled library to Native Adapter", final_rows(),
         "The barcode adapter's cohesive end hooks the motor-bearing Native Adapter; Quick T4 ligase acts for 20 min at room temperature."),
    ]
