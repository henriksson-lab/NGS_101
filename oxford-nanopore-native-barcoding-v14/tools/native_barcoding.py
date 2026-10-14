"""Molecular model for Oxford Nanopore SQK-NBD114.24 native barcoding."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from endprep import dA_tailed_scene, repair_and_dA_tail
from nanopore_library import (LA_TOP, NATIVE_ADAPTER_BOTTOM, motor_sequencing_plan,
                              no_primer_ending, role, single_strand_rows)
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role as molecular_role

TITLE = "Oxford Nanopore Native Barcoding V14 — SQK-NBD114.24"
NOTES = "01_native-barcoding-v14.html"
SOURCE = ('Commercial defining source: Oxford Nanopore '
          '<a href="https://nanoporetech.com/document/ligation-sequencing-gdna-native-barcoding-v14-sqk-nbd114-24">'
          'NBE_9169_v114</a>.')
SUMMARY = ("PCR-free sample multiplexing: repair and dA-tail each DNA sample, ligate a "
           "dT-tailed native barcode, pool, then ligate the motor-bearing Native Adapter "
           "through the barcode adapter's cohesive end.")
CAVEAT = ("ONT publishes both Native Adapter strands. Barcode-adapter flanks, the "
          "cohesive hook, complete motor-bearing assembly and motor attachment remain "
          "undisclosed.")

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
FINAL_LIBRARY = Construct([
    role("motor-loaded Native Adapter assembly", "[motor assembly]", "r2"),
    Segment("Native Adapter published top strand", LA_TOP, "r2"),
    *BARCODE_PAYLOAD,
], name="SQK-NBD114.24 pore-entry strand")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "Native Adapter"
SEQUENCING_ENDING = no_primer_ending("Native Adapter")
SEQUENCING_UNAVAILABLE = ("Nanopore sequencing is motor-driven rather than primed; both "
                          "NA strands are published, but the complete motor-bearing "
                          "assembly and barcode-adapter cohesive hook are not.")
FINAL_CAPTION = ("One productive pore-entry orientation. ** labels the adapter-hook "
                 "attachment and the barcode-adapter-to-dA-tailed-insert ligation.")
END_PREP_PROGRAM = ((20, 5), (65, 5))
BARCODE_LIGATION_MIN = 20
NATIVE_ADAPTER_LIGATION_MIN = 20
BARCODE_END = reported_adapter_end(
    "Native Barcode adapter", kind="cohesive", polarity="3-prime",
    terminal_sequence="T",
    statement="The vendor workflow specifies a dT-tailed barcode adapter for dA-tailed DNA.")
NATIVE_HOOK = AdapterLigationGeometry.unavailable(
    "The vendor identifies the barcode-to-Native-Adapter attachment but does not disclose the cohesive-hook bases.")


def oligos():
    yield oligo("Native Adapter top", [Segment("Native Adapter top", LA_TOP, "r2")])
    yield oligo("Native Adapter bottom", [Segment("Native Adapter bottom", NATIVE_ADAPTER_BOTTOM, "r2")])
    for name, forward in NATIVE_BARCODES.items():
        yield oligo(f"{name} forward", [Segment("native barcode", forward, "cbc",
                                                feature=BARCODE_FEATURE)])
        yield oligo(f"{name} reverse", [Segment("native barcode reverse", revcomp(forward),
                                                "cbc", feature=BARCODE_FEATURE)])


def final_rows():
    return single_strand_rows(FINAL_LIBRARY, (
        ("Native Adapter published top strand", "Native Adapter cohesive hook", "attachment"),
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


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "genomic DNA samples", MaterialKind.MOLECULAR,
                  state=molecule("genomic DNA", (molecular_role("genomic DNA"),)),
                  primary=True)
    repair = g.input("repair_mix", "FFPE repair and end-prep enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed genomic DNA",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Repair DNA ends and add one 3-prime dA",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    barcode = g.input("barcode_adapter", "dT-tailed Native Barcode adapter",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (barcode known; adapter flanks proprietary)")
    ligase_1 = g.input("barcode_ligase", "barcode ligase", MaterialKind.CATALYST)
    barcoded = g.product("barcoded", "native-barcoded DNA samples",
                         MaterialKind.MOLECULAR, inferred=True,
                         source=SOURCE + " (adapter flanks proprietary)")
    g.react("barcode", "Ligate one Native Barcode adapter per sample",
            StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (barcode, InputRole.STRUCTURAL),
             (ligase_1, InputRole.CATALYST)), (barcoded,), geometry=BARCODE_END)
    pool_condition = g.input("pool_condition", "equimolar pooling condition",
                             MaterialKind.CONDITION)
    pooled = g.product("pooled", "equimolar pool of native-barcoded samples",
                       MaterialKind.MOLECULAR, state=barcoded.state)
    excluded = g.product("excluded", "barcode-library material excluded from the pool",
                         MaterialKind.WASTE)
    g.react("pool", "Quantify and combine barcode libraries equimolarly",
            StepKind.SAMPLING,
            ((barcoded, InputRole.TEMPLATE), (pool_condition, InputRole.CONDITION)),
            (pooled, excluded), changes_molecule=False)
    g.discard(excluded, "not included in the equimolar sequencing pool")
    native = g.input("native_adapter", "motor-loaded Native Adapter",
                     MaterialKind.ADAPTER, inferred=True,
                     source=SOURCE + " (NA strands published; cohesive hook/motor assembly undisclosed)")
    ligase_2 = g.input("native_ligase", "Native Adapter ligase",
                       MaterialKind.CATALYST)
    library = g.product("library", "motor-adapted native-barcode library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly undisclosed)")
    g.react("native_adapter_ligation", "Attach Native Adapter to the pooled library",
            StepKind.LIGATION,
            ((pooled, InputRole.TEMPLATE), (native, InputRole.STRUCTURAL),
             (ligase_2, InputRole.CATALYST)), (library,), geometry=NATIVE_HOOK)
    g.sequence(library, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair_mix", "end_prep", "barcode_adapter",
                              "barcode_ligase", "barcode", "pool_condition", "pool",
                              "native_adapter", "native_ligase",
                              "native_adapter_ligation", "library"), SOURCE)
    return g.require_publishable()
