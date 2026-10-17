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
          'NBE_9169_v114 rev X (15 Jul 2026)</a>; barcode flanks and Native Adapter strands from the '
          '<a href="https://nanoporetech.com/document/chemistry-technical-document">Chemistry Technical '
          'Document</a> CHTD_500 v1 rev AV, sections 14–15.')
SUMMARY = ("PCR-free sample multiplexing: repair and dA-tail each DNA sample, clean up and "
           "normalise to equimolar mass, ligate a dT-tailed native barcode, stop with EDTA, "
           "pool and clean up, then ligate the motor-bearing Native Adapter through the "
           "barcode adapter's cohesive end.")
CAVEAT = ("ONT publishes both Native Adapter strands and, as analysis context, the "
          "native-barcode flanks (Chemistry Technical Document section 14). Drawing those "
          "flanks as the ligated barcode-adapter strand, the final flank T as the dT tail, "
          "and the 6-nt 3′ overhang of the published Native Adapter duplex as the cohesive "
          "hook are model derivations, not a disclosed oligo inventory. The motor-bearing "
          "assembly and motor attachment remain undisclosed.")

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
# CHTD_500 v1 rev AV section 14 ("Barcode flanking sequences"), written 5'->3' around
# the forward barcode and its reverse complement respectively.
NB_FLANK_FORWARD = ("AAGGTTAA", "CAGCACCT")
NB_FLANK_REVERSE = ("GGTGCTG", "TTAACCTTAGCAAT")
# Derived: the reverse 5' flank is the reverse complement of the forward 3' flank minus
# its final T, so that T is drawn as the single dT tail; the reverse 3' flank extends six
# bases past the forward 5' flank, and those bases are complementary to the 3' overhang
# left on the published Native Adapter top strand by its shorter bottom strand.
if revcomp(NB_FLANK_FORWARD[1][:-1]) != NB_FLANK_REVERSE[0]:
    raise ValueError("published barcode flanks no longer imply a single dT tail")
NA_HOOK_OVERHANG = LA_TOP[len(NATIVE_ADAPTER_BOTTOM):]
if (revcomp(NATIVE_ADAPTER_BOTTOM) != LA_TOP[:len(NATIVE_ADAPTER_BOTTOM)]
        or NB_FLANK_REVERSE[1] != revcomp(NA_HOOK_OVERHANG + NB_FLANK_FORWARD[0])):
    raise ValueError("published flank and Native Adapter strands no longer form the hook")
INSERT = Construct([Segment("genomic DNA", "X" * 42, placeholder=True)], name="genomic DNA")
END_PREP = repair_and_dA_tail(INSERT)
BARCODE_FEATURE = feature("native_barcode", "sample_index", "whitelist",
                          whitelist="Oxford Nanopore NB01–NB24 published barcode set")
BARCODE_PAYLOAD = [
    Segment("barcode outer flank (published)", NB_FLANK_FORWARD[0], "r1"),
    Segment("NB01 representative barcode", NB01, "cbc", feature=BARCODE_FEATURE),
    Segment("barcode inner flank (published)", NB_FLANK_FORWARD[1][:-1], "r1"),
    Segment("dT:dA ligation junction", NB_FLANK_FORWARD[1][-1], bottom="A"),
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
                          "NA strands and the barcode flanks are published, but the "
                          "complete motor-bearing assembly is not, and the hook pairing "
                          "is derived rather than disclosed.")
FINAL_CAPTION = ("One productive pore-entry orientation. ** labels the derived "
                 "Native-Adapter-to-barcode hook ligation and the "
                 "barcode-adapter-to-dA-tailed-insert ligation. Flank bases are printed "
                 "by ONT as barcode context; their placement in this strand and the dT "
                 "split are derived.")
END_PREP_PROGRAM = ((20, 5), (65, 5))
BARCODE_LIGATION_MIN = 20
NATIVE_ADAPTER_LIGATION_MIN = 20
BARCODE_END = reported_adapter_end(
    "Native Barcode adapter", kind="cohesive", polarity="3-prime",
    terminal_sequence="T",
    statement="The vendor workflow specifies a dT-tailed barcode adapter for dA-tailed DNA.")
NATIVE_HOOK = AdapterLigationGeometry.unavailable(
    "The vendor states that the barcode adapter's cohesive end hooks the Native Adapter "
    "but prints no hook bases. Derived candidate only: the published NA top strand "
    f"keeps a 3-prime {NA_HOOK_OVERHANG} overhang beyond the published bottom strand, "
    f"and the published reverse barcode flank ends in its complement "
    f"{revcomp(NA_HOOK_OVERHANG)}.")


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
        ("Native Adapter published top strand", "barcode outer flank (published)",
         "derived hook ligation"),
        ("barcode inner flank (published)", "dT:dA ligation junction", "ligation"),
    ))


def sections():
    barcode = Construct(BARCODE_PAYLOAD, name="barcoded DNA end")
    return [
        ("Repair and dA-tail each genomic DNA sample", dA_tailed_scene(END_PREP).rows(),
         "FFPE repair and Ultra II end prep: 20 °C for 5 min, then 65 °C for 5 min; "
         "AMPure XP cleanup with 80% ethanol washes; Qubit; equimolar mass taken forward."),
        ("Ligate one dT-tailed Native Barcode per sample", single_strand_rows(barcode,
         (("barcode inner flank (published)", "dT:dA ligation junction", "ligation"),)),
         "Blunt/TA Ligase Master Mix acts for 20 min at room temperature; EDTA stops the "
         "reaction; all barcoded samples are pooled and cleaned with 0.4X AMPure XP and "
         "Short Fragment Buffer washes."),
        ("Ligate the pooled library to Native Adapter", final_rows(),
         "The barcode adapter's cohesive end hooks the motor-bearing Native Adapter; Quick "
         "T4 DNA Ligase acts for 20 min at room temperature; AMPure XP cleanup with LFB or "
         "SFB washes (no ethanol) and Elution Buffer."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("dna", "genomic DNA sample (400 ng for >4 barcodes, 1000 ng for ≤4)",
                  MaterialKind.MOLECULAR,
                  state=molecule("genomic DNA", (molecular_role("genomic DNA"),)),
                  primary=True)
    repair = g.input("repair_mix", "NEBNext FFPE DNA Repair Mix and Ultra II End-prep enzymes",
                     MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed genomic DNA",
                       MaterialKind.MOLECULAR, state=END_PREP)
    g.react("end_prep", "Repair DNA ends and add one 3-prime dA (20 °C 5 min, 65 °C 5 min)",
            StepKind.END_REPAIR,
            ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (tailed,))
    beads_1 = g.input("end_prep_beads", "15 µl AMPure XP Beads, 80% ethanol washes, water elution",
                      MaterialKind.SELECTION)
    tailed_clean = g.product("tailed_clean", "bead-purified dA-tailed DNA",
                             MaterialKind.MOLECULAR, state=END_PREP)
    end_prep_waste = g.product("end_prep_waste", "end-prep supernatant and ethanol washes",
                               MaterialKind.WASTE)
    g.react("end_prep_cleanup", "Clean up each end-prep reaction on AMPure XP beads",
            StepKind.CLEANUP,
            ((tailed, InputRole.TEMPLATE), (beads_1, InputRole.SELECTION)),
            (tailed_clean, end_prep_waste))
    g.discard(end_prep_waste, "unbound end-prep reaction components removed with the supernatant")
    normalise = g.input("equimolar_rule", "Qubit quantification; equimolar mass per sample",
                        MaterialKind.CONDITION)
    normalised = g.product("normalised", "equimolar end-prepped DNA for barcode ligation",
                           MaterialKind.MOLECULAR, state=END_PREP)
    leftover = g.product("leftover", "end-prepped DNA not taken into barcode ligation",
                         MaterialKind.WASTE)
    g.react("normalise", "Take an equimolar mass of each sample into barcode ligation",
            StepKind.SAMPLING,
            ((tailed_clean, InputRole.TEMPLATE), (normalise, InputRole.CONDITION)),
            (normalised, leftover), changes_molecule=False)
    g.discard(leftover, "not carried into the barcode ligation")
    barcode = g.input("barcode_adapter", "dT-tailed Native Barcode adapter (one of NB01–NB24)",
                      MaterialKind.ADAPTER, inferred=True,
                      source=SOURCE + " (barcode and flanks printed; adapter duplex derived)")
    ligase_1 = g.input("barcode_ligase", "NEB Blunt/TA Ligase Master Mix", MaterialKind.CATALYST)
    barcoded = g.product("barcoded", "native-barcoded DNA sample",
                         MaterialKind.MOLECULAR, inferred=True,
                         source=SOURCE + " (adapter duplex derived from printed flanks)")
    g.react("barcode", "Ligate one Native Barcode adapter per sample (20 min, room temperature)",
            StepKind.LIGATION,
            ((normalised, InputRole.TEMPLATE), (barcode, InputRole.STRUCTURAL),
             (ligase_1, InputRole.CATALYST)), (barcoded,), geometry=BARCODE_END)
    edta = g.input("edta", "EDTA (4 µl per well)", MaterialKind.CONDITION)
    stopped = g.product("stopped", "EDTA-stopped barcode ligation",
                        MaterialKind.MOLECULAR, state=barcoded.state, inferred=True,
                        source=SOURCE + " (adapter duplex derived from printed flanks)")
    g.react("stop", "Stop each barcode ligation with EDTA", StepKind.QUENCHING,
            ((barcoded, InputRole.TEMPLATE), (edta, InputRole.CONDITION)), (stopped,),
            changes_molecule=False)
    others = g.input("other_barcoded", "other EDTA-stopped barcode ligations "
                     "(different NB barcodes, processed in parallel)", MaterialKind.MOLECULAR)
    pooled = g.product("pooled", "pool of all barcoded samples",
                       MaterialKind.MOLECULAR, state=barcoded.state, inferred=True,
                       source=SOURCE + " (adapter duplex derived from printed flanks)")
    g.react("pool", "Pool the entire volume of every barcoded sample", StepKind.POOL,
            ((stopped, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)), (pooled,),
            changes_molecule=False)
    beads_2 = g.input("pool_beads", "0.4X AMPure XP Beads, 700 µl SFB washes, "
                      "water elution at 37 °C", MaterialKind.SELECTION)
    pooled_clean = g.product("pooled_clean", "bead-purified barcoded pool",
                             MaterialKind.MOLECULAR, state=barcoded.state, inferred=True,
                             source=SOURCE + " (adapter duplex derived from printed flanks)")
    pool_waste = g.product("pool_waste", "pool supernatant (unligated barcode adapters, "
                           "ligase, EDTA) and SFB washes", MaterialKind.WASTE)
    g.react("pool_cleanup", "Clean up the pooled barcode ligations", StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (beads_2, InputRole.SELECTION)),
            (pooled_clean, pool_waste))
    g.discard(pool_waste, "unbound material removed with the bead supernatant and washes")
    native = g.input("native_adapter", "motor-loaded Native Adapter",
                     MaterialKind.ADAPTER, inferred=True,
                     source=SOURCE + " (NA strands published; motor assembly undisclosed)")
    ligase_2 = g.input("native_ligase", "NEBNext Quick Ligation Reaction Buffer and Quick T4 DNA Ligase",
                       MaterialKind.CATALYST)
    library = g.product("library", "motor-adapted native-barcode library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly undisclosed)")
    g.react("native_adapter_ligation", "Attach Native Adapter to the pooled library "
            "(20 min, room temperature)", StepKind.LIGATION,
            ((pooled_clean, InputRole.TEMPLATE), (native, InputRole.STRUCTURAL),
             (ligase_2, InputRole.CATALYST)), (library,), geometry=NATIVE_HOOK)
    beads_3 = g.input("library_beads", "20 µl AMPure XP Beads, 125 µl LFB or SFB washes "
                      "(not ethanol), Elution Buffer at 37 °C", MaterialKind.SELECTION)
    library_clean = g.product("library_clean", "bead-purified native-barcode library",
                              MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                              source=SOURCE + " (complete motor-bearing assembly undisclosed)")
    library_waste = g.product("library_waste", "ligation supernatant (unligated Native "
                              "Adapter) and LFB/SFB washes", MaterialKind.WASTE)
    g.react("library_cleanup", "Clean up the adapted library; LFB enriches ≥3 kb, SFB keeps all sizes",
            StepKind.CLEANUP,
            ((library, InputRole.TEMPLATE), (beads_3, InputRole.SELECTION)),
            (library_clean, library_waste))
    g.discard(library_waste, "unbound material removed with the bead supernatant and washes")
    g.sequence(library_clean, motor_sequencing_plan(
        FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("dna", "repair_mix", "end_prep", "end_prep_beads",
                              "end_prep_cleanup", "equimolar_rule", "normalise",
                              "barcode_adapter", "barcode_ligase", "barcode", "edta", "stop",
                              "other_barcoded", "pool", "pool_beads", "pool_cleanup",
                              "native_adapter", "native_ligase", "native_adapter_ligation",
                              "library_beads", "library_cleanup", "library_clean"), SOURCE)
    return g.require_publishable()
