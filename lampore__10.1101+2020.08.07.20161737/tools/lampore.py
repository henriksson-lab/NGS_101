"""Source-backed models for named isothermal-amplification sequencing protocols."""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
from batch_ngs import nextera_library, truseq_amplicon_library, truseq_library
from chemdraw import Construct, Segment, feature, oligo, revcomp
from isothermal import (lamp_product, lamp_topology, mda_product, proximity_extension_pair,
                        proximity_extension_product, spia_product)
from nanopore_library import (motor_library, motor_sequencing_plan, no_primer_ending,
                              role as nanopore_role)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
import seqprimers as sp
from proximity_topology import ProximityExtensionTopology



LAMPORE_FIP_FEATURE = feature(
    "lampore_fip_barcode", "sample_index", "whitelist",
    whitelist="eight commercial LamPORE FIP barcode variants")
LAMPORE_RAPID_FEATURE = feature(
    "lampore_rapid_barcode", "sample_index", "whitelist",
    whitelist="twelve Oxford Nanopore rapid barcode variants")
LAMPORE_LIBRARY = motor_library([
    nanopore_role("rapid barcode", "[rapid barcode]", "cbc",
                  feature=LAMPORE_RAPID_FEATURE),
    nanopore_role("transposase attachment tag", "[transposase attachment tag]", "me"),
    nanopore_role("barcoded LAMP concatemer", "[repeated target + FIP barcode]", "cbc",
                  feature=LAMPORE_FIP_FEATURE),
], adapter="rapid sequencing adapter", name="LamPORE sequenceable strand")


def _lampore_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("LamPORE", authority=source)
    sample = g.input("sample", "extracted sample RNA", MaterialKind.BIOLOGICAL,
                     state=molecule("sample RNA", (role("viral and human RNA"),)),
                     primary=True)
    primers = g.input("lamp_primers", "four multiplexed six-primer sets with one FIP barcode",
                      MaterialKind.OLIGO)
    enzymes = g.input("lamp_enzymes", "reverse transcriptase and strand-displacing polymerase",
                      MaterialKind.CATALYST)
    lamp = g.product("lamp_product", "sample-barcoded multiplex RT-LAMP concatemers",
                     MaterialKind.MOLECULAR,
                     state=lamp_product("ORF1a, E, N or ACTB target-derived DNA",
                                        "repeated FIP sample barcode"))
    g.react("rt_lamp", "Amplify four targets and repeatedly copy the FIP barcode",
            StepKind.EXTENSION,
            ((sample, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
             (enzymes, InputRole.CATALYST)), (lamp,))
    pool_condition = g.input("pool_condition", "pool eight FIP-barcoded reactions by column",
                             MaterialKind.CONDITION)
    pooled = g.product("pooled_lamp", "column pool of eight FIP-barcoded reactions",
                       MaterialKind.MOLECULAR, state=lamp.state)
    g.react("pool", "Pool the eight FIP-barcoded reactions in each column", StepKind.GENERIC,
            ((lamp, InputRole.TEMPLATE), (pool_condition, InputRole.CONDITION)), (pooled,),
            changes_molecule=False)
    rapid = g.input("rapid_barcode_adapter", "one of twelve loaded rapid barcode transposomes",
                    MaterialKind.ADAPTER, inferred=True,
                    source=source + " (commercial barcode and tag bases not disclosed)")
    transposase = g.input("transposase", "rapid-barcode transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "rapid-barcoded LAMP fragments", MaterialKind.MOLECULAR,
                       state=molecule("rapid-barcoded LAMP fragment", (
                           role("rapid barcode", tag="cbc", molecular_feature=LAMPORE_RAPID_FEATURE),
                           role("transposase attachment tag"),
                           role("target and repeated FIP barcode"))), inferred=True,
                       source=source + " (commercial tag bases not disclosed)")
    waste = g.product("off_products", "non-sequenceable transposition products",
                      MaterialKind.WASTE)
    g.react("rapid_barcode", "Fragment and add the second, rapid barcode",
            StepKind.TAGMENTATION,
            ((pooled, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST)), (tagged, waste))
    g.discard(waste, "does not carry the productive rapid-adapter attachment end")
    motor = g.input("motor_adapter", "motor-loaded rapid sequencing adapter",
                    MaterialKind.ADAPTER, inferred=True,
                    source=source + " (motor-adapter bases not disclosed)")
    library = g.product("library", "motor-adapted LamPORE library", MaterialKind.LIBRARY,
                        state=LAMPORE_LIBRARY, inferred=True,
                        source=source + " (commercial adapter bases not disclosed)")
    g.react("attach_motor", "Attach the motor-loaded sequencing adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((tagged, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL)), (library,))
    g.sequence(library, motor_sequencing_plan(
        LAMPORE_LIBRARY,
        limitation="Commercial rapid barcode, attachment tag and motor-adapter bases are unpublished.",
        read_uses={"Nanopore strand": "extract_then_include"}))
    g.declare_reviewed_scope((
        "sample", "lamp_primers", "lamp_enzymes", "rt_lamp", "pool_condition", "pool",
        "rapid_barcode_adapter", "rapid_barcode", "transposase", "motor_adapter",
        "attach_motor", "library"), source)
    return g.require_publishable()



source = ('<a href="https://doi.org/10.1101/2020.08.07.20161737">'
          'James et al., <i>medRxiv</i> (2020)</a>.')
spec = SimpleNamespace(
    TITLE="LamPORE — barcoded RT-LAMP with Nanopore readout",
    NOTES="01_lampore.html", SOURCE=source,
    SUMMARY=("Four primer sets amplify three SARS-CoV-2 targets and an ACTB "
             "control. One of eight FIP barcodes labels each sample before "
             "pooling; one of twelve rapid barcodes labels each pooled column "
             "before motor-driven Nanopore sequencing."),
    CAVEAT=("The defining source publishes the targets and barcode logic but not "
            "the commercial FIP-barcode or rapid-adapter bases. LAMP products are "
            "heterogeneous branched concatemers, so the drawing preserves roles "
            "and repeated content rather than inventing one linear sequence."),
    TOPOLOGY_PANELS=lamp_topology(
        "ORF1a, E, N or ACTB target-derived DNA",
        "repeated FIP sample barcode").panels(),
    oligos=lambda: (), FINAL_LIBRARY=LAMPORE_LIBRARY,
    FINAL_CAPTION=("One productive pore-entry strand after rapid transposase "
                   "barcoding and motor-adapter attachment. Bracketed elements are "
                   "source-established roles whose bases are not disclosed."),
    SEQ_PRIMERS=(), SEQUENCING_ENDING=no_primer_ending("rapid adapter"),
    SEQUENCING_UNAVAILABLE=("Nanopore sequencing is motor-driven; the commercial "
                            "rapid barcode, attachment tag and motor-adapter bases "
                            "are not published."),
    READ_SPECS=(sp.read_spec(
        "Sequencing output", None,
        "rapid barcode + repeated FIP sample barcode + ORF1a, E, N or ACTB target-derived sequence",
        "extract_then_include"),))


def protocol_graph() -> ProtocolGraph:
    return _lampore_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
