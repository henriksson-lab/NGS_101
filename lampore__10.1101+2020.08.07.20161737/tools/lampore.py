"""Source-backed models for named isothermal-amplification sequencing protocols."""
from __future__ import annotations

from types import SimpleNamespace

from chemdraw import feature
from isothermal import lamp_product, lamp_topology
from nanopore_library import (motor_library, motor_sequencing_plan, no_primer_ending,
                              role as nanopore_role)
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, StepKind,
                            Terminal)
from protocol_model import molecule, role
import seqprimers as sp



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
    primers = g.input("lamp_primers",
                      "10x pool of four six-primer sets with one barcoded FIP per target "
                      "(HPLC-purified; actin primers at lower concentration)",
                      MaterialKind.OLIGO, source=PREPRINT)
    enzymes = g.input("lamp_enzymes", "2x LAMP Master Mix (NEB E1700)",
                      MaterialKind.CATALYST, source=PREPRINT)
    lamp_condition = g.input("lamp_condition",
                             "25 µl 2x mix + 5 µl primer pool + 20 µl RNA; 65 °C for 35 min, then 80 °C "
                             "for 5 min (JCM summarises 65–80 °C, 40 min)",
                             MaterialKind.CONDITION, source=PREPRINT)
    lamp = g.product("lamp_product", "sample-barcoded multiplex RT-LAMP concatemers",
                     MaterialKind.MOLECULAR,
                     state=lamp_product("ORF1a, E, N or ACTB target-derived DNA",
                                        "repeated FIP sample barcode"))
    g.react("rt_lamp", "Amplify four targets and repeatedly copy the FIP barcode",
            StepKind.EXTENSION,
            ((sample, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
             (enzymes, InputRole.CATALYST), (lamp_condition, InputRole.CONDITION)), (lamp,),
            source=PREPRINT)
    pool_condition = g.input("pool_condition",
                             "pool the eight FIP-barcoded reactions of each plate column "
                             "(12 column pools)",
                             MaterialKind.CONDITION, source=PREPRINT)
    pooled = g.product("pooled_lamp", "column pool of eight FIP-barcoded reactions",
                       MaterialKind.MOLECULAR, state=lamp.state)
    g.react("pool", "Pool the eight FIP-barcoded reactions in each column", StepKind.GENERIC,
            ((lamp, InputRole.TEMPLATE), (pool_condition, InputRole.CONDITION)), (pooled,),
            changes_molecule=False)
    rapid = g.input("rapid_barcode_adapter",
                    "one of twelve SQK-RBK004 rapid barcodes (2.5 µl with 1 µl column pool "
                    "in 10 µl)",
                    MaterialKind.ADAPTER, inferred=True,
                    source=PREPRINT + " (commercial barcode and tag bases not disclosed)")
    transposase = g.input("transposase", "rapid-barcode transposase (in the kit barcode)",
                          MaterialKind.CATALYST,
                          source=PREPRINT + " (Introduction: the ONT rapid kit uses a "
                          "transposase); " + JCM)
    rapid_condition = g.input("rapid_condition", "30 °C for 2 min, then 80 °C for 2 min",
                              MaterialKind.CONDITION, source=PREPRINT)
    tagged = g.product("tagged", "rapid-barcoded LAMP fragments", MaterialKind.MOLECULAR,
                       state=molecule("rapid-barcoded LAMP fragment", (
                           role("rapid barcode", tag="cbc", molecular_feature=LAMPORE_RAPID_FEATURE),
                           role("transposase attachment tag"),
                           role("target and repeated FIP barcode"))), inferred=True,
                       source=source + " (commercial tag bases not disclosed)")
    off_products = g.product("off_products",
                             "LAMP fragments without a rapid-barcode attachment end",
                             MaterialKind.MOLECULAR, inferred=True,
                             source=PREPRINT + " (such fragments are not described; inferred "
                             "from transposase fragmentation)")
    g.react("rapid_barcode", "Fragment and add the second, rapid barcode",
            StepKind.TAGMENTATION,
            ((pooled, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL),
             (transposase, InputRole.CATALYST), (rapid_condition, InputRole.CONDITION)),
            (tagged, off_products))
    g.mark_terminal(Terminal(off_products.id, Disposition.UNRESOLVED,
        "Fragments lacking a rapid-barcode attachment end are not assigned a physical fate: "
        "all reactions are pooled and the 0.8x AMPure step is not described as removing them. "
        "The preprint states only that all LamPORE reads carry an ONT barcode at one end "
        "(Barcode demultiplexing). Its ~70% read selection concerns reads lacking a LAMP "
        "barcode, a different, computational exclusion, and does not assign this fate."))
    pooled_all = g.product("pooled_library", "single pool of all rapid-barcoded columns",
                           MaterialKind.MOLECULAR, state=tagged.state, inferred=True,
                           source=PREPRINT)
    g.react("pool_all", "Pool all rapid-barcoded columns into one library without normalization",
            StepKind.GENERIC, ((tagged, InputRole.TEMPLATE),), (pooled_all,),
            changes_molecule=False, source=PREPRINT + "; " + JCM)
    beads = g.input("bead_cleanup",
                    "0.8x AMPure beads, fresh 80% ethanol wash, elution in 15 µl EB",
                    MaterialKind.CONDITION, source=PREPRINT)
    cleaned = g.product("cleaned_library", "bead-purified pooled library",
                        MaterialKind.MOLECULAR, state=tagged.state, inferred=True,
                        source=PREPRINT)
    cleanup_waste = g.product("cleanup_waste", "material not retained by the bead cleanup",
                              MaterialKind.WASTE, source=PREPRINT)
    g.react("cleanup", "0.8x AMPure cleanup of the pooled library", StepKind.CLEANUP,
            ((pooled_all, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
            (cleaned, cleanup_waste), source=PREPRINT)
    g.discard(cleanup_waste, "not retained by the magnetic bead cleanup")
    motor = g.input("motor_adapter", "1 µl rapid adapter (RAP) with 11 µl eluate",
                    MaterialKind.ADAPTER, inferred=True,
                    source=PREPRINT + " (motor-adapter bases not disclosed)")
    motor_condition = g.input("motor_condition", "5 min at room temperature",
                              MaterialKind.CONDITION, source=PREPRINT)
    library = g.product("library", "motor-adapted LamPORE library", MaterialKind.LIBRARY,
                        state=LAMPORE_LIBRARY, inferred=True,
                        source=source + " (commercial adapter bases not disclosed)")
    g.react("attach_motor", "Attach the motor-loaded sequencing adapter",
            StepKind.ADAPTER_ATTACHMENT,
            ((cleaned, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
             (motor_condition, InputRole.CONDITION)), (library,),
            source=PREPRINT,
            reason=("Preprint Methods add RAP after the AMPure cleanup; the attachment "
                    "chemistry and adapter bases are commercial and not disclosed."))
    g.sequence(library, motor_sequencing_plan(
        LAMPORE_LIBRARY,
        limitation=("Commercial rapid barcode, attachment tag and motor-adapter bases are "
                    "unpublished. Preprint Methods: one MinION flow cell, 1 h (JCM: on a "
                    "GridION, 60 min). Calls: invalid <50 classified reads; negative <20, "
                    "inconclusive 20–49, positive >=50 SARS-CoV-2 reads. Every read carries an ONT "
                    "barcode at one end; selecting reads that also contain a LAMP barcode "
                    "retains about 70% of reads (preprint Barcode demultiplexing)."),
        read_uses={"Nanopore strand": "extract_then_include"}))
    g.declare_reviewed_scope((
        "sample", "lamp_primers", "lamp_enzymes", "lamp_condition", "rt_lamp",
        "pool_condition", "pool", "rapid_barcode_adapter", "rapid_barcode", "transposase",
        "rapid_condition", "pool_all", "bead_cleanup", "cleanup", "motor_adapter",
        "motor_condition",
        "attach_motor", "library"), source)
    return g.require_publishable()



PREPRINT = ('<a href="https://doi.org/10.1101/2020.08.07.20161737">James et al., '
            '<i>medRxiv</i> (2020) v3</a>, Methods 1')
JCM = ('<a href="https://doi.org/10.1128/JCM.03271-20">Peto et al., <i>J Clin Microbiol</i> '
       '(2021)</a>, Methods (LamPORE protocol 1 July 2020, v1 rev4)')
source = ('<a href="https://doi.org/10.1101/2020.08.07.20161737">'
          'James et al., <i>medRxiv</i> (2020)</a> (defining preprint v3, Methods 1–2) and '
          + JCM + '.')
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
