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



OLINK_FBC = feature("olink_forward_assay_barcode", "feature_barcode", "whitelist",
                    group="olink_assay", part="forward",
                    whitelist="proprietary Olink Explore assay barcode pairs")
OLINK_RBC = feature("olink_reverse_assay_barcode", "feature_barcode", "whitelist",
                    group="olink_assay", part="reverse",
                    whitelist="proprietary Olink Explore assay barcode pairs")
OLINK_SAMPLE = feature("olink_sample_index", "sample_index", "whitelist",
                       whitelist="proprietary Olink Explore sample indexes")

# Olink's WO2021191442A1 sequence listing publishes the adapter portions carried by
# its PCR primers (SEQ ID 1 and 2), but not the complete primers.  Validate them
# against the independently sourced Illumina constants rather than silently treating
# a similarly named adapter as identical.
OLINK_P5_ADAPTER = "AATGATACGGCGACCACCGA"
OLINK_P7_ADAPTER = "CAAGCAGAAGACGGCATACGAGAT"
if not il.TRUSEQ_P5_FULL.startswith(OLINK_P5_ADAPTER):
    raise ValueError("Olink patent P5 sequence is inconsistent with canonical Illumina end")
if OLINK_P7_ADAPTER != il.P7:
    raise ValueError("Olink patent P7 sequence is inconsistent with canonical Illumina end")

# Wik et al. report a 148-bp product and a 66-cycle Read 1 that starts immediately
# after Rd1SP and covers every inline identifier.  The complete canonical Illumina
# P5/Rd1SP end is 58 nt; together with the patent's 24-nt P7 end, this accounts for
# the product exactly.  Keep this arithmetic beside construction so a future adapter
# revision cannot silently turn the labelled 148-bp molecule into a different length.
OLINK_PRODUCT_BP = 148
OLINK_READ_CYCLES = 66
OLINK_P5_TO_RD1_LINKER = il.TRUSEQ_P5_FULL[
    len(OLINK_P5_ADAPTER):-len(il.TRUSEQ_READ1)]
if not OLINK_P5_TO_RD1_LINKER:
    raise ValueError("canonical Olink P5/Rd1SP end has no intervening sequence")
OLINK_INLINE_BP = OLINK_PRODUCT_BP - len(il.TRUSEQ_P5_FULL) - len(OLINK_P7_ADAPTER)
if OLINK_INLINE_BP != OLINK_READ_CYCLES:
    raise ValueError("published Olink product length and Read 1 length are inconsistent")

OLINK_LIBRARY = Construct([
    Segment("P5 adapter (patent SEQ ID 1)", OLINK_P5_ADAPTER, "p5"),
    Segment("canonical Illumina sequence between P5 and Rd1SP",
            OLINK_P5_TO_RD1_LINKER, "p5"),
    Segment("Read 1 sequencing-primer site", il.TRUSEQ_READ1, "r1"),
    role("forward assay barcode", tag="cbc", molecular_feature=OLINK_FBC),
    role("probe-arm hybridization junction"),
    role("reverse assay barcode", tag="cbc", molecular_feature=OLINK_RBC),
    role("inline sample index", tag="cbc", molecular_feature=OLINK_SAMPLE),
    Segment("P7 reverse complement (patent SEQ ID 2)",
            revcomp(OLINK_P7_ADAPTER), "p7"),
], name="148-bp Olink Explore PEA library")
OLINK_SEQ_PRIMERS = (sp.TRUSEQ["R1"],)


def _olink_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("Olink Explore / NGS-PEA", authority=source)
    sample = g.input("sample", "serially diluted plasma or serum sample",
                     MaterialKind.BIOLOGICAL,
                     state=molecule("sample", (role("target proteins"),)), primary=True)
    antibodies = g.input("antibody_probes", "matched antibody probe pair",
                         MaterialKind.SELECTION)
    probe_oligos = g.input("probe_oligos", "conjugated forward and reverse PEA oligos",
                           MaterialKind.OLIGO)
    bound = g.product("bound_pair", "target-bound matched PEA probe pair",
                      MaterialKind.MOLECULAR, state=proximity_extension_pair())
    unbound = g.product("unbound_probes", "unbound and mismatched probe population",
                        MaterialKind.WASTE)
    g.react("bind", "Bind matched oligo-antibody probes to the same protein",
            StepKind.BINDING,
            ((sample, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION),
             (probe_oligos, InputRole.STRUCTURAL)), (bound, unbound))
    g.discard(unbound, "Only correctly proximal probe pairs generate amplifiable reporters.")

    polymerase = g.input("extension_polymerase", "PEA extension polymerase",
                         MaterialKind.CATALYST)
    reporter = g.product("reporter", "proximity-extension reporter",
                         MaterialKind.MOLECULAR, state=proximity_extension_product())
    g.react("extend", "Extend across the hybridized probe arms",
            StepKind.EXTENSION,
            ((bound, InputRole.TEMPLATE), (polymerase, InputRole.CATALYST)), (reporter,))

    pcr1_f = g.input("pcr1_f", "universal PCR1 forward primer", MaterialKind.OLIGO)
    pcr1_r = g.input("pcr1_r", "universal PCR1 reverse primer", MaterialKind.OLIGO)
    pcr1_enzyme = g.input("pcr1_enzyme", "Olink PCR1 enzyme", MaterialKind.CATALYST)
    preamp = g.product("preamp", "25-cycle preamplified PEA reporters",
                       MaterialKind.MOLECULAR, state=reporter.state)
    g.react("pcr1", "Preamplify every extended reporter with universal primers",
            StepKind.PCR,
            ((reporter, InputRole.TEMPLATE), (pcr1_f, InputRole.STRUCTURAL),
             (pcr1_r, InputRole.STRUCTURAL), (pcr1_enzyme, InputRole.CATALYST)),
            (preamp,), geometry=PCRGeometry.unavailable(
                OLINK_LIBRARY, reason="Olink publishes primer roles but not PCR1 oligo bases."))

    pool_rule = g.input("block_pooling", "equal-volume abundance-block pooling rule",
                        MaterialKind.CONDITION)
    pooled = g.product("pooled_preamp", "sample-level pool of four abundance blocks",
                       MaterialKind.MOLECULAR, state=preamp.state)
    g.react("pool_blocks", "Pool PCR1 products from abundance blocks A–D",
            StepKind.GENERIC,
            ((preamp, InputRole.TEMPLATE), (pool_rule, InputRole.CONDITION)), (pooled,),
            changes_molecule=False)

    pcr2_f = g.input("pcr2_f", "PCR2 forward primer (same species as PCR1)",
                     MaterialKind.OLIGO)
    pcr2_r = g.input("pcr2_r", "sample-indexed P7 reverse primer", MaterialKind.OLIGO)
    pcr2_enzyme = g.input("pcr2_enzyme", "Olink PCR2 enzyme", MaterialKind.CATALYST)
    library = g.product("library", "148-bp Olink Explore sequencing library",
                        MaterialKind.LIBRARY, state=OLINK_LIBRARY)
    g.react("pcr2", "Add the inline sample index and P7 end in ten PCR cycles",
            StepKind.PCR,
            ((pooled, InputRole.TEMPLATE), (pcr2_f, InputRole.STRUCTURAL),
             (pcr2_r, InputRole.STRUCTURAL), (pcr2_enzyme, InputRole.CATALYST)),
            (library,), geometry=PCRGeometry.unavailable(
                OLINK_LIBRARY, reason="Olink publishes PCR2 architecture but not index-primer bases."))
    g.sequence(library, SequencingPlan(
        OLINK_LIBRARY, OLINK_SEQ_PRIMERS, ("Read 1",), {"Read 1": 66},
        {"Read 1": "extract_and_demultiplex"}))
    g.declare_reviewed_scope((
        "sample", "antibody_probes", "probe_oligos", "bind",
        "extension_polymerase", "extend", "pcr1_f", "pcr1_r", "pcr1_enzyme",
        "pcr1", "block_pooling", "pool_blocks", "pcr2_f", "pcr2_r",
        "pcr2_enzyme", "pcr2", "library"), source)
    return g.require_publishable()




source = ('<a href="https://doi.org/10.1016/j.mcpro.2021.100168">'
          'Wik et al., <i>Molecular &amp; Cellular Proteomics</i> (2021)</a>.')
spec = SimpleNamespace(
    TITLE="Olink Explore — proximity extension assay with NGS readout",
    NOTES="01_olink-explore.html", SOURCE=source,
    SUMMARY=("Matched antibody–oligo pairs co-bind a protein, hybridize in "
             "proximity and template polymerase extension. Two PCR stages retain "
             "the paired assay barcodes and add an inline sample index for counting."),
    CAVEAT=("Olink publishes the 148-bp element order and run layout. Complete "
            "probe, assay-barcode, hybridization-arm, sample-index and PCR-primer "
            "sequences remain unverified in the source audit. Only the standard "
            "Illumina end is rendered base by base; the other regions remain "
            "structured roles."),
    TOPOLOGY_PANELS=ProximityExtensionTopology().panels(),
    oligos=lambda: (), FINAL_LIBRARY=OLINK_LIBRARY,
    FINAL_CAPTION=("Published 148-bp architecture: the canonical overlapping P5 / "
                   "Read 1 site, a 66-nt proprietary inline region containing both "
                   "assay barcodes, the copied probe junction and sample index, then "
                   "P7."),
    SEQ_PRIMERS=OLINK_SEQ_PRIMERS, ACQUIRED_ROLES=("Read 1",),
    READ_LENGTHS={"Read 1": OLINK_READ_CYCLES},
    READ_USES={"Read 1": "extract_and_demultiplex"},
    SHOW_READ_CONTENT=False,
    READ_SPECS=(sp.read_spec(
        "Read 1", OLINK_READ_CYCLES,
        "forward assay barcode + reverse assay barcode + inline sample index",
        "extract_and_demultiplex"),),
    SEQUENCING_INTRO=("A standard Read 1 primer binds the published Rd1SP site. "
                      "The NovaSeq run acquires one 66-cycle read and no separate "
                      "index read; assay and sample identities are inline."))


def protocol_graph() -> ProtocolGraph:
    return _olink_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
