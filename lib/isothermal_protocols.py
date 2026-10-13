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


LAMP_F3 = "A+AC+AC+AA+GC+TTTCGGCAG"
LAMP_B3 = "G+AA+AT+TT+GG+ATCTTTGTCATCC"
LAMP_FIP_LEFT = "TGCGGCCAATGTTTGTAATCAG"
LAMP_FIP_BARCODE = "AGCAAACTAC"  # first validated Tier-1 barcode
LAMP_FIP_RIGHT = "CCAAGGAAATTTTGGGGAC"
LAMP_FIP = LAMP_FIP_LEFT + LAMP_FIP_BARCODE + LAMP_FIP_RIGHT
LAMP_BIP = "CGCATTGGCATGGAAGTCACTTTGATGGCACCTGTGTAG"
LAMP_LF = "TTCCTTGTCTGATTAGTTC"
LAMP_LB = "ACCTTCGGGAACGTGGTT"
LAMP_PCR_FWD_SITE = "AACGCTGAAGCGCTGGGGGCAAA"
LAMP_PCR_REV_SITE = "GTTTGTAATCAGTTCCTTGTCTG"
LAMP_PCR_FWD = il.TRUSEQ_READ1 + LAMP_PCR_FWD_SITE
# The source prints this first-stage tail with TGACT rather than canonical GTGACT.
LAMP_PCR_REV_ARM = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
LAMP_PCR_REV = LAMP_PCR_REV_ARM + LAMP_PCR_REV_SITE
LAMP_POOL_I5 = "TCGTGGAGCG"
LAMP_POOL_I7 = "CGCTCAGTTC"
LAMP_POOL_FWD = il.P5 + LAMP_POOL_I5 + il.TRUSEQ_READ1
LAMP_POOL_REV = il.P7 + LAMP_POOL_I7 + il.TRUSEQ_READ2

LAMP_BARCODE_FEATURE = feature(
    "lamp_seq_fip_barcode", "sample_index", "whitelist",
    whitelist="LAMP-Seq Supplementary Data 1 Tier-1 or Tier-2 FIP barcodes")

# A LAMP reaction contains several branched and inverted isomers. This construct is only
# the source-established PCR-selected interval: exact landing sites and the embedded FIP
# barcode are retained, while their intervening order within a selected concatemer is not
# guessed from the six-primer reaction.
LAMP_SELECTED_INSERT = [
    Segment("PCR-N forward landing region", LAMP_PCR_FWD_SITE),
    role("LAMP-derived N-gene sequence before the recovered FIP barcode"),
    Segment("FIP sample barcode (representative)", LAMP_FIP_BARCODE, "cbc",
            feature=LAMP_BARCODE_FEATURE),
    role("LAMP-derived N-gene sequence after the recovered FIP barcode"),
    Segment("PCR-N reverse landing region reverse complement",
            revcomp(LAMP_PCR_REV_SITE)),
]
LAMP_PCR_TEMPLATE = Construct(LAMP_SELECTED_INSERT, name="PCR-selected LAMP concatemer unit")
LAMP_PREAMP = Construct([
    Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
    *LAMP_SELECTED_INSERT,
    Segment("first-PCR reverse arm reverse complement", revcomp(LAMP_PCR_REV_ARM), "r2"),
], name="first-stage LAMP-Seq amplicon")
LAMP_LIBRARY, LAMP_SEQ_PRIMERS = truseq_amplicon_library(
    LAMP_SELECTED_INSERT, "LAMP-Seq dual-index library", i5_nt=10, i7_nt=10)


def _lamp_oligos():
    feature_seg = Segment("FIP sample barcode", LAMP_FIP_BARCODE, "cbc",
                          feature=LAMP_BARCODE_FEATURE)
    yield oligo("N-F3-LNA", [Segment("LNA-modified N-F3", LAMP_F3, placeholder=True)])
    yield oligo("N-B3-LNA", [Segment("LNA-modified N-B3", LAMP_B3, placeholder=True)])
    yield oligo("N-FIP-Barcode (representative Tier-1 barcode)", [
        Segment("F1c", LAMP_FIP_LEFT), feature_seg, Segment("F2", LAMP_FIP_RIGHT)])
    yield oligo("N-BIP", [Segment("B1c + B2", LAMP_BIP)])
    yield oligo("N-LF", [Segment("loop forward", LAMP_LF)])
    yield oligo("N-LB", [Segment("loop backward", LAMP_LB)])
    yield oligo("PCR-N-fwd-10", [Segment("Read 1 tail + N target", LAMP_PCR_FWD)])
    yield oligo("PCR-N-rev-11", [Segment("Read 2 tail + N target", LAMP_PCR_REV)])
    yield oligo("Pool_fwd_1", [Segment("P5 + 10-nt i5 + Read 1", LAMP_POOL_FWD)])
    yield oligo("Pool_rev_1", [Segment("P7 + 10-nt i7 + Read 2", LAMP_POOL_REV)])


def _lamp_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("LAMP-Seq", authority=source)
    sample = g.input("sample", "swab lysate containing SARS-CoV-2 RNA",
                     MaterialKind.BIOLOGICAL,
                     state=molecule("swab lysate", (role("viral RNA"),)), primary=True)
    heat = g.input("heat", "95 °C lysis and inactivation condition", MaterialKind.CONDITION)
    lysate = g.product("heated_lysate", "heat-inactivated RNA-containing lysate",
                       MaterialKind.MOLECULAR, state=sample.state)
    g.react("inactivate", "Heat-inactivate the swab lysate", StepKind.THERMAL_TREATMENT,
            ((sample, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (lysate,))

    primers = g.input("lamp_primers", "six-primer set with barcoded FIP",
                      MaterialKind.OLIGO, state=LAMP_PCR_TEMPLATE)
    bst = g.input("bst", "strand-displacing Bst polymerase", MaterialKind.CATALYST)
    lamp = g.product("lamp_product", "barcoded branched RT-LAMP products",
                     MaterialKind.MOLECULAR, state=lamp_product("SARS-CoV-2 N-gene DNA"))
    g.react("rt_lamp", "Reverse-transcribe and amplify by loop-mediated strand displacement",
            StepKind.EXTENSION,
            ((lysate, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
             (bst, InputRole.CATALYST)), (lamp,))

    dilution = g.input("dilution", "1:40,000 pre-PCR dilution", MaterialKind.CONDITION)
    pooled = g.product("pooled_lamp", "pooled and diluted barcoded LAMP products",
                       MaterialKind.MOLECULAR, state=lamp.state)
    g.react("pool_dilute", "Pool barcoded reactions and dilute before PCR",
            StepKind.GENERIC,
            ((lamp, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)), (pooled,),
            changes_molecule=False)

    pcr_f = g.input("pcr_f", "PCR-N-fwd-10", MaterialKind.OLIGO)
    pcr_r = g.input("pcr_r", "PCR-N-rev-11", MaterialKind.OLIGO)
    polymerase = g.input("pcr_polymerase", "NEBNext PCR master mix", MaterialKind.CATALYST)
    preamp = g.product("preamp", "first-stage LAMP-Seq amplicon",
                       MaterialKind.MOLECULAR, state=LAMP_PREAMP)
    geometry1 = PCRGeometry.exact(
        LAMP_PCR_TEMPLATE,
        PCRPrimer("PCR-N-fwd-10", LAMP_PCR_FWD, anneal_nt=len(LAMP_PCR_FWD_SITE)),
        PCRPrimer("PCR-N-rev-11", LAMP_PCR_REV, anneal_nt=len(LAMP_PCR_REV_SITE)),
        statement="Ludwig et al. Supplementary Data 1 publishes both target-tailed PCR primers.")
    g.react("target_pcr", "Select one N-gene interval from the heterogeneous LAMP products",
            StepKind.PCR,
            ((pooled, InputRole.TEMPLATE), (pcr_f, InputRole.STRUCTURAL),
             (pcr_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (preamp,), geometry=geometry1)

    pool_f = g.input("pool_f", "Pool_fwd_1 representative indexed primer", MaterialKind.OLIGO)
    pool_r = g.input("pool_r", "Pool_rev_1 representative indexed primer", MaterialKind.OLIGO)
    final = g.product("library", "dual-index LAMP-Seq library", MaterialKind.LIBRARY,
                      state=LAMP_LIBRARY)
    geometry2 = PCRGeometry.exact(
        LAMP_PREAMP,
        PCRPrimer("Pool_fwd_1", LAMP_POOL_FWD, anneal_nt=len(il.TRUSEQ_READ1)),
        PCRPrimer("Pool_rev_1", LAMP_POOL_REV, anneal_nt=len(LAMP_PCR_REV_ARM)),
        statement="Supplementary Data 1 publishes 192 pairs of 10-nt pool-index primers; pair 1 is shown.")
    g.react("index_pcr", "Install the 10-nt pool indexes and complete the Illumina ends",
            StepKind.PCR,
            ((preamp, InputRole.TEMPLATE), (pool_f, InputRole.STRUCTURAL),
             (pool_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (final,), geometry=geometry2)
    g.sequence(final, SequencingPlan(
        LAMP_LIBRARY, tuple(LAMP_SEQ_PRIMERS),
        ("Read 1", "Index 1 (i7)", "Index 2 (i5)"),
        {"Read 1": 100, "Index 1 (i7)": 8, "Index 2 (i5)": 8},
        {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex"}))
    g.declare_reviewed_scope((
        "sample", "heat", "inactivate", "lamp_primers", "bst", "rt_lamp",
        "dilution", "pool_dilute", "pcr_f", "pcr_r", "pcr_polymerase",
        "target_pcr", "pool_f", "pool_r", "index_pcr", "library"), source)
    return g.require_publishable()


OLINK_FBC = feature("olink_forward_assay_barcode", "feature_barcode", "whitelist",
                    group="olink_assay", part="forward",
                    whitelist="proprietary Olink Explore assay barcode pairs")
OLINK_RBC = feature("olink_reverse_assay_barcode", "feature_barcode", "whitelist",
                    group="olink_assay", part="reverse",
                    whitelist="proprietary Olink Explore assay barcode pairs")
OLINK_SAMPLE = feature("olink_sample_index", "sample_index", "whitelist",
                       whitelist="proprietary Olink Explore sample indexes")
OLINK_LIBRARY = Construct([
    Segment("P5", il.P5, "p5"),
    Segment("Read 1 sequencing-primer site", il.TRUSEQ_READ1, "r1"),
    role("forward assay barcode", tag="cbc", molecular_feature=OLINK_FBC),
    role("probe-arm hybridization junction"),
    role("reverse assay barcode", tag="cbc", molecular_feature=OLINK_RBC),
    role("inline sample index", tag="cbc", molecular_feature=OLINK_SAMPLE),
    Segment("P7 reverse complement", il.P7_RC, "p7"),
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


OVATION_LIBRARY, OVATION_PRIMERS = truseq_library([
    role("fragment of Ribo-SPIA amplified cDNA")
], "Ovation-derived Illumina library", dual_index=False, inferred_adapters=True)


def _ovation_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("Ovation RNA-Seq V2 / Ribo-SPIA", authority=source)
    rna = g.input("rna", "total RNA", MaterialKind.BIOLOGICAL,
                  state=molecule("total RNA", (role("RNA transcript"),)), primary=True)
    first_primer = g.input("first_primer", "chimeric DNA/RNA first-strand primer",
                           MaterialKind.OLIGO, inferred=True,
                           source=source + " (primer bases proprietary)")
    rt = g.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA/RNA hybrid",
                      MaterialKind.MOLECULAR,
                      state=molecule("first-strand product", (
                          role("RNA template"), role("first-strand cDNA"),
                          role("chimeric primer handle"))))
    g.react("first_strand_rt", "Prime RNA and synthesize first-strand cDNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (first_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    second_primer = g.input("second_primer", "chimeric DNA/RNA second-strand primer",
                            MaterialKind.OLIGO, inferred=True,
                            source=source + " (primer bases proprietary)")
    polymerase = g.input("dna_polymerase", "DNA polymerase", MaterialKind.CATALYST)
    duplex = g.product("cdna_duplex", "double-stranded cDNA with a DNA/RNA SPIA end",
                       MaterialKind.MOLECULAR,
                       state=molecule("SPIA-competent cDNA", (
                           role("double-stranded transcript copy"),
                           role("DNA/RNA chimeric SPIA priming end"))))
    g.react("second_strand", "Synthesize second-strand cDNA and install the SPIA end",
            StepKind.EXTENSION,
            ((first, InputRole.TEMPLATE), (second_primer, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (duplex,))
    spia_primer = g.input("spia_primer", "chimeric SPIA primer", MaterialKind.OLIGO,
                          inferred=True, source=source + " (primer bases proprietary)")
    rnase_h = g.input("rnase_h", "RNase H", MaterialKind.CATALYST)
    amplified = g.product("spia_cdna", "linearly amplified Ribo-SPIA cDNA population",
                          MaterialKind.MOLECULAR, state=spia_product())
    g.react("spia", "Repeatedly cleave, reprime and displace linear cDNA copies",
            StepKind.EXTENSION,
            ((duplex, InputRole.TEMPLATE), (spia_primer, InputRole.STRUCTURAL),
             (rnase_h, InputRole.CATALYST), (polymerase, InputRole.CATALYST)), (amplified,))
    shear = g.input("shear_condition", "mechanical shearing condition", MaterialKind.CONDITION)
    fragmented = g.product("fragmented_cdna", "sheared amplified cDNA", MaterialKind.MOLECULAR,
                           state=molecule("sheared cDNA", (role("transcript-derived cDNA fragment"),)))
    g.react("fragment", "Shear the amplified cDNA for library construction",
            StepKind.FRAGMENTATION,
            ((amplified, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragmented,))
    adapters = g.input("library_reagents", "Encore Rapid Illumina library reagents",
                       MaterialKind.ADAPTER, inferred=True,
                       source=source + " (commercial oligo sequences not reported)")
    library = g.product("library", "Illumina library from Ribo-SPIA cDNA",
                        MaterialKind.LIBRARY, state=OVATION_LIBRARY, inferred=True,
                        source=source + " (canonical shell inferred from compatibility)")
    g.react("make_library", "End-repair, adapter-add and enrich the sheared cDNA",
            StepKind.ADAPTER_ATTACHMENT,
            ((fragmented, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL)), (library,))
    g.sequence(library, SequencingPlan(
        OVATION_LIBRARY, OVATION_PRIMERS,
        ("Read 1", "Index 1 (i7)", "Read 2"),
        {"Read 1": 50, "Index 1 (i7)": 8, "Read 2": 50},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"}))
    g.declare_reviewed_scope((
        "rna", "first_primer", "reverse_transcriptase", "first_strand_rt",
        "second_primer", "dna_polymerase", "second_strand", "spia_primer", "rnase_h",
        "spia", "shear_condition", "fragment", "library_reagents", "make_library",
        "library"), source)
    return g.require_publishable()


REPLIG_LIBRARY, REPLIG_PRIMERS = nextera_library([
    role("fragment sampled from phi29 MDA product")
], "REPLI-g / Nextera library")


def _replig_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("REPLI-g Single Cell", authority=source)
    cell = g.input("cell", "single cell", MaterialKind.BIOLOGICAL,
                   state=molecule("single cell", (role("diploid genomic DNA"),)), primary=True)
    denature = g.input("lysis_denaturation", "lysis and alkaline denaturation reagents",
                       MaterialKind.CONDITION)
    dna = g.product("denatured_dna", "lysed and denatured single-cell genomic DNA",
                    MaterialKind.MOLECULAR,
                    state=molecule("denatured genomic DNA", (role("single genomic strands"),)))
    g.react("lyse_denature", "Lyse the cell and denature genomic DNA",
            StepKind.EXTRACTION, ((cell, InputRole.TEMPLATE),
                                  (denature, InputRole.CONDITION)), (dna,))
    random = g.input("random_primers", "proprietary exonuclease-resistant random hexamers",
                     MaterialKind.OLIGO, inferred=True,
                     source=source + " (primer sequences proprietary)")
    phi29 = g.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    mda = g.product("mda_product", "branched multiple-displacement amplification product",
                    MaterialKind.MOLECULAR, state=mda_product())
    g.react("mda", "Perform processive random-primed phi29 strand-displacement amplification",
            StepKind.EXTENSION,
            ((dna, InputRole.TEMPLATE), (random, InputRole.STRUCTURAL),
             (phi29, InputRole.CATALYST)), (mda,))
    tn5_adapters = g.input("tn5_adapters", "loaded Nextera adapter oligos",
                           MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "Nextera-tagged fragments of MDA DNA",
                       MaterialKind.MOLECULAR,
                       state=molecule("tagged MDA fragment", (
                           role("left mosaic end"), role("MDA-derived genomic fragment"),
                           role("right mosaic end"))))
    waste = g.product("tagging_waste", "non-library tagmentation products", MaterialKind.WASTE)
    g.react("tagment", "Tagment the amplified DNA", StepKind.TAGMENTATION,
            ((mda, InputRole.TEMPLATE), (tn5_adapters, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, waste))
    g.discard(waste, "lacks two productive heterologous Nextera ends")
    p5 = g.input("p5_primer", "Nextera i5/P5 PCR primer", MaterialKind.OLIGO)
    p7 = g.input("p7_primer", "Nextera i7/P7 PCR primer", MaterialKind.OLIGO)
    pcr = g.input("pcr_enzyme", "library PCR polymerase", MaterialKind.CATALYST)
    library = g.product("library", "dual-index Nextera library", MaterialKind.LIBRARY,
                        state=REPLIG_LIBRARY)
    g.react("index_pcr", "Complete flow-cell ends and add dual indexes by PCR",
            StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (p7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                REPLIG_LIBRARY,
                reason="The QIAGEN protocol delegates index-primer identities to the selected Nextera kit."))
    g.sequence(library, SequencingPlan(
        REPLIG_LIBRARY, REPLIG_PRIMERS, tuple(p.role for p in REPLIG_PRIMERS),
        {"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 150},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    g.declare_reviewed_scope((
        "cell", "lysis_denaturation", "lyse_denature", "random_primers", "phi29",
        "mda", "mda_product",
        "tn5_adapters", "tn5", "tagment", "p5_primer", "p7_primer", "pcr_enzyme",
        "index_pcr", "library"), source)
    return g.require_publishable()


PICOPLEX_I5 = feature("picoplex_i5", "sample_index", "whitelist",
                      whitelist="PicoPLEX DNA-seq 48-pair Dual Index Plate")
PICOPLEX_I7 = feature("picoplex_i7", "sample_index", "whitelist",
                      whitelist="PicoPLEX DNA-seq 48-pair Dual Index Plate")
PICOPLEX_LIBRARY = Construct([
    role("Illumina-compatible P5-side end"),
    role("i5 sample index", tag="cbc", length_bp=8, molecular_feature=PICOPLEX_I5),
    role("read-one common handle"), role("single-cell genomic insert"),
    role("read-two common handle"),
    role("i7 sample index", tag="cbc", length_bp=8, molecular_feature=PICOPLEX_I7),
    role("Illumina-compatible P7-side end"),
], name="PicoPLEX DNA-seq dual-index library")
SUREPLEX_LIBRARY, SUREPLEX_PRIMERS = nextera_library([
    role("fragment of SurePlex-amplified genomic DNA")
], "SurePlex / VeriSeq PGS library")


def _picoplex_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("PicoPLEX DNA-seq", authority=source)
    cell = g.input("cell", "one to ten cells or 6–60 pg purified DNA",
                   MaterialKind.BIOLOGICAL,
                   state=molecule("low-input sample", (role("genomic DNA"),)), primary=True)
    lysis = g.input("lysis_reagent", "Cell Extraction Buffer and Extraction Enzyme",
                    MaterialKind.CONDITION)
    dna = g.product("released_dna", "released genomic DNA", MaterialKind.MOLECULAR,
                    state=molecule("released genomic DNA", (role("genomic DNA"),)))
    g.react("lyse", "Lyse cells and release genomic DNA", StepKind.EXTRACTION,
            ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (dna,))
    quasi = g.input("quasi_random_primers", "proprietary quasi-random preamplification primers",
                    MaterialKind.OLIGO, inferred=True,
                    source=source + " (primer sequences proprietary)")
    displacement = g.input("displacement_polymerase", "strand-displacing polymerase",
                           MaterialKind.CATALYST)
    preamp = g.product("preamp_product", "linear displacement DOP-PCR preamplification product",
                       MaterialKind.MOLECULAR,
                       state=molecule("PicoPLEX preamplification population", (
                           role("quasi-random common handle"),
                           role("sampled genomic tract"),
                           role("displaced extension products"))))
    g.react("preamp", "Perform multiple quasi-random priming and linear displacement cycles",
            StepKind.EXTENSION,
            ((dna, InputRole.TEMPLATE), (quasi, InputRole.STRUCTURAL),
             (displacement, InputRole.CATALYST)), (preamp,))
    i5 = g.input("i5_primer", "8-nt i5 indexed amplification primer", MaterialKind.OLIGO,
                 inferred=True, source=source + " (complete primer bases proprietary)")
    i7 = g.input("i7_primer", "8-nt i7 indexed amplification primer", MaterialKind.OLIGO,
                 inferred=True, source=source + " (complete primer bases proprietary)")
    pcr = g.input("pcr_enzyme", "Amplification Mix polymerase", MaterialKind.CATALYST)
    library = g.product("library", "dual-index Illumina-compatible PicoPLEX library",
                        MaterialKind.LIBRARY, state=PICOPLEX_LIBRARY, inferred=True,
                        source=source + " (complete library bases proprietary)")
    g.react("index_pcr", "Exponentially amplify and install the unique dual indexes",
            StepKind.PCR,
            ((preamp, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
             (i7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                PICOPLEX_LIBRARY,
                reason="The manual withholds quasi-random/common-handle and complete index-primer bases."))
    g.sequence(library, SequencingPlan(
        limitation="Illumina compatibility and dual indexes are published, but primer landing bases are not.",
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        read_lengths={"Index 1 (i7)": 8, "Index 2 (i5)": 8},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    g.declare_reviewed_scope((
        "cell", "lysis_reagent", "lyse", "quasi_random_primers",
        "displacement_polymerase", "preamp", "preamp_product", "i5_primer", "i7_primer",
        "pcr_enzyme", "index_pcr", "library"), source)
    return g.require_publishable()


def _sureplex_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("SurePlex with VeriSeq PGS library preparation", authority=source)
    cell = g.input("cell", "single cell", MaterialKind.BIOLOGICAL,
                   state=molecule("single cell", (role("genomic DNA"),)), primary=True)
    lysis = g.input("lysis_reagent", "SurePlex cell-extraction reagents",
                    MaterialKind.CONDITION)
    dna = g.product("released_dna", "released genomic DNA", MaterialKind.MOLECULAR,
                    state=molecule("released genomic DNA", (role("genomic DNA"),)))
    g.react("lyse", "Lyse the cell and release genomic DNA", StepKind.EXTRACTION,
            ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (dna,))
    primers = g.input("preamp_primers", "proprietary SurePlex preamplification primers",
                      MaterialKind.OLIGO, inferred=True,
                      source=source + " (primer sequences proprietary)")
    enzyme = g.input("preamp_enzyme", "SurePlex preamplification enzyme",
                     MaterialKind.CATALYST)
    wga = g.product("wga_product", "SurePlex whole-genome amplification product",
                    MaterialKind.MOLECULAR,
                    state=molecule("SurePlex WGA population", (
                        role("common preamplification handle"),
                        role("sampled genomic tract"), role("displaced extension products"))))
    g.react("wga", "Preamplify and exponentially amplify the single-cell genome",
            StepKind.EXTENSION,
            ((dna, InputRole.TEMPLATE), (primers, InputRole.STRUCTURAL),
             (enzyme, InputRole.CATALYST)), (wga,))
    adapters = g.input("tn5_adapters", "VeriSeq PGS loaded adapter transposome",
                       MaterialKind.ADAPTER)
    tn5 = g.input("tn5", "VeriSeq PGS transposase", MaterialKind.CATALYST)
    tagged = g.product("tagged", "tagmented SurePlex DNA", MaterialKind.MOLECULAR,
                       state=molecule("tagged SurePlex fragment", (
                           role("left mosaic end"), role("SurePlex-amplified genomic fragment"),
                           role("right mosaic end"))))
    waste = g.product("tagging_waste", "non-library tagmentation products", MaterialKind.WASTE)
    g.react("tagment", "Fragment and tag SurePlex-amplified DNA", StepKind.TAGMENTATION,
            ((wga, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
             (tn5, InputRole.CATALYST)), (tagged, waste))
    g.discard(waste, "lacks two productive heterologous VeriSeq/Nextera ends")
    p5 = g.input("p5_primer", "S503 or S504 i5/P5 index primer", MaterialKind.OLIGO)
    p7 = g.input("p7_primer", "N701–N712 i7/P7 index primer", MaterialKind.OLIGO)
    polymerase = g.input("pcr_enzyme", "Nextera PCR Master Mix", MaterialKind.CATALYST)
    library = g.product("library", "dual-index single-read VeriSeq PGS library",
                        MaterialKind.LIBRARY, state=SUREPLEX_LIBRARY)
    g.react("index_pcr", "Amplify the insert and add the selected indexes",
            StepKind.PCR,
            ((tagged, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
             (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(
                SUREPLEX_LIBRARY,
                reason="The guide identifies the commercial S5xx/N7xx primers by reagent name."))
    roles = ("Read 1", "Index 1 (i7)", "Index 2 (i5)")
    primers_used = tuple(p for p in SUREPLEX_PRIMERS if p.role in roles)
    g.sequence(library, SequencingPlan(
        SUREPLEX_LIBRARY, primers_used, roles,
        {"Index 1 (i7)": 8, "Index 2 (i5)": 8},
        {"Read 1": "include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex"}))
    g.declare_reviewed_scope((
        "cell", "lysis_reagent", "lyse", "preamp_primers", "preamp_enzyme", "wga",
        "wga_product", "tn5_adapters", "tn5", "tagment", "p5_primer", "p7_primer",
        "pcr_enzyme", "index_pcr", "library"), source)
    return g.require_publishable()


def get(key: str):
    if key == "olink-explore":
        source = ('<a href="https://doi.org/10.1016/j.mcpro.2021.100168">'
                  'Wik et al., <i>Molecular &amp; Cellular Proteomics</i> (2021)</a>.')
        spec = SimpleNamespace(
            TITLE="Olink Explore — proximity extension assay with NGS readout",
            NOTES="01_olink-explore.html", SOURCE=source,
            SUMMARY=("Matched antibody–oligo pairs co-bind a protein, hybridize in "
                     "proximity and template polymerase extension. Two PCR stages retain "
                     "the paired assay barcodes and add an inline sample index for counting."),
            CAVEAT=("Olink publishes the 148-bp element order and run layout but not the "
                    "probe, assay-barcode, hybridization-arm or index sequences. Those "
                    "regions are therefore structured roles, not invented bases."),
            TOPOLOGY_PANELS=ProximityExtensionTopology().panels(),
            oligos=lambda: (), FINAL_LIBRARY=OLINK_LIBRARY,
            FINAL_CAPTION=("Published 148-bp architecture: P5, Read 1 site, paired forward "
                           "and reverse assay barcodes separated by the copied probe junction, "
                           "inline sample index and P7."),
            SEQ_PRIMERS=OLINK_SEQ_PRIMERS, ACQUIRED_ROLES=("Read 1",),
            READ_LENGTHS={"Read 1": 66},
            READ_USES={"Read 1": "extract_and_demultiplex"},
            SHOW_READ_CONTENT=False,
            READ_SPECS=(sp.read_spec(
                "Read 1", 66,
                "forward assay barcode + reverse assay barcode + inline sample index",
                "extract_and_demultiplex"),),
            SEQUENCING_INTRO=("A standard Read 1 primer binds the published Rd1SP site. "
                              "The NovaSeq run acquires one 66-cycle read and no separate "
                              "index read; assay and sample identities are inline."))
        spec.protocol_graph = lambda: _olink_graph(source)
        return spec
    if key == "lampore":
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
        spec.protocol_graph = lambda: _lampore_graph(source)
        return spec
    if key == "ovation-rna-seq-v2":
        source = ('<a href="https://doi.org/10.1093/nar/gkr547">'
                  'Tariq et al., <i>Nucleic Acids Research</i> (2011)</a> and the '
                  'NuGEN Ovation RNA-Seq V2 user guide.')
        spec = SimpleNamespace(
            TITLE="Ovation RNA-Seq V2 — Ribo-SPIA whole-transcriptome amplification",
            NOTES="01_ovation-rna-seq-v2.html", SOURCE=source,
            SUMMARY=("Chimeric DNA/RNA primers initiate first- and second-strand cDNA. "
                     "RNase H repeatedly exposes SPIA priming sites, and a strand-displacing "
                     "polymerase repeatedly displaces cDNA copies before "
                     "conventional Illumina library construction."),
            CAVEAT=("NuGEN did not publish the chimeric primer bases. The defining paper "
                    "used the Ovation product with a separate Encore Rapid Illumina "
                    "library workflow; canonical TruSeq-compatible ends below are marked "
                    "inferred rather than attributed to the Ribo-SPIA kit itself."),
            oligos=lambda: (), FINAL_LIBRARY=OVATION_LIBRARY,
            FINAL_CAPTION=("Sequencing library made from sheared Ribo-SPIA cDNA. The "
                           "transcript-derived insert is source-established; exact commercial "
                           "adapter/index oligos were not disclosed in the defining paper."),
            SEQ_PRIMERS=OVATION_PRIMERS,
            ACQUIRED_ROLES=("Read 1", "Index 1 (i7)", "Read 2"),
            READ_LENGTHS={"Read 1": 50, "Index 1 (i7)": 8, "Read 2": 50},
            READ_USES={"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"},
            SEQUENCING_INTRO=("The study acquired paired 50-nt Illumina reads. Standard "
                              "primer placement is computed from the represented final shell."))
        spec.protocol_graph = lambda: _ovation_graph(source)
        return spec
    if key == "repli-g-single-cell":
        source = ('QIAGEN <a href="https://www.qiagen.com/en-US/resources/download/'
                  'KitHandbook/hb-1306-003-hb-replig-singlecell-0424-ww">'
                  'REPLI-g Advanced DNA Single Cell Kit Handbook</a> and its dedicated '
                  'Nextera NGS library-preparation protocol.')
        spec = SimpleNamespace(
            TITLE="REPLI-g Single Cell — phi29 MDA followed by Nextera library prep",
            NOTES="01_repli-g-single-cell.html", SOURCE=source,
            SUMMARY=("A single cell is lysed and its genomic DNA alkaline-denatured. "
                     "proprietary random hexamers prime processive phi29 multiple-displacement "
                     "amplification; the official downstream protocol tagmentates the amplified "
                     "DNA and completes dual-index Nextera ends by PCR."),
            CAVEAT=("Random-primer sequences and formulation are proprietary. MDA creates "
                    "a branched population rather than a single linear amplicon; only the "
                    "subsequent Nextera library molecule is drawn linearly."),
            oligos=lambda: (), FINAL_LIBRARY=REPLIG_LIBRARY,
            FINAL_CAPTION=("Dual-index Nextera library from one fragment sampled from the "
                           "branched MDA population. Adapter geometry is supplied by the "
                           "separate QIAGEN-authorized Nextera workflow."),
            SEQ_PRIMERS=REPLIG_PRIMERS,
            ACQUIRED_ROLES=tuple(p.role for p in REPLIG_PRIMERS),
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8,
                          "Index 2 (i5)": 8, "Read 2": 150},
            READ_USES={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                       "Index 2 (i5)": "demultiplex", "Read 2": "include"},
            SEQUENCING_INTRO=("The kit does not prescribe a read length; 2 × 150 with dual "
                              "8-cycle indexes is shown as a conventional Nextera run plan."))
        spec.protocol_graph = lambda: _replig_graph(source)
        return spec
    if key == "sureplex":
        source = ('Illumina <a href="https://support.illumina.com/content/dam/'
                  'illumina-support/documents/documentation/chemistry_documentation/'
                  'veriseq-pgs/sureplex-protocol-reference-guide-15053626-02.pdf">'
                  'SurePlex Summary Protocol</a> and '
                  '<a href="https://support.illumina.com/content/dam/illumina-support/'
                  'documents/documentation/chemistry_documentation/veriseq-pgs/'
                  'veriseq-pgs-library-prep-reference-guide-15052877-05.pdf">'
                  'VeriSeq PGS Library Prep Reference Guide</a>.')
        primers = tuple(p for p in SUREPLEX_PRIMERS if p.role != "Read 2")
        spec = SimpleNamespace(
            TITLE="SurePlex — single-cell WGA with VeriSeq PGS library preparation",
            NOTES="01_sureplex.html", SOURCE=source,
            SUMMARY=("SurePlex lyses a single cell and performs proprietary preamplification "
                     "followed by whole-genome amplification. The named VeriSeq PGS workflow "
                     "then tagmentates 1 ng of unpurified product and adds S5xx/N7xx indexes "
                     "by limited-cycle PCR for single-read MiSeq sequencing."),
            CAVEAT=("SurePlex preamplification-primer bases are proprietary. Its WGA is "
                    "therefore shown by established molecular roles; the subsequent "
                    "VeriSeq/Nextera library shell and sequencing-primer sites are explicit."),
            oligos=lambda: (), FINAL_LIBRARY=SUREPLEX_LIBRARY,
            FINAL_CAPTION=("Dual-index VeriSeq PGS library derived from SurePlex WGA DNA. "
                           "The official application acquires Read 1 and index reads, not Read 2."),
            SEQ_PRIMERS=primers,
            ACQUIRED_ROLES=("Read 1", "Index 1 (i7)", "Index 2 (i5)"),
            READ_LENGTHS={"Index 1 (i7)": 8, "Index 2 (i5)": 8},
            READ_USES={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                       "Index 2 (i5)": "demultiplex"},
            SEQUENCING_INTRO=("VeriSeq PGS specifies single-read acquisition with one or "
                              "two index reads; the sample-derived insert begins in Read 1."))
        spec.protocol_graph = lambda: _sureplex_graph(source)
        return spec
    if key == "picoplex-dna-seq":
        source = ('Takara Bio <a href="https://takara.co.kr/file/manual/pdf/'
                  'PicoPLEX-DNA-seq-Manual-QAM-106-002-updated%28RB3381%29.pdf">'
                  'PicoPLEX DNA-seq Kit Instruction Manual QAM-106-002</a>.')
        spec = SimpleNamespace(
            TITLE="PicoPLEX DNA-seq — displacement DOP-PCR with dual indexing",
            NOTES="01_picoplex-dna-seq.html", SOURCE=source,
            SUMMARY=("Single-cell DNA is lysed, repeatedly copied by proprietary "
                     "quasi-random priming and strand displacement, then exponentially "
                     "amplified with one of 48 unique dual-index pairs to make an "
                     "Illumina-compatible library in the same vessel."),
            CAVEAT=("The manual defines the three modules and 8-nt i5/i7 index plate but "
                    "withholds quasi-random primer, common-handle and complete adapter "
                    "sequences. Their roles are shown explicitly without assigning bases."),
            oligos=lambda: (), FINAL_LIBRARY=PICOPLEX_LIBRARY,
            FINAL_CAPTION=("Source-declared Illumina-compatible dual-index library. "
                           "Bracketed regions are proprietary; the manual does not expose "
                           "enough bases to compute exact sequencing-primer landing sites."),
            SEQ_PRIMERS=(),
            SEQUENCING_ENDING=("The finished dual-index library clusters on an Illumina "
                               "flow cell and uses the instrument's paired-read and index "
                               "sequencing primers; exact landing bases are undisclosed."),
            SEQUENCING_UNAVAILABLE=("The commercial common handles and full adapter "
                                    "sequences are proprietary."),
            READ_SPECS=(sp.read_spec("Read 1", None, "genomic insert from one end", "include"),
                        sp.read_spec("Index 1 (i7)", 8, "i7 sample index", "demultiplex"),
                        sp.read_spec("Index 2 (i5)", 8, "i5 sample index", "demultiplex"),
                        sp.read_spec("Read 2", None, "genomic insert from the opposite end", "include")))
        spec.protocol_graph = lambda: _picoplex_graph(source)
        return spec
    if key != "lamp-seq":
        raise KeyError(key)
    source = ('<a href="https://doi.org/10.1038/s41587-021-00966-9">'
              'Ludwig et al., <i>Nature Biotechnology</i> (2021)</a> and Supplementary Data 1.')
    spec = SimpleNamespace(
        TITLE="LAMP-Seq — barcoded RT-LAMP with Illumina readout",
        NOTES="01_lamp-seq.html", SOURCE=source,
        SUMMARY=("A sample-specific 10-nt barcode is embedded in the forward inner LAMP "
                 "primer. Saturated RT-LAMP reactions are pooled, strongly diluted, and "
                 "receive a second pair of 10-nt pool indexes during Illumina library PCR."),
        CAVEAT=("LAMP produces a population of branched, stem-looped concatemers rather "
                "than one linear molecule. The page therefore shows the source-established "
                "primer and PCR-selected regions without inventing a single universal "
                "ordering for every LAMP intermediate."),
        TOPOLOGY_PANELS=lamp_topology(
            "SARS-CoV-2 N-gene DNA", "10-nt FIP sample barcode").panels(),
        oligos=_lamp_oligos, FINAL_LIBRARY=LAMP_LIBRARY,
        FINAL_CAPTION=("Dual-index Illumina amplicon selected from the heterogeneous LAMP "
                       "product. Exact PCR landing sites and a representative published "
                       "FIP barcode are shown; variable concatemer-derived sequence remains "
                       "role-labelled."),
        SEQ_PRIMERS=LAMP_SEQ_PRIMERS,
        ACQUIRED_ROLES=("Read 1", "Index 1 (i7)", "Index 2 (i5)"),
        READ_LENGTHS={"Read 1": 100, "Index 1 (i7)": 8, "Index 2 (i5)": 8},
        READ_USES={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex"},
        SEQUENCING_INTRO=("The published MiSeq/iSeq run acquires 100 cycles of Read 1 and "
                          "8 cycles for each pool index; no Read 2 is acquired."))
    spec.protocol_graph = lambda: _lamp_graph(source)
    return spec
