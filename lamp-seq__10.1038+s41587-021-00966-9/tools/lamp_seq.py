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
# Supplementary Data 1 "Pool barcoding primers" (row Pool_fwd_1/Pool_rev_1) and the
# Methods two-step PCR print 23- and 24-nt partial sequencing-primer tails, not full
# TruSeq arms: the first PCR supplies the remaining Read 1 / Read 2 bases.
LAMP_POOL_FWD_TAIL = il.TRUSEQ_READ1[:23]   # ACACTCTTTCCCTACACGACGCT
LAMP_POOL_REV_TAIL = il.TRUSEQ_READ2[:24]   # GTGACTGGAGTTCAGACGTGTGCT
LAMP_POOL_FWD = il.P5 + LAMP_POOL_I5 + LAMP_POOL_FWD_TAIL
LAMP_POOL_REV = il.P7 + LAMP_POOL_I7 + LAMP_POOL_REV_TAIL
# The 5' G of the pool reverse tail completes canonical GTGACT on the TGACT first-PCR arm.
if LAMP_POOL_REV_TAIL[1:] != LAMP_PCR_REV_ARM[:len(LAMP_POOL_REV_TAIL) - 1]:
    raise ValueError("pool reverse tail does not overlap the first-PCR reverse arm")

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
    yield oligo("Pool_fwd_1", [Segment("P5", il.P5, "p5"),
                               Segment("i5 pool index", LAMP_POOL_I5, "cbc"),
                               Segment("partial Read 1 tail", LAMP_POOL_FWD_TAIL, "r1")])
    yield oligo("Pool_rev_1", [Segment("P7", il.P7, "p7"),
                               Segment("i7 pool index", LAMP_POOL_I7, "cbc"),
                               Segment("partial Read 2 tail", LAMP_POOL_REV_TAIL, "r2")])


def _lamp_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("LAMP-Seq", authority=source)
    sample = g.input("sample", "oropharyngeal swab in QuickExtract (+2 ng/µl pX330 carrier DNA)",
                     MaterialKind.BIOLOGICAL,
                     state=molecule("swab lysate", (role("viral RNA"),)), primary=True)
    heat = g.input("heat", "≥10 min room temperature, then 95 °C 5 min", MaterialKind.CONDITION)
    lysate = g.product("heated_lysate", "heat-inactivated RNA-containing lysate",
                       MaterialKind.MOLECULAR, state=sample.state)
    g.react("inactivate", "Heat-inactivate the swab lysate (then 30 min carbon/ion-exchange beads and HCl pH adjustment, not drawn separately)", StepKind.THERMAL_TREATMENT,
            ((sample, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (lysate,))

    fip = g.input("lamp_fip", "barcoded C-FIP (representative Tier-1 barcode AGCAAACTAC)",
                  MaterialKind.OLIGO, state=LAMP_FIP)
    lnas = g.input("lamp_outer", "C-F3-LNA and C-B3-LNA outer primers (LNA residues not drawn)",
                   MaterialKind.OLIGO)
    others = g.input("lamp_bip_loops", "C-BIP, C-LF and C-LB",
                     MaterialKind.OLIGO, state=(LAMP_BIP, LAMP_LF, LAMP_LB))
    bst = g.input("bst", "NEB LAMP master mix (E1700L) plus Bst 3.0 (M0374L), 65 °C 80 min",
                  MaterialKind.CATALYST)
    lamp = g.product("lamp_product", "barcoded branched RT-LAMP products",
                     MaterialKind.MOLECULAR, state=lamp_product("SARS-CoV-2 N-gene DNA"))
    g.react("rt_lamp", "Reverse-transcribe and amplify by loop-mediated strand displacement",
            StepKind.EXTENSION,
            ((lysate, InputRole.TEMPLATE), (fip, InputRole.STRUCTURAL),
             (lnas, InputRole.STRUCTURAL), (others, InputRole.STRUCTURAL),
             (bst, InputRole.CATALYST)), (lamp,))

    dilution = g.input("dilution", "1:100,000 pre-PCR dilution (1:1,000,000 in first PCR)",
                       MaterialKind.CONDITION)
    pooled = g.product("pooled_lamp", "pooled and diluted barcoded LAMP products",
                       MaterialKind.MOLECULAR, state=lamp.state)
    g.react("pool_dilute", "Pool barcoded reactions and dilute before PCR",
            StepKind.GENERIC,
            ((lamp, InputRole.TEMPLATE), (dilution, InputRole.CONDITION)), (pooled,),
            changes_molecule=False)

    pcr_f = g.input("pcr_f", "PCR-N-fwd-10 (Methods: PCR-C-fwd)", MaterialKind.OLIGO,
                    state=LAMP_PCR_FWD)
    pcr_r = g.input("pcr_r", "PCR-N-rev-11 (Methods: PCR-C-rev)", MaterialKind.OLIGO,
                    state=LAMP_PCR_REV)
    polymerase = g.input("pcr_polymerase", "NEBNext 2× Master Mix", MaterialKind.CATALYST)
    preamp = g.product("preamp", "first-stage LAMP-Seq amplicon",
                       MaterialKind.MOLECULAR, state=LAMP_PREAMP)
    geometry1 = PCRGeometry.exact(
        LAMP_PCR_TEMPLATE,
        PCRPrimer("PCR-N-fwd-10", LAMP_PCR_FWD, anneal_nt=len(LAMP_PCR_FWD_SITE)),
        PCRPrimer("PCR-N-rev-11", LAMP_PCR_REV, anneal_nt=len(LAMP_PCR_REV_SITE)),
        statement="Ludwig et al. Supplementary Data 1 publishes both target-tailed PCR primers.")
    g.react("target_pcr", "Select one N-gene interval by 18-cycle PCR",
            StepKind.PCR,
            ((pooled, InputRole.TEMPLATE), (pcr_f, InputRole.STRUCTURAL),
             (pcr_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (preamp,), geometry=geometry1)

    pool_f = g.input("pool_f", "Pool_fwd_1 representative indexed primer", MaterialKind.OLIGO,
                     state=LAMP_POOL_FWD)
    pool_r = g.input("pool_r", "Pool_rev_1 representative indexed primer", MaterialKind.OLIGO,
                     state=LAMP_POOL_REV)
    indexed = g.product("indexed_pcr", "pool-indexed LAMP-Seq PCR product", MaterialKind.LIBRARY,
                        state=LAMP_LIBRARY)
    geometry2 = PCRGeometry.exact(
        LAMP_PREAMP,
        PCRPrimer("Pool_fwd_1", LAMP_POOL_FWD, anneal_nt=len(LAMP_POOL_FWD_TAIL)),
        PCRPrimer("Pool_rev_1", LAMP_POOL_REV, anneal_nt=len(LAMP_POOL_REV_TAIL) - 1),
        statement="Supplementary Data 1 publishes 192 pairs of 10-nt pool-index primers; pair 1 is shown.")
    g.react("index_pcr", "Install pool indexes by a second 18-cycle PCR",
            StepKind.PCR,
            ((preamp, InputRole.TEMPLATE), (pool_f, InputRole.STRUCTURAL),
             (pool_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (indexed,), geometry=geometry2)
    gel = g.input("gel_columns", "2% agarose E-gel band excision, two Qiagen silica spin columns",
                  MaterialKind.SELECTION)
    final = g.product("library", "dual-index LAMP-Seq library", MaterialKind.LIBRARY,
                      state=LAMP_LIBRARY)
    rejected = g.product("gel_rejects", "off-size gel material, column flow-through and washes",
                         MaterialKind.WASTE)
    g.react("gel_purify", "Pool PCR products, cut the library band and column-purify twice",
            StepKind.SELECTION, ((indexed, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
            (final, rejected))
    g.discard(rejected, "outside the excised library band or not retained by the columns")
    g.sequence(final, SequencingPlan(
        LAMP_LIBRARY, tuple(LAMP_SEQ_PRIMERS),
        ("Read 1", "Index 1 (i7)", "Index 2 (i5)"),
        {"Read 1": 100, "Index 1 (i7)": 8, "Index 2 (i5)": 8},
        {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
         "Index 2 (i5)": "demultiplex"}))
    g.declare_reviewed_scope((
        "sample", "heat", "inactivate", "lamp_fip", "lamp_outer", "lamp_bip_loops", "bst",
        "rt_lamp", "dilution", "pool_dilute", "pcr_f", "pcr_r", "pcr_polymerase",
        "target_pcr", "pool_f", "pool_r", "index_pcr", "indexed_pcr", "gel_columns",
        "gel_purify", "library"), source)
    return g.require_publishable()




source = ('<a href="https://doi.org/10.1038/s41587-021-00966-9">'
          'Ludwig et al., <i>Nature Biotechnology</i> (2021)</a> and Supplementary Data 1.')
spec = SimpleNamespace(
    TITLE="LAMP-Seq — barcoded RT-LAMP with Illumina readout",
    NOTES="01_lamp-seq.html", SOURCE=source,
    SUMMARY=("This page follows the paper's two-PCR LAMP-Seq workflow. A sample-specific "
             "10-nt barcode is embedded in the forward inner LAMP primer. Saturated "
             "reactions are pooled, diluted 1:100,000, then receive two 10-nt pool "
             "indexes in the second PCR."),
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
    SEQUENCING_INTRO=("The MiSeq Nano and iSeq runs printed in the open-source protocol acquire "
                      "100 cycles of Read 1 and 8 cycles for each pool index, so only the first "
                      "8 nt of each 10-nt pool index are read; no Read 2 is acquired. The "
                      "commercial two-step branch drawn here names MiSeq or iSeq without "
                      "separate cycle numbers."))


def protocol_graph() -> ProtocolGraph:
    return _lamp_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
