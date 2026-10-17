"""PicoPLEX DNA-seq (Rubicon Genomics manual QAM-106-002) protocol declaration."""
from __future__ import annotations

from types import SimpleNamespace

from chemdraw import Construct, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
import seqprimers as sp


MANUAL = "PicoPLEX DNA-seq Instruction Manual QAM-106-002"

# Figure 4 (manual p. 7) and Figure 7 (Appendix 2, p. 17): rows A–H carry i5
# D501–D508, columns 1–6 carry i7 D701–D706; wells A7–H12 are blank.  The manual
# prints each 8-nt index as below but does not say which strand or read orientation
# the printed text represents, so these are not placed into a library sequence.
I5_INDEXES = {
    "D501": "TATAGCCT", "D502": "ATAGAGGC", "D503": "CCTATCCT", "D504": "GGCTCTGA",
    "D505": "AGGCGAAG", "D506": "TAATCTTA", "D507": "CAGGACGT", "D508": "GTACTGAC",
}
I7_INDEXES = {
    "D701": "ATTACTCG", "D702": "TCCGGAGA", "D703": "CGCTCATT",
    "D704": "GAGATTCC", "D705": "ATTCAGAA", "D706": "GAATTCGT",
}
# "Analysis of the Sequencing Reads" (p. 15): the first 11 cycles of each read
# contain quasi-random bases; the manual says to trim the initial 14 bases.
QUASI_RANDOM_CYCLES = 11
TRIM_BASES = 14

PICOPLEX_I5 = feature("picoplex_i5", "sample_index", "whitelist",
                      whitelist=f"{MANUAL} Figure 4/7 i5 D501–D508")
PICOPLEX_I7 = feature("picoplex_i7", "sample_index", "whitelist",
                      whitelist=f"{MANUAL} Figure 4/7 i7 D701–D706")
PICOPLEX_LIBRARY = Construct([
    role("Illumina-compatible P5-side end"),
    role("i5 sample index", tag="cbc", length_bp=8, molecular_feature=PICOPLEX_I5),
    role("read-one common handle"),
    role("quasi-random primer-derived bases", length_bp=QUASI_RANDOM_CYCLES),
    role("single-cell genomic insert"),
    role("quasi-random primer-derived bases (opposite end)", length_bp=QUASI_RANDOM_CYCLES),
    role("read-two common handle"),
    role("i7 sample index", tag="cbc", length_bp=8, molecular_feature=PICOPLEX_I7),
    role("Illumina-compatible P7-side end"),
], name="PicoPLEX DNA-seq dual-index library")


def _picoplex_graph(source: str) -> ProtocolGraph:
    g = ProtocolGraph("PicoPLEX DNA-seq", authority=source)
    cell = g.input("cell", "one to ten cells (≤2.5 µL PBS) or 6–60 pg purified DNA in Cell Extraction Buffer",
                   MaterialKind.BIOLOGICAL,
                   state=molecule("low-input sample", (role("genomic DNA"),)), primary=True)
    lysis = g.input("lysis_reagent",
                    "Cell Extraction Buffer + Extraction Enzyme Dilution Buffer/Cell Extraction Enzyme "
                    "master mix; 75 °C 10 min, 95 °C 4 min, 22 °C hold (10 µL)",
                    MaterialKind.CONDITION)
    dna = g.product("released_dna", "released genomic DNA", MaterialKind.MOLECULAR,
                    state=molecule("released genomic DNA", (role("genomic DNA"),)))
    g.react("lyse", "Lyse cells and release genomic DNA", StepKind.EXTRACTION,
            ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (dna,))
    quasi = g.input("quasi_random_primers",
                    "proprietary quasi-random primers in Pre-Amp Buffer",
                    MaterialKind.OLIGO, inferred=True,
                    source=source + " (primer sequences proprietary)")
    preamp_enzyme = g.input("preamp_enzyme",
                            "Pre-Amp Enzyme (identity and strand-displacement activity not stated)",
                            MaterialKind.CATALYST)
    preamp_cycling = g.input("preamp_cycling",
                             "same tube, 15 µL: 95 °C 2 min; 12 × (95 °C 15 s, 15 °C 50 s, "
                             "25 °C 40 s, 35 °C 30 s, 65 °C 40 s, 75 °C 40 s); 4 °C hold",
                             MaterialKind.CONDITION)
    preamp = g.product("preamp_product",
                       "pre-amplification library (Figure 1: accumulating hairpin library)",
                       MaterialKind.MOLECULAR,
                       state=molecule("PicoPLEX preamplification population", (
                           role("quasi-random primer common tail"),
                           role("quasi-random primer-derived bases"),
                           role("sampled genomic tract"),
                           role("complement of quasi-random primer and tail"))))
    g.react("preamp", "Multiple cycles of quasi-random priming and linear amplification",
            StepKind.EXTENSION,
            ((dna, InputRole.TEMPLATE), (quasi, InputRole.STRUCTURAL),
             (preamp_enzyme, InputRole.CATALYST), (preamp_cycling, InputRole.CONDITION)),
            (preamp,))
    i5 = g.input("i5_primer", "8-nt i5 indexed amplification primer (Dual Index Plate row, D501–D508)",
                 MaterialKind.OLIGO,
                 inferred=True, source=source + " (index printed; complete primer bases proprietary)")
    i7 = g.input("i7_primer", "8-nt i7 indexed amplification primer (Dual Index Plate column, D701–D706)",
                 MaterialKind.OLIGO,
                 inferred=True, source=source + " (index printed; complete primer bases proprietary)")
    pcr = g.input("pcr_enzyme", "Amplification Buffer + Amplification Enzyme", MaterialKind.CATALYST)
    pcr_cycling = g.input("pcr_cycling",
                          "same tube, 50 µL (optional EvaGreen/fluorescein): 95 °C 4 min; "
                          "4 × (95 °C 20 s, 63 °C 25 s, 72 °C 40 s); "
                          "7 × (95 °C 20 s, 72 °C 55 s); 4 °C hold",
                          MaterialKind.CONDITION)
    library = g.product("library", "dual-index Illumina-compatible PicoPLEX library",
                        MaterialKind.LIBRARY, state=PICOPLEX_LIBRARY, inferred=True,
                        source=source + " (complete library bases proprietary)")
    g.react("index_pcr", "Exponentially amplify and install the unique dual indexes",
            StepKind.PCR,
            ((preamp, InputRole.TEMPLATE), (i5, InputRole.STRUCTURAL),
             (i7, InputRole.STRUCTURAL), (pcr, InputRole.CATALYST),
             (pcr_cycling, InputRole.CONDITION)), (library,),
            geometry=PCRGeometry.unavailable(
                PICOPLEX_LIBRARY,
                reason="The manual withholds quasi-random/common-handle and complete index-primer bases."))
    pool = g.product("pooled_libraries", "equal-volume pool of uniquely dual-indexed libraries",
                     MaterialKind.LIBRARY, state=PICOPLEX_LIBRARY, inferred=True)
    others = g.input("other_libraries",
                     "other PicoPLEX libraries, each with a different dual-index pair (up to 48 in all)",
                     MaterialKind.LIBRARY, state=PICOPLEX_LIBRARY, inferred=True)
    g.react("pool", "Pool 10 µL of each indexed library", StepKind.POOL,
            ((library, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)), (pool,),
            changes_molecule=False)
    ampure = g.input("ampure", "AMPure XP 1:1 (100 µL + 100 µL pool), two 80% ethanol washes, "
                     "50 µL TE elution", MaterialKind.SELECTION)
    clean = g.product("purified_pool", "AMPure-purified pooled PicoPLEX library",
                      MaterialKind.LIBRARY, state=PICOPLEX_LIBRARY, inferred=True)
    waste = g.product("ampure_supernatant", "supernatant and washes with unincorporated primers",
                      MaterialKind.WASTE)
    g.react("ampure_cleanup", "Purify the pool on AMPure XP beads", StepKind.CLEANUP,
            ((pool, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (clean, waste),
            changes_molecule=False)
    g.discard(waste, "unincorporated primers and reagents removed with the bead supernatant")
    g.sequence(clean, SequencingPlan(
        limitation="Illumina compatibility and dual indexes are published, but primer landing bases are not.",
        acquired_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"),
        read_lengths={"Index 1 (i7)": 8, "Index 2 (i5)": 8},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    g.annotate(SourceDiscrepancy(
        "purified_pool", "Quasi-random read prefix: 11 cycles stated, 14 bases trimmed",
        (SourceClaim(f"{MANUAL}, Analysis of the Sequencing Reads (p. 15)",
                     f"the first {QUASI_RANDOM_CYCLES} cycles of each read contain quasi-random bases"),
         SourceClaim(f"{MANUAL}, Analysis of the Sequencing Reads (p. 15)",
                     f"trim the initial {TRIM_BASES} bases from each read or begin at base 15")),
        f"The library shows {QUASI_RANDOM_CYCLES}-nt quasi-random primer-derived regions at both "
        f"insert ends; the read table follows the manual's {TRIM_BASES}-base trim.",
        "Both numbers appear in the same paragraph of the manual.",
        "The manual does not explain the extra three trimmed bases (for example a fixed primer "
        "segment or a safety margin)."))
    g.annotate(SourceDiscrepancy(
        "cell", "Purified-DNA input range is worded two ways",
        (SourceClaim(f"{MANUAL}, Starting Material: Genomic DNA (p. 6)",
                     "small amounts (less than 6 pg to 60 pg) of purified genomic DNA"),
         SourceClaim(f"{MANUAL}, Cell Lysis Protocol step 1 (p. 8)",
                     "cells (1 to 10) or DNA (6-60 pg)")),
        "The model uses the protocol-step range 6–60 pg; the positive control uses 15 pg.",
        "Both statements are in the same manual revision.",
        "Whether inputs below 6 pg are supported is not stated elsewhere in the manual."))
    g.declare_reviewed_scope((
        "cell", "lysis_reagent", "lyse", "quasi_random_primers", "preamp_enzyme",
        "preamp_cycling", "preamp", "preamp_product", "i5_primer", "i7_primer",
        "pcr_enzyme", "pcr_cycling", "index_pcr", "library", "pool", "pooled_libraries",
        "other_libraries", "ampure", "ampure_cleanup", "purified_pool"), source)
    return g.require_publishable()


def _index_rows() -> tuple[str, ...]:
    i5 = " · ".join(f"{k} {v}" for k, v in I5_INDEXES.items())
    i7 = " · ".join(f"{k} {v}" for k, v in I7_INDEXES.items())
    return (f"<b>i5 indexes (plate rows A–H):</b> {i5}",
            f"<b>i7 indexes (plate columns 1–6):</b> {i7}",
            "<i>As printed in manual Figures 4 and 7. The manual does not say which strand "
            "these sequences are written on, and it gives no full primer sequences.</i>")


source = ('Rubicon Genomics (now Takara Bio) <a href="https://takara.co.kr/file/manual/pdf/'
          'PicoPLEX-DNA-seq-Manual-QAM-106-002-updated%28RB3381%29.pdf">'
          'PicoPLEX DNA-seq Instruction Manual QAM-106-002</a> (Cat# R300381).')
spec = SimpleNamespace(
    TITLE="PicoPLEX DNA-seq — quasi-random preamplification with dual indexing",
    NOTES="01_picoplex-dna-seq.html", SOURCE=source,
    SUMMARY=("Single-cell DNA is lysed, copied through 12 cycles of proprietary "
             "quasi-random priming and linear amplification, then exponentially "
             "amplified in the same tube with one of 48 dual-index pairs. The "
             "indexed libraries are pooled and purified once on AMPure XP beads."),
    CAVEAT=("The manual prints the 8-nt i5/i7 index sequences but withholds the "
            "quasi-random primer, common-handle and complete adapter sequences. "
            "Those regions are shown by role, without bases."),
    oligos=_index_rows, FINAL_LIBRARY=PICOPLEX_LIBRARY,
    FINAL_CAPTION=("Illumina-compatible dual-index library as the manual describes it. "
                   "Bracketed regions are proprietary. The manual does not give enough "
                   "bases to compute where the sequencing primers land. Libraries are "
                   "~300–1000 bp in total (~200–900 bp inserts)."),
    SEQ_PRIMERS=(),
    SEQUENCING_ENDING=("The purified pool is sequenced with standard Illumina reagents "
                       "for multiplexed libraries (starting point: 16 pM on MiSeq v3 "
                       "with at least 5% PhiX). The exact primer landing bases are not "
                       "disclosed, and the manual sets no read lengths."),
    SEQUENCING_UNAVAILABLE=("The commercial common handles and full adapter "
                            "sequences are proprietary."),
    READ_SPECS=(sp.read_spec("Read 1", None,
                             f"first {QUASI_RANDOM_CYCLES} cycles quasi-random primer bases "
                             f"(manual: trim {TRIM_BASES}), then genomic insert", "include"),
                sp.read_spec("Index 1 (i7)", 8, "i7 sample index (D701–D706)", "demultiplex"),
                sp.read_spec("Index 2 (i5)", 8, "i5 sample index (D501–D508)", "demultiplex"),
                sp.read_spec("Read 2", None,
                             f"first {QUASI_RANDOM_CYCLES} cycles quasi-random primer bases "
                             f"(manual: trim {TRIM_BASES}), then genomic insert from the opposite end",
                             "include")))


def protocol_graph() -> ProtocolGraph:
    return _picoplex_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
