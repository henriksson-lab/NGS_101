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


def protocol_graph() -> ProtocolGraph:
    return _picoplex_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
