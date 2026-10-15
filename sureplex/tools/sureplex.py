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



SUREPLEX_LIBRARY, SUREPLEX_PRIMERS = nextera_library([
    role("fragment of SurePlex-amplified genomic DNA")
], "SurePlex / VeriSeq PGS library")



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


def protocol_graph() -> ProtocolGraph:
    return _sureplex_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
