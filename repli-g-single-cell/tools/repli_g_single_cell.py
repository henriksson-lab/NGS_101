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


def protocol_graph() -> ProtocolGraph:
    return _replig_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
