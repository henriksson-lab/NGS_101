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
], "REPLI-g / Nextera library", inferred_adapters=True)


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
    random = g.input("random_primers", "exonuclease-resistant MDA primers (current product page: random hexamers)",
                     MaterialKind.OLIGO, inferred=True,
                     source=source + " (exact primer bases and modifications unresolved)")
    phi29 = g.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    mda = g.product("mda_product", "branched multiple-displacement amplification product",
                    MaterialKind.MOLECULAR, state=mda_product())
    g.react("mda", "Perform processive random-primed phi29 strand-displacement amplification",
            StepKind.EXTENSION,
            ((dna, InputRole.TEMPLATE), (random, InputRole.STRUCTURAL),
             (phi29, InputRole.CATALYST)), (mda,))
    blank = g.input("no_template_control", "no-template control (parallel tube)",
                    MaterialKind.CONDITION)
    dimer_seed = g.product("primer_dimer_seed", "annealed MDA-primer dimer in blank control",
                           MaterialKind.MOLECULAR,
                           state=molecule("primer-dimer seed", (role("intermolecular primer pairing"),)))
    g.react("control_primer_pair", "Allow MDA primers to pair in the no-template control",
            StepKind.HYBRIDIZATION,
            ((random, InputRole.TEMPLATE), (random, InputRole.STRUCTURAL),
             (blank, InputRole.CONDITION)), (dimer_seed,))
    primer_dimer = g.product("primer_dimer_dna", "high-molecular-weight primer-dimer DNA in no-template control",
                             MaterialKind.WASTE)
    g.react("control_primer_dimer", "Extend primer dimers in the no-template control",
            StepKind.EXTENSION,
            ((dimer_seed, InputRole.TEMPLATE), (random, InputRole.STRUCTURAL),
             (phi29, InputRole.CATALYST)), (primer_dimer,))
    g.discard(primer_dimer, "negative-control artifact, not sampled for the library")
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
    g.react("index_pcr", "INFERRED — complete flow-cell ends and add dual indexes by PCR",
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
        "mda", "mda_product", "no_template_control", "control_primer_pair",
        "primer_dimer_seed", "control_primer_dimer",
        "primer_dimer_dna",
        "tn5_adapters", "tn5", "tagment", "p5_primer", "p7_primer", "pcr_enzyme",
        "index_pcr", "library"), source)
    return g.require_publishable()




source = ('QIAGEN <a href="https://www.ys-bio.cn/app/wp-content/uploads/2021/12/'
          'EN-REPLI-g-Single-Cell-Handbook.pdf">'
          'REPLI-g Single Cell Kit Handbook (October 2012; supplier mirror)</a> and the '
          '<a href="https://www.qiagen.com/en-US/resources/download/Protocols/'
          'rg35-nextera-ngs-library-preparation-from-dna-cdna-amplified-with-repli-g-kits">'
          'August 2014 RG35 Nextera protocol</a> for that named kit.')
spec = SimpleNamespace(
    TITLE="REPLI-g Single Cell — phi29 MDA followed by Nextera library prep",
    NOTES="01_repli-g-single-cell.html", SOURCE=source,
    SUMMARY=("A single cell is lysed and its genomic DNA alkaline-denatured. "
             "proprietary random hexamers prime processive phi29 multiple-displacement "
             "amplification; the official downstream protocol tagmentates the amplified "
             "DNA and completes dual-index Nextera ends by PCR."),
    CAVEAT=("The inspected 2012 handbook does not print the primer bases or "
            "exonuclease-resistant chemistry; exact kit-oligo composition remains open. "
            "The current 2024 handbook and full August 2014 RG35 PDF still need direct review. "
            "The August 2014 RG35 indexed original names this REPLI-g Single Cell kit; "
            "the November 2018 revision instead names Advanced DNA Single Cell. "
            "The handbook also reports primer-dimer DNA in a parallel no-template control. "
            "Precise Nextera adapter arms and index choices are inferred until the cited "
            "October 2012 Illumina guide is inspected. MDA creates "
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
