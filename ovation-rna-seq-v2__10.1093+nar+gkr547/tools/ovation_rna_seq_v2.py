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


def protocol_graph() -> ProtocolGraph:
    return _ovation_graph(source)

spec.protocol_graph = protocol_graph
globals().update(vars(spec))
