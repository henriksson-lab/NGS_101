"""Construct model for single-cell reduced-representation bisulfite sequencing.

The wet-lab order is from Guo et al. (Genome Research 2013). The paper identifies
standard premethylated indexed Illumina adapters; their single-index sequences are
published in Illumina's authoritative adapter-sequence document, which `lib/illumina.py`
carries, and `lib/rrbs.py` assembles every drawn base from those constants. Nothing here
derives from the scg_lib_structs drawn page: the paper names neither library PCR primer,
so both are sequence-less graph inputs and the PCR geometry stays unavailable.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from rrbs import (CELL_INDEX, INDEXED_FIXED_5, INDEX_NT, MSPI_OVERHANG,
                  MSPI_SITE, SEQ_PRIMERS, UNIVERSAL, a_tailed_fragment_scene,
                  adapter_scene, digested_fragment_scene, final_library,
                  indexed_segments, universal_segments)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1101/gr.161679.113"
    graph = ProtocolGraph("scRRBS", authority=source)
    genome = graph.input("genome", "single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    mspi = graph.input("mspi", "MspI restriction enzyme", MaterialKind.CATALYST)
    digested = graph.product("mspi_fragments", "MspI-digested genomic fragments",
                             MaterialKind.MOLECULAR, state=digested_fragment_scene())
    graph.react("mspi_digest", "Digest genomic DNA at C^CGG",
                StepKind.DIGESTION,
                ((genome, InputRole.TEMPLATE), (mspi, InputRole.CATALYST)), (digested,))
    klenow = graph.input("klenow", "Klenow exo- end-repair/dA-tailing enzyme",
                         MaterialKind.CATALYST)
    prepared = graph.product("a_tailed_fragments", "filled and dA-tailed MspI fragments",
                             MaterialKind.MOLECULAR, state=a_tailed_fragment_scene())
    graph.react("end_prepare", "Fill recessed ends and add one 3'-A",
                StepKind.END_REPAIR,
                ((digested, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)),
                (prepared,))
    adapters = graph.input("methylated_adapters", "canonical premethylated indexed TruSeq adapter illustration",
                           MaterialKind.ADAPTER, state=adapter_scene(), inferred=True)
    ligase = graph.input("ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_fragments", "adapter-ligated genomic fragments",
                            MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "INFERRED — ligate a compatible premethylated indexed Y-adapter",
                StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining paper names standard premethylated indexed Illumina adapters, but does not print their ligation-facing oligos or exact end geometry."))
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_library", "bisulfite-converted ligation products",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert adapter-ligated molecules",
                StepKind.CHEMICAL_CONVERSION,
                ((ligated, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    column = graph.input("zymo_column", "Zymo-Spin column with tRNA carrier",
                         MaterialKind.SELECTION)
    purified = graph.product("purified_converted", "desulfonated and column-purified converted DNA",
                             MaterialKind.MOLECULAR)
    column_waste = graph.product("column_waste", "unrecovered material from column purification",
                                 MaterialKind.WASTE)
    graph.react("column_purification", "On-column desulfonation and purification with tRNA carrier",
                StepKind.CLEANUP,
                ((converted, InputRole.TEMPLATE), (column, InputRole.SELECTION)),
                (purified, column_waste))
    graph.discard(column_waste, "not recovered in the column eluate; composition not represented")
    qp1 = graph.input("p5_arm_pcr_primer", "P5-arm library PCR primer, identity unresolved",
                      MaterialKind.OLIGO)
    qp2 = graph.input("p7_arm_pcr_primer", "P7-arm library PCR primer, identity unresolved",
                      MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "PfuTurbo Cx polymerase (first PCR)",
                             MaterialKind.CATALYST)
    final = final_library()
    first = graph.product("first_pcr_product", "first-round PCR product",
                          MaterialKind.MOLECULAR)
    reason = ("The inspected 2013 Methods passage does not assign library primer identities; "
              "linked originals remain uninspected. Pre-PCR molecular states are not reconstructed "
              "and the illustrated canonical Illumina arms remain inferred for scRRBS.")
    graph.react("library_pcr", "First PCR: 25 cycles with PfuTurbo Cx",
                StepKind.PCR,
                ((purified, InputRole.TEMPLATE), (qp1, InputRole.STRUCTURAL),
                 (qp2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (first,), geometry=PCRGeometry.unavailable(final, reason=reason))
    phusion = graph.input("phusion", "Phusion high-fidelity polymerase (second PCR)",
                         MaterialKind.CATALYST)
    amplified = graph.product("amplified", "second-round PCR product",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("second_pcr", "Second PCR: 22 cycles with Phusion",
                StepKind.PCR,
                ((first, InputRole.TEMPLATE), (qp1, InputRole.STRUCTURAL),
                 (qp2, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(final, reason=reason))
    gel = graph.input("size_gel", "12% native polyacrylamide TBE gel, 200–500 bp window",
                      MaterialKind.SELECTION)
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    offsize = graph.product("offsize", "products outside the recovered gel-size window",
                            MaterialKind.WASTE)
    graph.react("size_selection", "Resolve on native PAGE and recover 200–500 bp products",
                StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (library, offsize))
    graph.discard(offsize, "not recovered from the selected gel interval")
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Index 1 (i7)": INDEX_NT},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"},
        limitation="Paired-end HiSeq 2000 is source-stated; the illustrated six-base index design and its acquisition remain unverified for the study."))
    graph.declare_reviewed_scope((
        "genome", "mspi", "mspi_digest", "klenow", "end_prepare",
        "methylated_adapters", "ligase", "adapter_ligation", "bisulfite",
        "bisulfite_conversion", "zymo_column", "column_purification",
        "p5_arm_pcr_primer", "p7_arm_pcr_primer", "pcr_polymerase", "library_pcr",
        "phusion", "second_pcr", "size_gel", "size_selection", "library"), source)
    return graph.require_publishable()
