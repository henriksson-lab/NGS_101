"""Source-owned reaction graph for scCOOL-seq."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scnome-seq__10.7554+eLife.23203" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nomecool as N
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from pcr import PCRGeometry


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/cr.2017.82"
    graph = ProtocolGraph(N.SCCOOL, authority=source)
    converted = N.mark_and_convert_graph(graph)

    primer1 = graph.input("random_primer1", "biotinylated random primer 1",
                          MaterialKind.OLIGO, state=tuple(N.primer1_segments()))
    polymerase = graph.input("extension_polymerase", "Klenow fragment (3-prime to 5-prime exo-minus)",
                             MaterialKind.CATALYST)
    first = graph.product("first_strands", "biotinylated first-strand products",
                          MaterialKind.MOLECULAR)
    graph.react("first_priming", "Random-prime the first strand", StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (primer1, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (first,))

    exoi = graph.input("exonuclease_i", "Exonuclease I", MaterialKind.CATALYST)
    treated = graph.product("exoi_treated", "first-strand products after primer digestion",
                            MaterialKind.MOLECULAR)
    graph.react("primer_digestion", "Digest residual primers with Exonuclease I",
                StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (exoi, InputRole.CATALYST)), (treated,))
    ampure = graph.input("ampure", "Agencourt AMPure XP beads", MaterialKind.SELECTION)
    cleaned = graph.product("cleaned_first", "purified biotinylated first strands",
                            MaterialKind.MOLECULAR)
    cleanup_waste = graph.product("cleanup_waste", "material not recovered in first AMPure purification",
                                  MaterialKind.WASTE)
    graph.react("first_cleanup", "Purify first-strand products with AMPure XP",
                StepKind.CLEANUP,
                ((treated, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (cleaned, cleanup_waste))
    graph.discard(cleanup_waste, "removed during AMPure recovery; molecular composition not assigned")

    beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_strands", "bead-bound first strands",
                             MaterialKind.MOLECULAR)
    waste = graph.product("uncaptured_material", "uncaptured material", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated first strands and wash",
                StepKind.AFFINITY_CAPTURE,
                ((cleaned, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, waste))
    graph.discard(waste, "unbound material and removed converted templates are washed away")

    primer2 = graph.input("random_primer2", "random primer 2", MaterialKind.OLIGO,
                          state=tuple(N.primer2_segments()))
    template = graph.product("dual_handle_products", "dual-handle second-strand products",
                             MaterialKind.MOLECULAR)
    graph.react("second_priming", "Random-prime the captured second strand",
                StepKind.EXTENSION,
                ((captured, InputRole.TEMPLATE), (primer2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (template,))

    forward = graph.input("forward_pe1", "Forward PE1.0 PCR primer", MaterialKind.OLIGO,
                          inferred=True, source=source + " (sequence inferred from reagent name)")
    reverse = graph.input("indexed_reverse", "NEB indexed reverse PCR primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (exact primer sequence unavailable)")
    reviewed = ("genome", "mcvi", "gpc_marking", "bisulfite", "bisulfite_conversion",
                "random_primer1", "extension_polymerase", "first_priming",
                "streptavidin_beads", "capture", "random_primer2", "second_priming",
                "forward_pe1", "indexed_reverse")
    polymerase_pcr = graph.input("pcr_polymerase", "KAPA HiFi HotStart DNA polymerase",
                                MaterialKind.CATALYST)
    final = N.cool_library()
    amplified = graph.product("amplified", "amplified scCOOL-seq library",
                              MaterialKind.MOLECULAR, state=final, inferred=True)
    graph.react("indexed_pcr", "Amplify on the beads: 13 PCR cycles in main-article route",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase_pcr, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(final,
                    reason="Exact indexed PCR primer identities and pre-PCR molecular states remain unresolved."))
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            inferred=True)
    final_waste = graph.product("final_waste", "material not recovered by final AMPure purifications",
                                MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify amplified libraries twice with AMPure XP",
                StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (library, final_waste))
    graph.discard(final_waste, "not recovered in final library eluate; exact composition unrepresented")
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(N.SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 150, "Read 2": 150},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"},
        limitation="PCR primer identity and index width remain inferred; two 150-base reads are source-stated."))
    graph.annotate(SourceDiscrepancy(
        "indexed_pcr", "Published PCR settings differ",
        (SourceClaim("Main article, Methods: Single-cell COOL-seq library preparation (PMC5539349 JATS153)",
                     "13 cycles; Forward PE1.0 and Illumina reverse indexed primer"),
         SourceClaim("Supplementary Data S1, p3 PCR amplification steps1–2 (x22.pdf.txt114–120)",
                     "12 cycles; NEB Universal PCR Primer and PCR Index Primer")),
        "The graph follows the main article's 13-cycle route with named primers; it does not assign exact primer sequences.",
        "Both settings are printed in the published article/supplement.",
        "The sources do not resolve whether these are alternate implementations; no universal corrected cycle count is asserted."))
    graph.annotate(SourceDiscrepancy(
        "gpc_marking", "Published methyltransferase heat-inactivation durations differ",
        (SourceClaim("Main article, Methods: Preparation and in vitro methylation (JATS145)",
                     "65 C for 20 min"),
         SourceClaim("Supplementary Data S1, p1 Day1 step3 (x22.pdf.txt35)",
                     "65 C for 25 min")),
        "The compact graph leaves the heat-inactivation duration unassigned.",
        "Both published originals specify 65 C but give different durations.",
        "No corrected duration is established by this capped pass."))
    graph.declare_reviewed_scope((*reviewed, "exonuclease_i", "primer_digestion", "ampure",
                                 "first_cleanup", "pcr_polymerase", "indexed_pcr",
                                 "final_cleanup", "library"), source)
    return graph.require_publishable()
