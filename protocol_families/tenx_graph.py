"""Shared typed provenance for Chromium single-cell RNA library construction."""
from __future__ import annotations

from protocol_graph import (Disposition, InputRole, MaterialKind, StepKind, Terminal)


def partition_material(graph, material, reviewed, *, label="cells or nuclei",
                       unresolved_reason=None):
    """Partition one Chromium input and account for non-resolvable GEM products.

    By default those products are a physical discard.  A caller whose source pools
    every GEM without physical removal passes ``unresolved_reason`` (its own
    source-qualified explanation); the products then stay molecular and end at an
    UNRESOLVED terminal instead of being drawn as discarded waste.
    """
    gems = graph.input("gems", "Chromium GEM partitioning condition",
                       MaterialKind.CONDITION)
    occupied = graph.product("occupied", f"occupied GEMs containing {label}",
                             MaterialKind.MOLECULAR, state=material.state)
    failed = graph.product("failed", "empty, multiplet or unresolved GEM products",
                           MaterialKind.WASTE if unresolved_reason is None
                           else MaterialKind.MOLECULAR)
    graph.react("partition", f"Partition {label} with barcoded gel beads",
                StepKind.PARTITION,
                ((material, InputRole.TEMPLATE), (gems, InputRole.CONDITION)),
                (occupied, failed), changes_molecule=False)
    if unresolved_reason is None:
        graph.discard(failed, "does not retain a resolvable single-cell barcode")
    else:
        graph.mark_terminal(Terminal(failed.id, Disposition.UNRESOLVED, unresolved_reason))
    reviewed.extend(("gems", "partition"))
    return occupied
