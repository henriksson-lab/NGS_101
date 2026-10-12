"""Shared molecular geometry for Tn5 adapter loading and transfer."""
from __future__ import annotations

import nextera as nx
from chemdraw import Scene, Segment


def loaded_tn5_adapter(top: list[Segment], *, mosaic_end: str,
                       label: str = "Tn5 top oligo",
                       bottom_name: str = "Tn5MErev",
                       bottom_label: str = "Tn5 mosaic-end bottom",
                       bottom_sequence: str = nx.ME_RC,
                       mod5: str = "p", mod3: str = "") -> Scene:
    """Draw a loaded Tn5 adapter with only its terminal ME duplexed.

    The defining structural invariant is enforced here: the named mosaic end must be the
    3′-terminal segment of the top oligo and must carry the canonical 19-nt ME.  Any
    protocol-specific handle, barcode or connector therefore remains a real single-stranded
    overhang rather than accidentally appearing duplexed.
    """
    if not top:
        raise ValueError("a loaded Tn5 adapter needs a top oligo")
    if top[-1].name != mosaic_end:
        raise ValueError("the Tn5 mosaic end must be the 3′-terminal top segment")
    if top[-1].top != nx.ME:
        raise ValueError("the loaded Tn5 top oligo must end in the canonical 19-nt ME")
    if bottom_sequence != nx.ME_RC:
        raise ValueError("the Tn5 mosaic-end bottom must be the canonical reverse complement")
    sc = Scene()
    sc.strand("Tn5 top", top, label=label)
    sc.anneal("Tn5 bottom", [Segment(bottom_name, bottom_sequence, "me")],
              to="Tn5 top", pair=(bottom_name, mosaic_end), label=bottom_label,
              mod5=mod5, mod3=mod3)
    sc.labels("Tn5 top")
    return sc
