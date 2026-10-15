#!/usr/bin/env python3
"""Audit narrow molecular claims in the ten ownership pilots and boundary examples.

This is an independent analysis of ordinary Python protocol declarations, not a
source-transcription check or a gate for publishing a protocol. Use --details to see
every claim and its precise scope. A contradiction exits nonzero.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "lib"), str(ROOT / "tools")]

from chemistry_validation import Status, check_graph, warnings  # noqa: E402
from selftest_source_oligos import (clear_protocol_imports, graph_module_paths,
                                    load)  # noqa: E402


PILOTS = (
    "smart-seq2__10.1038+nmeth.2639",
    "paiso-seq__10.1038+s41467-019-13228-9",
    "sccat-seq__10.1038+s41467-018-08205-7",
    "scm-t-seq__10.1038+nmeth.3728",
    "scmt-seq__10.1186+s13059-016-0950-z",
    "scnmt-seq__10.1038+s41467-018-03149-4",
    "cel-seq2__10.1186+s13059-016-0938-8",
    "sci-atac-seq-2018__10.1038+nature25981",
    "splash__10.1016+j.molcel.2016.04.028",
    "sf4su-seq__10.1093+nar+gkag652",
)
BOUNDARIES = (
    "smmip__10.1101+gr.147686.112",  # hybridization and gap-sealing circle
    "fisseq__10.1126+science.1250212",  # circle and rolling-circle copy
    "small-seq__10.1038+nbt.3701",  # RNA ligation and dephosphorylation
    "malbac__10.1126+science.1229164",  # intramolecular fold with unknown bases
    "in-situ-hi-c__10.1016+j.cell.2014.11.021",  # restriction/fill/contact join
)


def audit(name: str):
    paths = graph_module_paths(ROOT / name)
    if not paths and name == "in-situ-hi-c__10.1016+j.cell.2014.11.021":
        paths = [ROOT / name / "tools" / "insitu_hic.py"]
    if len(paths) != 1:
        raise ValueError(f"{name}: expected one graph declaration, found {len(paths)}")
    path = paths[0]
    module = load(path)
    sys.path.insert(0, str(path.parent))
    try:
        graph = (module.protocol_graph() if hasattr(module, "protocol_graph")
                 else module.typed_graph())
        return graph, check_graph(graph)
    finally:
        sys.path.pop(0)
        clear_protocol_imports()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dirs", nargs="*", help="protocol directories; default: ten pilots and five boundary examples")
    parser.add_argument("--details", action="store_true", help="show each scoped claim")
    args = parser.parse_args()
    names = args.dirs or (*PILOTS, *BOUNDARIES)
    total = Counter()
    for name in names:
        graph, claims = audit(name)
        counts = Counter(claim.status.value for claim in claims)
        total.update(counts)
        summary = ", ".join(f"{status.value}={counts[status.value]}" for status in Status)
        print(f"{name}: {len(graph.steps)} steps; {summary}")
        for warning in warnings(graph):
            print(f"  WARNING {warning.target}: {warning.message}")
        for claim in claims:
            if args.details or claim.status is Status.CONTRADICTORY:
                print(f"  {claim.status.value.upper()} {claim.target}: "
                      f"{claim.question} — {claim.reason}")
    print("Total claims: " + ", ".join(f"{status.value}={total[status.value]}"
                                       for status in Status))
    return bool(total[Status.CONTRADICTORY.value])


if __name__ == "__main__":
    raise SystemExit(main())
