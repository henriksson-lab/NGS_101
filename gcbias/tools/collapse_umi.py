#!/usr/bin/env python3
"""UMI error-collapse, directional method -- the authors' UMIClusterer(threshold=2).

WHY THIS IS NOT OPTIONAL HERE. Every read whose UMI is miscalled lands on a fresh UMI
(13 nt = 67M values, so no collisions). Error UMIs are therefore proportional to READS,
and reads are GC-dependent -- so the inflation of the molecule count is GC-dependent too,
and it pushes molecules and reads apart in exactly the direction of the claimed
"CRISPR-MIP corrects GC bias" effect. At ~11.6 reads/UMI and 0.1-0.2%/base, 15-30% of
observed "molecules" are error artefacts against a claimed effect of 0.027 log2. Without
collapse the correction and its largest artefact are confounded.

METHOD (umi_tools "directional", reimplemented since umi_tools is not installed):
within each guide, sort UMIs by count descending; a node A absorbs B when
    hamming(A,B) <= threshold   AND   count_A >= 2*count_B - 1
Absorption is transitive through the graph, and each UMI is assigned once.
Pairwise comparison within a guide is used rather than neighbour enumeration: with tens
to hundreds of UMIs per guide that is far cheaper than generating 741 distance-2
neighbours per UMI.

Usage
    python3 collapse_umi.py <guide_umi.tsv> --threshold 2 --out-prefix <pref>
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from umimodel import gc_fraction  # noqa: E402


def hamming_le(a: str, b: str, k: int) -> bool:
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def directional(counts: dict[str, int], threshold: int) -> list[int]:
    """Return collapsed counts: one entry per surviving UMI group."""
    umis = sorted(counts, key=lambda u: (-counts[u], u))
    parent = {u: u for u in umis}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for i, a in enumerate(umis):
        ca = counts[a]
        for b in umis[i + 1:]:
            cb = counts[b]
            # counts are sorted DESCENDING, so cb falls as b advances and
            # (2*cb - 1) falls with it. An early failure is therefore followed by
            # successes -- breaking here would silently skip most true children.
            if ca < 2 * cb - 1:
                continue
            if find(a) == find(b):
                continue
            if hamming_le(a, b, threshold):
                parent[find(b)] = find(a)
    tot: dict[str, int] = defaultdict(int)
    for u in umis:
        tot[find(u)] += counts[u]
    return list(tot.values())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("guide_umi", type=Path)
    ap.add_argument("--threshold", type=int, default=2)
    ap.add_argument("--out-prefix", type=Path, default=None)
    a = ap.parse_args()

    per: dict[str, dict[str, int]] = defaultdict(dict)
    with a.guide_umi.open() as fh:
        fh.readline()
        for line in fh:
            g, u, n = line.rstrip("\n").rsplit("\t", 2)
            per[g][u] = int(n)
    before = sum(len(v) for v in per.values())
    reads = sum(sum(v.values()) for v in per.values())

    out: dict[str, list[int]] = {}
    for g, cd in per.items():
        out[g] = directional(cd, a.threshold) if len(cd) > 1 else list(cd.values())
    after = sum(len(v) for v in out.values())
    print(f"{a.guide_umi.name}")
    print(f"  guides {len(per):,}   reads {reads:,}")
    print(f"  UMIs {before:,} -> {after:,}  ({100 * (before - after) / before:.1f}% collapsed"
          f" at edit distance <= {a.threshold})")
    print(f"  reads per UMI {reads / before:.3f} -> {reads / after:.3f}")

    pref = a.out_prefix or a.guide_umi.with_name(
        a.guide_umi.name.replace(".guide_umi.tsv", "") + ".collapsed")
    pg = Path(f"{pref}.guide.tsv")
    with pg.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\tm1\tm2\n")
        for g, v in sorted(out.items()):
            d, r = len(v), sum(v)
            fh.write(f"{g}\t{gc_fraction(g):.4f}\t{r}\t{d}\t{r / d:.6f}\t"
                     f"{sum(x * x for x in v) / d:.6f}\n")
    print(f"  wrote {pg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
