#!/usr/bin/env python3
"""FASTQ -> a (guide, UMI) read-count table, plus a per-guide summary.

Reads the Schmierer author-submitted FASTQ, where the record looks like

    @K00110:115:HGLV2BBXX:2:1101:1824:1244 1:N:0:NATCGC+ATCACG
    NCACATTCGACTCATACAGG
                           ^^^^^^ i7 = the 6-nt RSL      ^^^^^^ i5 = sample index

so the guide comes from the sequence line and the UMI from the header. Writes TSV to
$CHEM_DATA; nothing lands in the repo.

Usage
    python3 extract.py $CHEM_DATA/schmierer/input.400000.fq.gz --sample-index ATCACG
    python3 extract.py ... --collapse 1      # adjacency-collapse UMIs (edit distance 1)

Outputs (next to the input, with suffixes):
    *.guide_umi.tsv   guide, umi, reads
    *.guide.tsv       guide, gc, reads, distinct_umi, m1, m2   <- what the model consumes
"""
from __future__ import annotations

import argparse
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from umimodel import gc_fraction  # noqa: E402

GUIDE_LEN = 20


def records(path: Path):
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield h.rstrip("\n"), seq


def parse_index(header: str) -> tuple[str, str]:
    """'... 1:N:0:NATCGC+ATCACG' -> ('NATCGC', 'ATCACG'). ('','') if absent."""
    tail = header.rsplit(":", 1)[-1] if ":" in header else ""
    if "+" not in tail:
        return "", ""
    i7, _, i5 = tail.partition("+")
    return i7.strip(), i5.strip()


def collapse_adjacent(counts: dict[str, int], max_dist: int = 1,
                      ratio: float = 3.0) -> dict[str, int]:
    """Fold rare UMIs into a more abundant neighbour -- the Schmierer authors' rule.

    Their script uses edit distance 1 and a threefold count difference. Implemented here as
    substitutions only (same length), which is what an Illumina index read produces.
    """
    if max_dist <= 0:
        return dict(counts)
    order = sorted(counts, key=lambda u: (-counts[u], u))
    out = dict(counts)
    for umi in reversed(order):
        if umi not in out:
            continue
        n = out[umi]
        best, best_n = None, 0
        for i in range(len(umi)):
            for b in "ACGT":
                if b == umi[i]:
                    continue
                cand = umi[:i] + b + umi[i + 1:]
                m = out.get(cand, 0)
                if m >= ratio * n and m > best_n:
                    best, best_n = cand, m
        if best is not None:
            out[best] += n
            del out[umi]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", type=Path)
    ap.add_argument("--sample-index", default=None,
                    help="keep only reads whose i5 equals this (e.g. ATCACG)")
    ap.add_argument("--umi-len", type=int, default=6)
    ap.add_argument("--collapse", type=int, default=0,
                    help="adjacency-collapse UMIs at this edit distance (0 = off)")
    ap.add_argument("--min-guide-reads", type=int, default=50,
                    help="drop guides below this many reads (error guides)")
    ap.add_argument("--whitelist", type=Path, default=None,
                    help="a *.guide.tsv whose first column lists the real library guides")
    ap.add_argument("--correct-guide", action="store_true",
                    help="reassign reads whose guide is 1 mismatch from exactly one "
                         "whitelist guide. WITHOUT this, GC-rich guides lose more reads to "
                         "sequencing error than GC-poor ones, which mimics amplification "
                         "bias; see the README.")
    ap.add_argument("--max-reads", type=int, default=None,
                    help="stop after this many FASTQ records. Use it to match depth "
                         "between samples: the naive reads-per-UMI statistic is "
                         "depth-dependent, so an unmatched comparison confounds "
                         "sequencing depth with whatever else differs.")
    ap.add_argument("--out-prefix", type=Path, default=None)
    a = ap.parse_args()

    pref = a.out_prefix or a.fastq.with_suffix("").with_suffix("")

    # Build the 1-mismatch correction map, if asked. A neighbour that is 1 mismatch from
    # two different library guides is ambiguous and is dropped rather than guessed.
    fix: dict[str, str] = {}
    white: set[str] = set()
    if a.correct_guide:
        if not a.whitelist:
            raise SystemExit("--correct-guide needs --whitelist")
        # Accept either a published whitelist (guide_id/guide/gene, from
        # download/get_libraries.py) or a *.guide.tsv from an earlier pass. Pick whichever
        # column actually holds 20-mers rather than trusting the header.
        import re as _re
        lines = [ln.rstrip("\n").split("\t") for ln in a.whitelist.open() if ln.strip()]
        body = lines[1:] if lines and not _re.fullmatch(
            r"[ACGT]{%d}" % GUIDE_LEN, lines[0][0].strip().upper()) else lines
        ncol = max(len(r) for r in body)
        best, bestn = 0, -1
        for c in range(ncol):
            n = sum(1 for r in body if len(r) > c
                    and _re.fullmatch(r"[ACGT]{%d}" % GUIDE_LEN, r[c].strip().upper()))
            if n > bestn:
                best, bestn = c, n
        white = {r[best].strip().upper() for r in body
                 if len(r) > best
                 and _re.fullmatch(r"[ACGT]{%d}" % GUIDE_LEN, r[best].strip().upper())}
        print(f"whitelist column {best} of {a.whitelist.name}")
        seen_twice: set[str] = set()
        for g in white:
            for i in range(GUIDE_LEN):
                for b in "ACGT":
                    if b == g[i]:
                        continue
                    nb = g[:i] + b + g[i + 1:]
                    if nb in white:
                        continue
                    if nb in fix and fix[nb] != g:
                        seen_twice.add(nb)
                    else:
                        fix[nb] = g
        for nb in seen_twice:
            fix.pop(nb, None)
        print(f"whitelist {len(white):,} guides; "
              f"1-mismatch map {len(fix):,} unambiguous neighbours "
              f"({len(seen_twice):,} ambiguous, dropped)")

    pair: dict[tuple[str, str], int] = defaultdict(int)
    corrected = 0
    kept = dropped_n = dropped_idx = total = 0
    for header, seq in records(a.fastq):
        if a.max_reads is not None and total >= a.max_reads:
            break
        total += 1
        i7, i5 = parse_index(header)
        if a.sample_index and i5 != a.sample_index:
            dropped_idx += 1
            continue
        umi = i7[:a.umi_len]
        guide = seq[:GUIDE_LEN]
        if "N" in umi or "N" in guide or len(guide) < GUIDE_LEN or len(umi) < a.umi_len:
            dropped_n += 1
            continue
        if white and guide not in white:
            g2 = fix.get(guide)
            if g2 is None:
                dropped_n += 1
                continue
            guide = g2
            corrected += 1
        pair[(guide, umi)] += 1
        kept += 1

    print(f"reads                 {total:,}")
    print(f"  dropped, wrong i5   {dropped_idx:,}")
    print(f"  dropped, N in guide/UMI {dropped_n:,}")
    print(f"  kept                {kept:,}")
    if white:
        print(f"  of which error-corrected {corrected:,} "
              f"({100 * corrected / max(1, kept):.2f}%)")

    per_guide: dict[str, dict[str, int]] = defaultdict(dict)
    for (g, u), n in pair.items():
        per_guide[g][u] = n
    print(f"distinct guides       {len(per_guide):,}")

    guides = {g: us for g, us in per_guide.items() if sum(us.values()) >= a.min_guide_reads}
    print(f"guides >= {a.min_guide_reads} reads   {len(guides):,}")

    if a.collapse:
        before = sum(len(u) for u in guides.values())
        guides = {g: collapse_adjacent(us, a.collapse) for g, us in guides.items()}
        after = sum(len(u) for u in guides.values())
        print(f"UMI collapse d<={a.collapse}      {before:,} -> {after:,} "
              f"({100 * (before - after) / before:.1f}% folded)")

    pu = Path(f"{pref}.guide_umi.tsv")
    with pu.open("w") as fh:
        fh.write("guide\tumi\treads\n")
        for g, us in sorted(guides.items()):
            for u, n in sorted(us.items(), key=lambda kv: -kv[1]):
                fh.write(f"{g}\t{u}\t{n}\n")

    pg = Path(f"{pref}.guide.tsv")
    with pg.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\tm1\tm2\n")
        for g, us in sorted(guides.items()):
            v = list(us.values())
            d = len(v)
            r = sum(v)
            m1 = r / d
            m2 = sum(x * x for x in v) / d
            fh.write(f"{g}\t{gc_fraction(g):.4f}\t{r}\t{d}\t{m1:.6f}\t{m2:.6f}\n")

    print(f"\nwrote {pu}\nwrote {pg}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
