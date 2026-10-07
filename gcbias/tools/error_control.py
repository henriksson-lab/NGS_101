#!/usr/bin/env python3
"""Does guide-read sequencing error vary with GC? A control for the GC-bias analysis.

Reads whose 20-nt guide carries a basecalling error do not match the library and are
discarded. If GC-rich guides accumulate more such errors, they lose more reads, which looks
exactly like reduced amplification. This measures that directly, by counting reads sitting
1 mismatch from each library guide.

The last column converts the error rate into the log2 read-count change it would cause on
its own. Compare that with the observed GC deficit -- but note it is an UPPER bound on the
confound, because errors remove distinct UMIs as well as reads, and the statistic of
interest (reads at fixed molecule count) is affected by the difference, not by the raw read
loss. The honest way to settle it is extract.py --correct-guide, which recovers the reads.

Usage
    python3 error_control.py $CHEM_DATA/schmierer/input.25000000.fq.gz \\
        $CHEM_DATA/schmierer/input.collapsed.guide.tsv
"""
import sys, gzip, math, statistics as st
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
sys.path.insert(0, str(_P(__file__).resolve().parents[1] / 'lib'))
sys.path.insert(0, str(_P(__file__).resolve().parents[2] / 'lib'))
from collections import Counter
from pathlib import Path
from gcbias import ols
from qualitative import read_guides, spearman

fq, gtsv = sys.argv[1], sys.argv[2]
lib = {r["guide"]: r["gc"] for r in read_guides(Path(gtsv))}
cnt = Counter()
with gzip.open(fq, "rt") as fh:
    while True:
        h = fh.readline()
        if not h: break
        s = fh.readline().rstrip("\n"); fh.readline(); fh.readline()
        if len(s) >= 20: cnt[s[:20]] += 1
print(f"distinct 20-mers {len(cnt):,}; library guides {len(lib):,}")

rows=[]
for g, gc in lib.items():
    exact = cnt.get(g, 0)
    if exact < 100: continue
    d1 = 0
    for i in range(20):
        for b in "ACGT":
            if b == g[i]: continue
            n = cnt.get(g[:i] + b + g[i+1:], 0)
            if n and n < exact * 0.2:      # neighbours far rarer than the parent = errors
                d1 += n
    rows.append(dict(gc=gc, exact=exact, d1=d1, frac=d1/(exact+d1)))
print(f"guides analysed {len(rows):,}")
gcv=[r["gc"] for r in rows]; fr=[r["frac"] for r in rows]
b,_,_,se = ols(gcv, fr)
print(f"\nfraction of a guide's reads that are 1-mismatch errors:")
print(f"  slope vs GC {b:+.5f} +/- {se:.5f}   spearman {spearman(gcv,fr):+.4f}")
print(f"\n{'GC band':14s} {'n':>6s} {'median error frac':>18s} {'rel. to mid band':>18s}")
bnds=[(0,.45),(.45,.55),(.55,.60),(.60,.70),(.70,1.01)]
meds=[]
for lo,hi in bnds:
    ch=[r for r in rows if lo<=r["gc"]<hi]
    meds.append(st.median([r["frac"] for r in ch]) if ch else None)
ref=meds[1]
for (lo,hi),m in zip(bnds,meds):
    ch=[r for r in rows if lo<=r["gc"]<hi]
    if m is None: continue
    # effect on log2 read count if these reads were lost
    delta = math.log2((1-m)/(1-ref))
    print(f"{lo:.2f}-{hi:.2f}     {len(ch):>6,} {m:>18.5f} {delta:>+18.5f}")
print("\n(last column = log2 read-count change this error rate alone would cause,")
print(" relative to the mid-GC band. Compare with the observed -0.055 at high GC.)")
