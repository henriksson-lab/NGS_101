#!/usr/bin/env python3
"""Reproduce every number in ANALYSIS.md from the two maps in this folder."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[0] / "crispr-mip__10.1101+2024.03.28.587082" / "tools"))

import crisprmip as cm
from crispr import SCAFFOLD_V1, clone_guide
from padlock import capture
from plasmid import BsmBI, Plasmid, amplify, digest, find_both, read_genbank

P1F = "AATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG"
P1R = "TCTACTATTCTTTCCCCTGCACTGTTGTGGGCGATGTGCGCTCTG"
P5 = ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
      "TTGTGGAAAGGACGAAACACCG")
P7 = ("CAAGCAGAAGACGGCATACGAGATATACTCAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
      "TCTACTATTCTTTCCCCTGCACTGT")
SPACER = "ATCGATCGATCGATCGATCG"

# The client supplied an "Addgene Brunello plasmid". Brunello ships in BOTH of these
# backbones -- #73179 in lentiCRISPRv2, #73178 in lentiGuide-Puro -- and the MIP probe
# only works in one of them. See ira1.md section 8.
# The two maps are SnapGene exports: third-party, so gitignored and not in a fresh
# clone. Say exactly what is missing and where it comes from rather than tracebacking.
MAPS = (("lentiCRISPR v2.dna", "lentiCRISPRv2  = Brunello pool #73179",
         "Addgene #73179 (Brunello in lentiCRISPRv2), full sequence, SnapGene export"),
        ("lentiGuide-Puro.dna", "lentiGuide-Puro = Brunello pool #73178",
         "Addgene #73178 (Brunello in lentiGuide-Puro), full sequence, SnapGene export"))
_absent = [(fn, src) for fn, _, src in MAPS if not (HERE / fn).exists()]
if _absent:
    print("cannot reproduce the analysis: the plasmid maps are not in the repository\n"
          "(they are third-party, so they are gitignored -- see to_debug/ira1.md)\n",
          file=sys.stderr)
    for fn, src in _absent:
        print(f"  missing to_debug/{fn}\n      obtain from: {src}", file=sys.stderr)
    sys.exit(1)

for fn, name, _src in MAPS:
    p = read_genbank(str(HERE / fn))
    cuts = digest(p, BsmBI)
    print(f"\n=== {name}: {len(p)} bp ===")
    print(f"  BsmBI sites {len(cuts)}"
          + (f", filler {cuts[1].top_cut - cuts[0].top_cut} bp -> EMPTY vector"
             if len(cuts) == 2 else ""))
    for label, tpl in (("empty ", p), ("cloned", clone_guide(p, SPACER))):
        a1 = amplify(tpl, P1F, P1R, min_anneal=18)
        if not a1:
            print(f"  {label}: PCR1 none")
            continue
        a2 = amplify(Plasmid("x", a1[0].seq, False, []), P5, P7, min_anneal=18)
        # only meaningful for the as-deposited map: clone_guide rotates the origin itself
        wrap = " (crosses the map origin)" if (a1[0].wraps_origin and tpl is p) else ""
        print(f"  {label}: PCR1 {a1[0].length:>5d} bp{wrap:<25s}"
              f"  PCR2 {a2[0].length if a2 else 0:>5d} bp")
        caps = capture(tpl, cm.probe(), max_gap=50000)
        print(f"          MIP gap {caps[0].gap if caps else 'NO CAPTURE':>5}"
              f"   ligation arm present: {bool(find_both(p.seq, cm.LIG_ARM, True))}")
    i = find_both(p.seq, SCAFFOLD_V1[:30], True)[0][0] + len(SCAFFOLD_V1)
    print(f"  after the scaffold: {p.sub(i, i + 40)}")
