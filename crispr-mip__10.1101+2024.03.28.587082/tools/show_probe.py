#!/usr/bin/env python3
"""Print the CRISPR-MIP probe, its capture and the resulting circle, to scale."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import crisprmip as cm
from crispr import clone_guide
from padlock import capture
from plasmid import read_genbank

VEC = HERE.parents[1] / "lenticrispr-gecko-screen__10.1126+science.1247005" / "ref" / "plasmids" / "addgene-52961_lentiCRISPRv2.gb"
SPACER = cm.EXAMPLE_SPACER              # a Brunello-style, non-G-initiated spacer

# Everything below is measured on the real vector map, so there is nothing to print without
# it. That map is third-party material and is not committed.
if not VEC.exists():
    print(f"cannot run: missing {VEC.relative_to(HERE.parents[1])}"
          " -- the lentiCRISPRv2 map (Addgene #52961). It is third-party material and is not"
          " committed; see lenticrispr-gecko-screen__10.1126+science.1247005/ref/plasmids/"
          "MANIFEST.md for how to re-obtain it.", file=sys.stderr)
    sys.exit(1)

probe = cm.probe()
v = clone_guide(read_genbank(str(VEC)), SPACER)
c = capture(v, probe)[0]
E, L, B = cm.EXT_ARM, cm.LIG_ARM, probe.backbone


def rule(t):
    print(f"\n{t}\n" + "=" * len(t))


rule(f"1. THE PROBE AS ORDERED  ({len(probe)} nt of ssDNA, one molecule)")
print(f"5'-p {L} {B} {E} -3'OH")
print(f"     {'^' * len(L)} {'^' * len(B)} {'^' * len(E)}")
print(f"     {'ligation arm, ' + str(len(L)) + ' nt':<{len(L)}} "
      f"{'backbone, ' + str(len(B)) + ' nt -- the Ns':<{len(B)}} "
      f"{'extension arm, ' + str(len(E)) + ' nt':<{len(E)}}")
print(f"     {'Tm %.1f C' % probe.arm_tms()[1]:<{len(L)}} "
      f"{'NOT complementary to the target':<{len(B)}} {'Tm %.1f C' % probe.arm_tms()[0]}")

rule(f"2. THE {len(B)}-nt BACKBONE  (Table S2)")
for n, l in (("UMI", cm.UMI_LEN), ("Read 2 site (TruSeq)", len(cm.READ2_SITE)),
             ("i7 index", cm.INDEX_LEN), ("Read 1 site (TruSeq, 3' part)", len(cm.READ1_SITE))):
    print(f"   {n:<30s} {l:>3d} nt")
print("\n   NOTE: P5/P7 are NOT in the backbone. The paper puts the PCR primer sites")
print("   INSIDE the captured sequence, so an unextended probe cannot amplify.")

rule("3. THE CAPTURE  (arms anneal; 112 nt of target lies between them)")
fill = c.fill
print(f"target top strand, 5'->3':")
print(f"   ...{E} {fill[:1]} {fill[1:21]} {fill[21:41]}...{fill[-8:]} {L}...")
print(f"      {'^' * len(E)} {'^':<1} {'^' * 20} {'^' * 20}   {'^' * 8} {'^' * len(L)}")
print(f"      {'ext arm binds':<{len(E)}} {'+1':<1} {'sgRNA spacer':<20} {'scaffold ->':<20}   {'...term':<8} {'lig arm binds'}")
print(f"\n      gap the polymerase must fill: {c.gap} nt"
      f"   (111 if the spacer already starts with G)")
print(f"      the probe itself anneals to the OPPOSITE strand, arching over this gap")

rule(f"4. THE CIRCLE  ({len(c.circle)} nt, covalently closed -- the only thing that survives)")
print(f"      .--> {E} [{c.gap} nt newly synthesised] {L} {B} --.")
print(f"      |    {'^' * len(E)}  {'^' * 28}  {'^' * len(L)} {'^' * len(B)}  |")
print(f"      |    {'ext arm':<{len(E)}}  {'copy of the target: +1 G, spacer,':<28}  "
      f"{'lig arm':<{len(L)}} {'the Ns':<{len(B)}}  |")
print(f"      |    {'':<{len(E)}}  {'scaffold, terminator':<28}  {'^ ligation':<{len(L)}} {'':<{len(B)}}  |")
print(f"      '{'-' * 100}--'")
print(f"\n      circle = probe ({len(probe)}) + gap ({c.gap}) = {len(c.circle)} nt")
print(f"      the UMI in those Ns is now attached to ONE captured genomic molecule,")
print(f"      before any amplification -- which is what makes the readout a count.")
