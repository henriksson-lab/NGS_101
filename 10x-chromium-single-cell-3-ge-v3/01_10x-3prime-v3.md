# 10x Chromium Single Cell 3' Gene Expression v3 research notes

The source audit and chemistry comparison are preserved in the combined
[v2/v3/v4 note](../10x-chromium-single-cell-3-ge-v2/01_10x-3prime-v2-v4.md#v3-cg000183-appendix-oligonucleotide-sequences).

The primary chemistry source is 10x Genomics user guide **CG000183 Rev A**.


## Source-led model corrections (2026-10-16)

🟢 Both the cDNA-amplification and sample-index-PCR primer architectures are
printed in the version-specific manufacturer source (v2: CG000108 Rev A Fig. 4;
v3: CG000183 Rev A Appendix pp. 54–55). The schematic now uses these sequences
in its PCR geometry instead of claiming that they were unavailable. Copied TSO
bases are derived from the printed TSO, and the post-template-switch state
includes that handle before recovery and amplification. These are repository
model corrections.

For v3, CG000183 Rev A Appendix p. 54 prints a 32-base adaptor ending `CAGTCA`,
while its ligation-product drawing ends `CAGTCAC`. The public page now shows
this source discrepancy. The model retains the printed 32-base input; the
sample-index reverse primer supplies the terminal C in the completed library.
The guide does not resolve which pre-PCR drawing reflects the supplied oligo.


## Scope of the adaptor discrepancy

The CG000183 terminal-base notice applies to the source-specific input adaptor
and its pre-index ligation state. An inheritance audit also inspected the
Cell Surface Protein, direct CRISPR-capture, ZipSeq and OAK graphs. The first
two displayed feature branches install their Read 2 ends by tag copying/PCR;
the latter two currently show final-library endpoints and do not assert the
32-base input-adaptor assembly. A final-library `CAGTCAC` arm alone is
compatible with both input lengths because index PCR completes it. Those
pages therefore do not receive this input-adaptor discrepancy notice merely
because their final libraries share canonical Illumina bases.


## Bounded source pass (2026-10-17, round breadth-first-06)

Re-read CG000183 Rev A Appendix oligonucleotide sequences (extracted lines
4868–5016, pp. 54–55), Steps 2.1–3.6 (TOC lines 98–132; 3.1 program lines
3750–3800; 3.2 lines 3890–3905; 3.4 lines 4028–4035; 3.6 lines 4225–4235) and
the sequencing table (lines 4370–4392), with the bundle scraper. All printed
GE oligos still match the model; read layout 28/8/0/91 confirmed. The two
existing published-discrepancy notices (poly(dT)VN label vs T30 appendix oligo;
32-base adaptor vs `CAGTCAC` ligation product) were rechecked and remain.

Model corrections:

- 🟢 The graph now carries the guide's bead cleanups: Dynabeads MyOne SILANE
  post GEM-RT cleanup (2.1), 0.6X SPRIselect cDNA cleanup (2.3), double-sided
  0.6X/0.8X after fragmentation/end repair/A-tailing (3.2), 0.8X post-ligation
  (3.4) and double-sided 0.6X/0.8X after sample-index PCR (3.6), with removed
  fractions declared as discarded waste; the sequenced library is the 3.6
  product.
- 🟢 The page now shows the second gel-bead primer printed on p. 54,
  `GTCAGATGTGTATAAGAGACAG` (Nextera Read 1) · 16 N barcode · 12 N UMI ·
  `TTGCTAGGACCGGCCTTAAAGC` (Capture Sequence 1), which the guide labels "Primer
  not used for 3' Gene Expression". It is not drawn as a reaction participant;
  whether it yields side products in GE GEMs is not stated.

Still open: RT primer and TSO graph inputs have no exact molecular state (the
shared Chromium operation does not accept one); later CG000183 revisions not
compared.

## Source-owned graph migration, 2026-10-18

Source-owned graph now declares partitioning, RT, template switching, recovery, enzymatic fragmentation/end preparation and sequencing locally. Exact printed PCR and dT-adapter architectures, all five cleanup branches, 28/8/91 acquired-cycle plan and both public source-discrepancy notices are preserved. First-strand cDNA is drawn single stranded. Empty/multiplet GEM outcomes are unresolved pooled products, not a physical discard justified only by barcode resolvability.

Bounded source scope: CG000183 Rev A source original Appendix pp54–55 (extracted lines4921–5078), Step3.1 vicinity lines3740–3800; existing cleanup declarations retained without expanding this bounded source reading. No accepted missing-information claim is added.

Remaining work: Exact transcript-dependent fragmentation/end-preparation product states and molecular reaction relations remain unresolved; Later CG000183 revisions not compared.
