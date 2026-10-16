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
