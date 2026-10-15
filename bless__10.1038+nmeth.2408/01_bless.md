# BLESS — research notes

Crosetto et al., *Nature Methods* (2013), [doi:10.1038/nmeth.2408](https://doi.org/10.1038/nmeth.2408), PMCID: PMC3651036.

## Sources read

- 🟢 Author-manuscript full text from PMC.
- 🟢 Supplementary Methods and Supplementary Table 5 from the PMC Cloud copy (`NIHMS449608-supplement-1.pdf`).

## Molecular interpretation

1. 🟢 Cells are fixed and gently lysed while nuclei remain intact. Native DSB ends are blunted and 5′ phosphorylated.
2. 🟢 The 91-nt proximal hairpin is phosphorylated, carries an internal biotin-dT and is ligated to DSBs in situ. Its terminal `TCGAGGTAGTA` sequence identifies a proximal end.
3. 🟢 DNA is extracted, fragmented and captured on streptavidin. This fragmentation creates the free distal end; it is not a second native break.
4. 🟢 Captured fragments are blunted and ligated to the distinct distal hairpin. The later L3 design leaves `TCGAGACGACG` as the distal tag.
5. 🟢 I-SceI digestion releases the hairpin-bound construct, proximal- and distal-specific primers amplify it, and XhoI digestion leaves the two short end tags.
6. 🟢 The paper's Methods, "Next-generation sequencing," specifies TruSeq DNA Sample Preparation Kit v2 (Illumina), without fragmentation or library size selection. Reads containing the proximal tag report the native DSB; reads containing the distal tag start at the later fragmentation end.
7. 🟡 The public diagram draws the final sequencing-kit shell in the familiar TruSeq orientation solely to locate the standard sequencing primers. The paper does not print the complete kit adapter oligos, so the outer shell is visibly inferred; the BLESS insert and tags are not.

## Exact oligos

Supplementary Table 5 gives L1, L3, P1 and P3. Supplementary Fig. 1 and the analysis section independently identify the retained 11-base proximal and distal tags. The internal proximal modification is written `[BtndT]` in the source and is rendered `/iBiodT/` on the page.

The first published experiments also used distal linker L2. The later L3 carries a separate distal barcode and is the form represented on the schematic because it is the specificity-improved BLESS construct. Supplementary Table 5 also prints P4/P5 for a Roche 454 branch; those are not omitted Illumina PCR primers. For the illustrated Illumina branch, L1/L3 and P1/P3 are complete through the BLESS insert. The downstream kit is explicitly TruSeq DNA Sample Preparation Kit v2; the paper does not print its adapter oligos or specify the selected index member.
