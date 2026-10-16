# Stereo-seq — research notes

## 1. What it is

Stereo-seq combines a patterned array of clonally amplified DNA nanoballs with in-situ poly(A) RNA capture. Each DNB's coordinate identifier is decoded before tissue placement; captured cDNA retains that identifier and a molecular identifier.

## 2. Sources read

- 🟢 Chen et al., *Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays*, Cell 185, 1777–1792.e21 (2022), DOI [10.1016/j.cell.2022.04.003](https://doi.org/10.1016/j.cell.2022.04.003), [published PDF at Cambridge](https://www.repository.cam.ac.uk/bitstreams/1d4e7717-3e50-4448-8bda-a35871dfbfef/download).
- 🟢 The published STAR Methods reagent table (local text lines 1365–1472) prints DNB precursor, fixed oligo, splint, CID sequencing primer, DNB F/R, capture oligo, modified TSO, cDNA PCR primer, Stereo-seq library F/R and splint, read1/read2 and MDA primer. Methods lines 1676–1725 give the chip, reverse transcription, Tn5 and library amplification order; lines 1748–1755 assign CID/MID/read roles.
- 🟢 The 2021 author preprint v3, Methods lines 900–924, independently prints the chip and capture oligos. Its library section names later primers without printing their bases.
- 🟢 The paper-linked [MOSTA detailed protocol](https://db.cngb.org/stomics/assets/html/stereo.seq.html) includes a chip/library procedure and oligo table. Its full post-ligation capture-probe construct is distinct from the short phosphorylated capture oligo printed by Cell.

## 3. Construct-changing steps

1. Generate coordinate-identifier circles, rolling-circle amplify them and pattern the resulting DNBs on the chip.
2. Decode CID-to-position correspondence by sequencing.
3. Ligate molecular-identifier/oligo-dT capture sequence to the decoded DNB oligo.
4. Capture tissue RNA and perform on-chip reverse transcription with template switching.
5. Release and amplify CID–MID-bearing cDNA, construct a DNBSEQ library, circularize and form sequencing DNBs.

## 4. Evidence boundary

🟢 The published Cell table prints the TSO, cDNA-PCR primer, library F/R and splint, and read primers. The in-house Tn5 transposome loading adapter and exact terminal library junction are not specified in the inspected article. 🔴 The named MGI platform guide and original supplemental information still need source review; absence is not accepted. The final panel keeps the unresolved junction inferred.

🟢 There is a published-source discrepancy in the CID precursor: Cell STAR Methods reagent table (local lines 1365–1369) prints `...TTGAACT...`; bioRxiv v3 Methods line 904 and the linked MOSTA detailed protocol oligo table print `...TTGNAACT...`, adding one degenerate N before the fixed handle. The model uses `...TTGNAACT...` because both author documents show that synthesis design; the precursor actually used in the published experiment remains uncertain. The Cell Methods (local line 1684) call the capture tail “22 nt poly-T,” while its reagent table (lines 1412–1414) prints 21 T bases followed by a V anchor. The model follows the printed oligo.
