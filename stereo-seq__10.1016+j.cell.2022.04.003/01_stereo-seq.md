# Stereo-seq — research notes

## 1. What it is

Stereo-seq combines a patterned array of clonally amplified DNA nanoballs with in-situ poly(A) RNA capture. Each DNB's coordinate identifier is decoded before tissue placement; captured cDNA retains that identifier and a molecular identifier.

## 2. Sources read

- 🟢 Chen et al., *Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays*, Cell 185, 1777–1792.e21 (2022), DOI [10.1016/j.cell.2022.04.003](https://doi.org/10.1016/j.cell.2022.04.003), [open paper](https://www.devo-evo.com/assets/PDFs/stereo_seq.pdf).
- 🟢 Main and STAR Methods establish chip-DNB generation and decoding, UMI/poly(dT) addition, tissue capture, RT, release, PCR and DNBSEQ read roles.

## 3. Construct-changing steps

1. Generate coordinate-identifier circles, rolling-circle amplify them and pattern the resulting DNBs on the chip.
2. Decode CID-to-position correspondence by sequencing.
3. Ligate molecular-identifier/oligo-dT capture sequence to the decoded DNB oligo.
4. Capture tissue RNA and perform on-chip reverse transcription with template switching.
5. Release and amplify CID–MID-bearing cDNA, construct a DNBSEQ library, circularize and form sequencing DNBs.

## 4. Evidence boundary

🟢 The source establishes topology, order and barcode roles. 🔴 Complete commercial chip oligos and current sequencing-primer bases are not published. The page explicitly draws these as inferred role placeholders and states why the final cPAS site cannot be shown base-for-base.
