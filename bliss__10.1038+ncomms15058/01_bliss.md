# BLISS — research notes

Defining source: Yan *et al.*, [10.1038/ncomms15058](https://doi.org/10.1038/ncomms15058).

## Evidence used

- 🟢 Fixed cells or tissue are immobilized, DSB ends are blunted, and a double-stranded adapter is ligated directly to breaks in situ.
- 🟢 The adapter contains a T7 promoter, RA5 sequence, sample barcode and an 8–12-nt UMI.
- 🟢 T7 in-vitro transcription linearly amplifies break-adjacent sequence before completion of the Illumina RNA-style library.
- 🟢 Supplementary Data 1 prints both strands of eleven classic BLISS adapters
  and ten longer-UMI adapter variants. The page uses classic adapter A1 as its
  representative exact duplex.
- 🟢 Supplementary Fig. 1a marks the sample-barcode end as the blunt ligation
  site and shows a 5′ phosphate on its bottom strand. The opposite bottom-strand end has
  a two-base 3′ overhang that prevents head-to-tail adapter concatemers.
- 🟡 The downstream commercial RNA-library finishing reagents are named by the
  procedure but are not part of Supplementary Data 1; the page does not infer their
  lot-specific ordering sequences from the exact break adapter.

## Protocol boundary

The page covers original BLISS. Its suspension implementation, sBLISS, is a separately named protocol and is not silently merged here.
