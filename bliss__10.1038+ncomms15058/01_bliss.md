# BLISS — research notes

Defining source: Yan *et al.*, [10.1038/ncomms15058](https://doi.org/10.1038/ncomms15058).

## Evidence used

- 🟢 Fixed cells or tissue are immobilized, DSB ends are blunted, and a double-stranded adapter is ligated directly to breaks in situ.
- 🟢 The adapter contains a T7 promoter, RA5 sequence, sample barcode and an 8–12-nt UMI.
- 🟢 T7 in-vitro transcription linearly amplifies break-adjacent sequence before completion of the Illumina RNA-style library.
- 🟢 The authors' [Protocol Exchange v2](https://assets-eu.researchsquare.com/files/nprot-5597/v2/17ba35a1-4de4-4f0a-8953-6004ecad1900.pdf?c=1631828212), pp. 3 and 10, specifies custom-synthesized TruSeq Small RNA RA3 ligation to IVT RNA, RTP reverse transcription, then RP1 and a chosen RPI index primer for PCR. The [Illumina sequence table](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm) prints RA3 `TGGAATTCTCGGGTGCCAAGG`, RTP `GCCTTGGCACCCGAGAATTCCA`, and RP1 `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`, plus the RPI family. The original method does not select one RPI member.
- 🟢 Supplementary Data 1 prints both strands of eleven classic BLISS adapters
  and ten longer-UMI adapter variants. The page uses classic adapter A1 as its
  representative exact duplex.
- 🟢 Supplementary Fig. 1a marks the sample-barcode end as the blunt ligation
  site and shows a 5′ phosphate on its bottom strand. The opposite bottom-strand end has
  a two-base 3′ overhang that prevents head-to-tail adapter concatemers.
- 🟡 RA3/RTP/RP1/RPI are a distinct library-finishing family; Supplementary Data 1
  lists the break adapters. The model shows RA3 ligation before RT. The method says
  these oligos were custom synthesized from TruSeq Small RNA sequences, but does not
  specify the custom RA3 end modifications or the RPI index chosen for a sample.

## Protocol boundary

The page covers original BLISS. Its suspension implementation, sBLISS, is a separately named protocol and is not silently merged here.
