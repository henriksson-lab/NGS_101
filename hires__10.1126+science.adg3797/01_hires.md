# HiRES — research notes

## Sources read

- 🟢 Liu et al. 2023, DOI `10.1126/science.adg3797`, article abstract and Methods metadata. The Science full text and supplement were not available through the source fetcher.
- 🟢 Chen et al. 2023, DOI `10.21769/BioProtoc.4886`, [author detailed protocol](https://bio-protocol.org/pdf/bio-protocol4886.pdf), Materials p. 4 and procedure pp. 8–11. Its abstract explicitly identifies the Science paper as the assay it implements.

## Construct path and oligos

🟢 The author protocol reverse-transcribes RNA **in situ before** MboI digestion and DNA proximity ligation. It then sorts individual nuclei into wells and jointly preamplifies contact DNA and transcript cDNA using GAT5-family primers. A Vazyme TruePrep V2 transposome fragments the mixed amplicons. RNA-P7 slightly enriches transcript reads before well-specific i5 and i7 dual indexes complete a **single joint library**. Contact and RNA reads are separated during analysis; they are not separate products of three cell-index ligation rounds. The previous page's three-round scheme was unsupported by this detailed source and has been replaced.

🟢 The four PAGE-purified primers printed in the Bio-protocol Materials section, with source line wraps rejoined, are:

- GAT5-RT: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGTATTGCGCAATGNNNNNNNNTTTTTTTTTTTTTTTVN-3′`.
- GAT5-7N: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGTNNNNNNN-3′`.
- GAT5: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGT-3′`.
- RNA-P7: `5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGTTGAGGTAGTATTGCGCAATG-3′`.

🟢 GAT5-RT contains an eight-base random tract and anchored poly(T); GAT5-7N has a seven-base random 3′ tract. The protocol reports paired-end 150-base NovaSeq reads. The authors' [analysis repository](https://github.com/skelviper/HiRES) at commit `498941193bdb879c68ebb5d11f1157270ffcaf66`, `HiRES_preprocess_pipeline/rules/HiRES_split.rules` lines 1–37, classifies paired reads using the RNA-arm motif `XGGTTGAGGTAGTATTGCGCAATG` and then filters RNA pairs for the random-eight/poly(T) motif. This supports the two computational read classes from one joint library.

🟢 The Vazyme TD501–503 manual (extracted PDF lines 145–178 and 295–311) says its tagmentation mix holds two adapter species at equal molar ratio and prints a standard finished TruePrep library with P5, i5, Read 1 arm, insert, Read 2 arm, i7 and P7. The [Illumina UD index reference](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Illumina-UDIndexes.htm) prints index adapter arms and plate-member sequences. 🟡 The schematic still marks the HiRES-specific shell inferred because the loaded transposome oligos, their exact terminal geometry, the RNA-P7 junction and chosen per-well i5/i7 members are not established by these sources. 🔴 The defining Science article and supplement remain unavailable in the source manifest.
