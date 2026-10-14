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

🟢 GAT5-RT contains an eight-base random tract and anchored poly(T); GAT5-7N has a seven-base random 3′ tract. The protocol reports paired-end 150-base NovaSeq reads. 🟡 The schematic uses inferred TruePrep/Illumina adapter shells because their exact oligo strands, terminal geometry, and chosen i5/i7 index member bases are absent from the author protocol. 🔴 The Science supplement remains unavailable in the source manifest, so any additional variants described there have not been audited.
