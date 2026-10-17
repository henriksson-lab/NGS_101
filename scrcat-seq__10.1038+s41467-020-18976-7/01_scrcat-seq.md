# scRCAT-seq

Defining paper: Liu et al., *Nature Communications* 11, 5148 (2020), DOI
[10.1038/s41467-020-18976-7](https://doi.org/10.1038/s41467-020-18976-7). Read from the PMC
open-access copy PMC7555861 (JATS XML and PDF).

## Sources read (round 2026-10-17-breadth-first-17)

- JATS `.xml.txt`: Results line 65 (method overview), Fig. 1a legend line 71, Methods line 149
  ("scRCAT-seq library construction for a single cell"), lines 169–171 (TSS and TES read
  processing). Rendered Fig. 1a on PDF page 3.
- Oligo scraper run over the bundle: library hits are the Methods line 149/169 sequences
  (and their PDF copies). Supplementary Table 7 (MOESM1 lines 3085–3131) holds cloning primers
  for novel-gene validation, and MOESM5 "Fig. 4b" lists 16-nt 10x-style cell barcodes; neither
  is library chemistry.
- Scope: the 10x Chromium 3′ v3 / 5′ application (Results line 113, Methods line 165) uses
  separate vendor-kit libraries analysed with the scRCAT-seq pipeline; it is not drawn here.
- Cited Smart-seq2 protocol (Picelli et al., *Nat. Protoc.* 2014, ref 21): not openly
  available (manifest of the SMART-seq2 bundle records a manual publisher URL); not read.

## Oligos (verbatim)

🟢 RT primer `5′-AAGCAGTGGTATCAACGCAGAGTN8[16 bps of cell barcode]T30VN-3′` (Methods line 149).
The N8 is the UMI; the 16-bp field is the cell barcode.

🟢 P5 index primer
`5′AATGATACGGCGACCACCGAGATCTACAC[8 bps of index]TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGGTATCAACGCAGAGT`
and P7 index primer `5′CAAGCAGAAGACGGCATACGAGAT[8 bps of index]GTCTCGTGGGCTCGG` (Methods line 149).

🟢 TSS reads are selected by the tag `5′-GTGGTATCAACGCAGAGTACATGGG-3′`, and the TSO sequence
`5′-GTGGTATCAACGCAGAGTACAT-3′` is trimmed (Methods line 169). The full TSO is not printed.

## How they interlock (🟡 derived; computed in the module)

The P5 index primer ends in `GTGGTATCAACGCAGAGT`, the last 18 bases of the RT-primer handle
and the first 18 bases of the TSO read part. One P5 primer therefore primes the tag at either
cDNA end, and the P7 index primer's `GTCTCGTGGGCTCGG` primes the Nextera XT Read 2 adaptor end.
Internal Tn5 fragments lack the tag site.

## Steps

1. 🟢 RT with SuperScript III (replacing SuperScript II of Smart-seq2) and the RT primer; ERCC spike-in processed in parallel.
2. 🟢 Template switching and PCR "according to the Smart-seq2 protocol" (ref 21; TSO and PCR primer not printed). Modelled as a warned boundary.
3. 🟢 Purification (method not stated); 0.1 ng cDNA tagmented with the Nextera XT kit.
4. 🟢 One end-capture PCR with the P5 and P7 index primers amplifies fragments from both cDNA ends ("fragments of both ends of the cDNA were selectively amplified"). Polymerase and cycle number are not stated.
5. 🟢 1.8× AMPure XP; E-Gel 2% SizeSelect, 200–300 bp; Bioanalyzer.
6. 🟢 Illumina X Ten, PE150.

Correction in this round: the earlier model split the library into separate TSS and TES PCR
branches with inferred primers. The Methods describe one PCR and one library; TSS and TES
reads are separated during processing (Methods lines 169–171).

## Final library and reads (🟡 derived)

- TES-end molecule: P5–[i5]–S5–ME–`GTGGTATCAACGCAGAGT`–N8–[cell barcode]–T30–VN–transcript 3′ end (antisense)–ME′–S7′–[i7]–P7′. Read 1 crosses the tag, UMI, cell barcode and oligo(dT). Read 2 enters from the Tn5 end toward the poly(A). This matches the processing rule (R1 with the 3′ tag, R2 with ≥10 A).
- TSS-end molecule: P5–[i5]–S5–ME–`GTGGTATCAACGCAGAGT`–`ACAT`–`GGG`–transcript 5′ end–ME′–S7′–[i7]–P7′. Read 1 begins with the printed TSS tag.
- The P7-side ME′ is the Nextera XT adaptor (kit not inspected) and is marked inferred. Standard Nextera primer sites are computed; the paper names only the platform and PE150.
- Results line 65 says the cell barcodes are "the same as sequencing indexes" added during amplification, while the RT primer also carries a 16-bp cell barcode. Both are kept as separate features.

## Discrepancy shown on the page

Results line 65 and the Fig. 1a legend (line 71) say the tag added to **both** ends contains
the UMI and cell identifier, and Fig. 1a draws CS2+TAG on the TSO. Methods line 169 identifies
TSS reads by `GTGGTATCAACGCAGAGTACATGGG` with no random or barcode positions, and only the RT
primer is printed with N8 and a cell barcode. The model draws UMI and cell barcode on the TES
end only. The full TSO is not printed, so whether a TSO-end UMI/cell identifier existed
remains open.

## Open questions

- Full Smart-seq2 TSO, preamplification primer and cycling (ref 21, not openly available).
- End-capture PCR polymerase and cycle number; index-read cycles; 16-bp cell-barcode list.
- Nextera XT adaptor original (kit document) not inspected.
