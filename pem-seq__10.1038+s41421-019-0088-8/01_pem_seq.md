# PEM-seq — research notes

Yin et al., *Cell Discovery* (2019), [doi:10.1038/s41421-019-0088-8](https://doi.org/10.1038/s41421-019-0088-8).

## Assay target

The sequenced material is human genomic DNA around a bait double-strand break, quantifying repair outcomes. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a locus bait primer, a random-barcode bridge duplex and Illumina tagging primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining open-access article and Supplementary Methods.
- 🟢 Defining-paper Supplementary Table S11, which prints the bridge strands, P5/P7 tagging primers, I7 series, and locus-specific biotin/nested primers with modifications.
- 🟢 Liu et al., *STAR Protocols* (2022), [doi:10.1016/j.xpro.2021.101088](https://doi.org/10.1016/j.xpro.2021.101088), the authors' stepwise PEM-seq procedure. The running text extracts normally; **Table 2 (the oligo table) does not**. It is a vector object: `pdftotext` returns only its caption, and the separately registered `fx6.gif` is a 100×160-px thumbnail. Read it by rendering page 11 of `PMC9019705_PMC9019705.pdf` at ≥200 dpi, e.g.
  `pdftoppm -r 400 -f 11 -l 11 -png $CHEM_DATA/sources/pem-seq__10.1038+s41421-019-0088-8/PMC9019705_PMC9019705.pdf out`.
  Table 1 (bait DSBs, page 10) has the same problem. The colour coding in Table 2 is load-bearing: the index bases are
  printed in red and the RMB tract orange and underlined, so a text-only reading cannot recover either boundary.

## Construct logic

1. 🟢 Cells are allowed to repair a recurrent bait DSB. Genomic DNA therefore already contains germline/perfect-rejoin, insertion, deletion, inversion, and translocation outcomes before library preparation.
2. 🟢 Extracted DNA is sonicated to 300–700 bp. A 5′-biotinylated locus primer then makes one extension product from each recoverable template. The single extension avoids the 80-cycle linear amplification used by LAM-HTGTS.
3. 🟢 Free biotin primer is removed, strands are denatured, and biotinylated extension products are immobilized on streptavidin beads.
4. 🟢 A 5′-phosphorylated bridge upper strand is ligated to the free 3′ end. Both bridge strands carry 3′ C7 amino blocks. The upper strand contains the pattern `NNNNTNNNANNNTNNNN`: fourteen random bases separated by fixed T/A/T bases.
5. 🟢 On-bead nested PCR uses a bait-side I5 primer and an I7-series bridge-side primer. The inline sample barcode is the **four** bases between the Read 1 tail and the locus-specific nested primer: `AGCG` (I5-R1), `GCCT` (R2), `AGGA` (R3), `TCAG` (R4). For the RAG1A/B example the nested primer is therefore `AGGATCTCACCCGGAACAGC` (20 nt), and the oligo reads `ACTCTTTCCCTACACGACGCTCTTCCGATCT` + `AGCG` + `AGGATCTCACCCGGAACAGC`. I7-B01 carries the six-base index `CGTGAT` in the ordering oligo; its reverse complement is present on the displayed library strand.
6. 🟢 A tagging PCR with the published P5-I5 and P7 oligos completes the flow-cell ends. The authors sequence 2×150 bp on HiSeq.

## Why the inline I5 barcode is 4 nt, not 5

The concatenated oligo cannot distinguish the parses on its own — `…CCGATCT` + `AGCGA` +
`GGATCTCACC…` and `…CCGATCT` + `AGCG` + `AGGATCTCACC…` are the same string — so this was
resolved two independent ways.

- 🟢 **Within Table S11.** The table prints I5-R1…R4 for three loci. Under the 4-nt parse
  each locus keeps one shared nested primer across all four barcodes: RAG1A/B
  `AGGATCTCACCCGGAACAGC`, MYC1 `AGAGCCATTCTCTGGCTCAG`, RAG1C
  `TATCAATATCCCACTGATGTATC`. Under the 5-nt parse the barcode set would have to be
  `AGCGA/GCCTA/AGGAA/TCAGA` at RAG1A/B and MYC1 but `AGCGT/GCCTT/AGGAT/TCAGT` at RAG1C —
  i.e. a sample barcode whose last base depends on which locus is being amplified. The
  4-nt parse is the only one that makes the barcode a property of the sample.
- 🟢 **STAR Protocols Table 2.** The rendered table (page 11, see above) prints the index
  bases in red for rows I5-Nested-1…4: exactly `AGCG`, `GCCT`, `AGGA`, `TCAG`, each
  following the complete `…CTTCCGATCT` tail. Row I5-Nested-1-Myc6 shows the same red
  `AGCG` followed by `GGAAACCAGAGGGAATCCTC`, the Myc6 nested primer, in black.

The library's bases are unchanged by this correction; only the feature boundary, the
segment length and the Read 1 cycle span move by one base.

## Cross-source conflict: the RMB length

🟡 The two sources disagree about the bridge upper strand, and the model deliberately
follows the defining paper.

| Source | Bridge-adapter-upper random region | Random bases |
|---|---|---|
| Defining paper, Suppl. Table S11 | `NNNNTNNNANNNTNNNN` | 14 |
| STAR Protocols Table 2 | `NNNTNNNNTNNN` | 10 |

The STAR table contradicts its own prose, which says the adapter contains "a 14-bp RMB"
(and its Step 7 narrative repeats it); the defining paper also compares 14-bp with 7-bp
RMB libraries in its supplementary figures. Table S11 is internally consistent with both
papers' text, so `PEM_BA_UP` keeps the 14-base, T/A/T-separated form. The STAR Table 2
printing is recorded here as an error in the later paper, not as a variant to model; a
comment in `tools/chem.py` says so at the constant.

## Represented branch

The page shows the RAG1A/B bait with I5-R1 and I7-B01. Supplementary Table S11 provides alternative bait primers, four I5 barcode/nested variants for each locus, and B/R I7 series. These alternatives change target or sample identity without changing the molecular workflow.

## Read interpretation

Read 1 enters from the nested bait primer and covers the known bait sequence before reaching the repair junction. Read 2 enters from the bridge side and encounters the structured random molecular barcode. The model records the four random tracts as parts of one RMB group, excluding the fixed T/A/T separators from the 14-base UMI length.

PEM-seq quantifies outcomes recovered from the chosen bait side. A repair product that deletes either required primer-binding site is not recoverable, and translocation abundance also reflects three-dimensional contact probability; the library should not be read as a direct frequency measurement for every off-target locus.
