# PEM-seq — research notes

Yin et al., *Cell Discovery* (2019), [doi:10.1038/s41421-019-0088-8](https://doi.org/10.1038/s41421-019-0088-8).

## Assay target

The sequenced material is human genomic DNA around a bait double-strand break, quantifying repair outcomes. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a locus bait primer, a random-barcode bridge duplex and Illumina tagging primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining open-access article and Supplementary Methods.
- 🟢 Defining-paper Supplementary Table S11, which prints the bridge strands, P5/P7 tagging primers, I7 series, and locus-specific biotin/nested primers with modifications.
- 🟢 Liu et al., *STAR Protocols* (2022), [doi:10.1016/j.xpro.2021.101088](https://doi.org/10.1016/j.xpro.2021.101088), the authors' stepwise PEM-seq procedure and oligo table.

## Construct logic

1. 🟢 Cells are allowed to repair a recurrent bait DSB. Genomic DNA therefore already contains germline/perfect-rejoin, insertion, deletion, inversion, and translocation outcomes before library preparation.
2. 🟢 Extracted DNA is sonicated to 300–700 bp. A 5′-biotinylated locus primer then makes one extension product from each recoverable template. The single extension avoids the 80-cycle linear amplification used by LAM-HTGTS.
3. 🟢 Free biotin primer is removed, strands are denatured, and biotinylated extension products are immobilized on streptavidin beads.
4. 🟢 A 5′-phosphorylated bridge upper strand is ligated to the free 3′ end. Both bridge strands carry 3′ C7 amino blocks. The upper strand contains the pattern `NNNNTNNNANNNTNNNN`: fourteen random bases separated by fixed T/A/T bases.
5. 🟢 On-bead nested PCR uses a bait-side I5 primer and an I7-series bridge-side primer. For the RAG1A/B I5-R1 example, the five bases `AGCGA` between the Read 1 tail and common nested genomic primer act as the inline sample barcode. I7-B01 carries the six-base index `CGTGAT` in the ordering oligo; its reverse complement is present on the displayed library strand.
6. 🟢 A tagging PCR with the published P5-I5 and P7 oligos completes the flow-cell ends. The authors sequence 2×150 bp on HiSeq.

## Represented branch

The page shows the RAG1A/B bait with I5-R1 and I7-B01. Supplementary Table S11 provides alternative bait primers, four I5 barcode/nested variants for each locus, and B/R I7 series. These alternatives change target or sample identity without changing the molecular workflow.

## Read interpretation

Read 1 enters from the nested bait primer and covers the known bait sequence before reaching the repair junction. Read 2 enters from the bridge side and encounters the structured random molecular barcode. The model records the four random tracts as parts of one RMB group, excluding the fixed T/A/T separators from the 14-base UMI length.

PEM-seq quantifies outcomes recovered from the chosen bait side. A repair product that deletes either required primer-binding site is not recoverable, and translocation abundance also reflects three-dimensional contact probability; the library should not be read as a direct frequency measurement for every off-target locus.
