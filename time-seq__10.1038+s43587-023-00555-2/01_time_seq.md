# TIME-seq

Griffin *et al.*, “TIME-seq reduces time and cost of DNA methylation measurement for
epigenetic clock construction,” *Nature Aging* (2024),
[doi:10.1038/s43587-023-00555-2](https://doi.org/10.1038/s43587-023-00555-2).

## Sources read

- 🟢 The author manuscript, especially “TIME-seq library preparation,” “Sequencing,”
  and Extended Data Figure 2.
- 🟢 The detailed Supplementary Information, sections II and IV.
- 🟢 Supplementary Table 5, including the Tn5 adaptors, hybridization blockers,
  indexed PCR primers and custom sequencing primers.

## Construct-changing steps

1. 🟢 Each 100-ng genomic-DNA sample is tagmented for 15 min at 55 °C with a uniquely
   indexed TIME-seq transposome. Adaptor B carries a five-base internal sample barcode;
   methylated adaptor A and indexed adaptor B are separately annealed to the
   5′-phosphorylated reverse mosaic-end oligo when transposomes are assembled.
2. 🟢 Samples are stopped, combined by pool and purified. Klenow fragment (3′→5′ exo−)
   performs 30 min of end repair at 37 °C using a dNTP mix in which 5-methyl-dCTP replaces
   dCTP. This protects newly copied adapter/gap cytosines against bisulfite conversion.
3. 🟢 Blocking oligos suppress adapter capture while biotinylated RNA baits hybridize to
   target loci for 4 h at 65 °C. RNA:DNA hybrids are captured on streptavidin beads,
   stringently washed and eluted as single-stranded DNA with NaOH.
4. 🟢 Captured DNA is converted with the EpiTect Fast Bisulfite kit, then amplified for
   19 cycles with Q5U and indexed P5/P7 primers. The PCR indices distinguish pools; the
   five-base transposome index distinguishes samples within a pool.

## Oligos and read layout

🟢 Supplementary Table 5 publishes 64 indexed B adaptors (TS301–TS364). The schematic
uses TS301 (`AAAGA`) as one explicit representative and N501/N701 as representative PCR
indices. It preserves `/5Phos/` and `/iMe-dC/` exactly as written in the source.

🟢 TIME-seq requires its published custom Read 2 and Index Read 1 primers. Read 1 and
Index 2 use the corresponding Nextera sites. The economical MiSeq/NextSeq program is
Read 1 145–153 cycles, i7 8 cycles, optional i5 8 cycles and Read 2 5 cycles. Thus the
short Read 2 reports the internal Tn5 sample barcode.

## Scope

The page documents library construction and sequencing geometry, not epigenetic-clock
training, benchmark performance or bait-design results. Bait sequences are target-set
specific and are not reproduced as “key oligos”; the adapter-blocking oligos are.
