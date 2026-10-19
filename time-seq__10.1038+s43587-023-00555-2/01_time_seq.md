# TIME-seq

Griffin *et al.*, “TIME-seq reduces time and cost of DNA methylation measurement for
epigenetic clock construction,” *Nature Aging* (2024),
[doi:10.1038/s43587-023-00555-2](https://doi.org/10.1038/s43587-023-00555-2).

## Sources read

- 🟢 The author manuscript (PMC11332592), especially “TIME-seq library preparation”
  (`PMC11332592_PMC11332592.xml.txt` lines 189–201), “Sequencing” (line 205), the
  Results overview (line 117) and Extended Data Figure 2 legend (line 281).
- 🟢 The publisher Supplementary Information (`43587_2023_555_MOESM1_ESM.pdf.txt`),
  sections II (adaptor annealing and transposome loading, lines 79–100) and IV
  (sequencing, lines 152–163).
- 🟢 Supplementary Table 5 (`43587_2023_555_MOESM7_ESM.xlsx.txt`), including the Tn5
  adaptors (TS301–TS364), hybridization blockers, indexed PCR primers and custom
  sequencing primers.
- The bioRxiv v3 preprint is superseded by the publication and is not used for values.

## Construct-changing steps

1. 🟢 Each 100-ng genomic-DNA sample is tagmented for 15 min at 55 °C with a uniquely
   indexed TIME-seq transposome. Adaptor B carries a five-base internal sample barcode;
   methylated adaptor A and indexed adaptor B are separately annealed to the
   5′-phosphorylated reverse mosaic-end oligo when transposomes are assembled.
   100 ng unmethylated lambda DNA is tagmented with each pool as a bisulfite-conversion
   control. STOP buffer (guanidine thiocyanate) is added and incubated 5 min at 55 °C.
2. 🟢 Samples are combined by pool with Zymo DNA binding buffer and purified on a Clean &
   Concentrator-25 column. Klenow fragment (3′→5′ exo−) in NEBuffer 2
   performs 30 min of end repair at 37 °C using a dNTP mix in which 5-methyl-dCTP replaces
   dCTP; a second column and SpeedVac concentrate the pool to 5 µL. 🟡 Methyl-dCTP
   protects newly copied adapter/gap cytosines against bisulfite conversion.
3. 🟢 DNA with Cot-1 DNA and the A/B blocking oligos is denatured 5 min at 95 °C;
   biotinylated RNA baits (750 ng) hybridize to target loci for 4 h at 65 °C. RNA:DNA hybrids are captured on streptavidin beads,
   stringently washed and eluted as single-stranded DNA with NaOH.
4. 🟢 Captured DNA is converted with the EpiTect Fast Bisulfite kit, then amplified for
   19 cycles with Q5U and indexed P5/P7 primers, then cleaned with 1.8× SPRI beads. The PCR indices distinguish pools; the
   five-base transposome index distinguishes samples within a pool.

## Oligos and read layout

🟢 Supplementary Table 5 publishes 64 indexed B adaptors (TS301–TS364). The schematic
uses TS301 (`AAAGA`) as one explicit representative and N501/N701 as representatives of
the published N501–N506 (P5/i5) and N701–N704 (P7/i7) PCR primers. 🟢 Each N70x primer is
P7, the 8-nt i7 pool index, the 23-nt custom sequencing handle `AAGCAGTGGTATCAACGCAGATC`
and the 14-nt B-adaptor prefix `TGGGTGGAGGGTGG`; the custom Read 2 primer is handle + prefix
and the custom Index 1 primer is its reverse complement (Table 5, sheets “Indexed PCR
Primers” and “Sequencing Primers”). It preserves `/5Phos/` and `/iMe-dC/` exactly as written
in the source.

🟡 Published discrepancy (computed from Supplementary Table 5, sheet “Hyb Blocking Oligos”,
checked in the original workbook): the “Adaptor A Blocking Oligo”
`CTGTCTCTTATACACATCTCCGAGCCCACGAGAC` is the reverse complement of the Nextera Tn5ME-B arm
`GTCTCGTGGGCTCGG` plus the mosaic end, not of TIME-seq methylated adaptor A (arm
`TCGTCGGCAGCGTC`); only its 19-nt mosaic-end portion pairs with adaptor A. The B blocker is
the exact reverse complement of adaptor B (`HHHHH` opposite `DDDDD`). The model keeps the
published sequence as a hybridization input and draws no blocker duplex; which sequence was
intended is unresolved. Minor table inconsistencies, also shown on the page: the P5 structure
row has a 7-N i5 placeholder under an “(8-nt)” header while all listed indexes have 8 nt, and
the adaptor B structure row labels the 14-nt prefix “Adaptor A”.

🟢 TIME-seq requires its published custom Read 2 and Index Read 1 primers. Read 1 and
Index 2 use the corresponding Nextera sites. The economical MiSeq/NextSeq program is
Read 1 145–153 cycles, i7 8 cycles, optional i5 8 cycles and Read 2 5 cycles. Thus the
short Read 2 reports the internal Tn5 sample barcode.

## Derived observations and open questions

- 🟡 Computed from Supplementary Table 5: none of the 64 five-base TS301–TS364 barcodes
  contains C, nor does the 14-base B prefix, so the Tn5 sample barcode read by the
  5-cycle Read 2 is unaffected by bisulfite. Adaptor B is listed without `/iMe-dC/`;
  its mosaic-end segment contains one C. The schematic draws the product whose B-side
  mosaic end and barcode complement were copied by Klenow with 5-methyl-dCTP; the
  complementary strand-of-origin product (B adaptor top strand converted at that C) is
  not drawn. Neither source discusses this strand.
- 🟡 Tagmentation with equal A and indexed B adaptor loading (1:1:2 A:B:glycerol,
  Supplementary Note II) is inferred also to yield fragments carrying the same adaptor at
  both ends. 🔴 Neither the Methods nor Supplementary Note II removes or quantifies them;
  the graph keeps them as an unresolved terminal rather than a physical discard (no
  removal step precedes PCR; whether they are captured, amplified or sequenced is not
  reported).
- Sample B in the graph uses TS302 only to show two distinct published barcodes; actual
  barcode-to-sample assignments are study specific.
- Bait sequences, bait production (Supplementary Note III) and Supplementary Tables 2–4
  were not re-read in this capped pass.

## Scope

The page documents library construction and sequencing geometry, not epigenetic-clock
training, benchmark performance or bait-design results. Bait sequences are target-set
specific and are not reproduced as “key oligos”; the adapter-blocking oligos are.
