# TEA-seq — source notes

Defining paper: Swanson et al., *eLife* (2021), [doi:10.7554/eLife.63632](https://doi.org/10.7554/eLife.63632).

## Sources inspected

The original article Methods (TEA-seq and sequencing), Supplementary file 5,
10x CG000338 Rev A p.80, CG000373 Rev A pp.2–3 and the authors' detailed
bench protocol step 80 were read after running the oligo scraper. Independent
located checks and remaining source work are recorded in the source-review ledger.
The detailed protocol and vendor sequencing note were followed from the paper's
methods and sequencing instructions.

## Construction and reaction order

- 🟢 TEA-seq combines TotalSeq-A antibody tags with the Multiome workflow.
  Tagmentation and GEM partitioning precede shared reverse transcription/barcoding,
  cleanup and preamplification with added ADT-Rev-AMP. An aliquot enters ATAC indexing;
  another enters cDNA/ADT amplification and size separation. GEX preparation then
  includes fragmentation, end repair, A-tailing, adapter ligation and TT indexing.
- 🟢 CG000338 Rev A p.80 prints distinct gel-bead scaffolds: partial TruSeq Read 1,
  cell barcode N16, UMI N12 and T30VN for RNA; P5, cell barcode N16 and CGCGTCTG
  for ATAC. It also prints the complete dual-index GEX and single-sample-index ATAC
  products. Each is declared locally; a shared endpoint recipe previously assigned
  inappropriate handles/index architecture. That was a repository modelling error.
- 🟢 TEA Methods explicitly assign Supplementary file 5's SI-PCR-Oligo and ADT-i7
  to ADT indexing, despite the workbook's ICICLE heading. Verbatim oligos:

  - SI-PCR-Oligo: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`
  - ADT-Rev-AMP: `CCTTGGCACCCGAGAATTCC`
  - ADT-i7: `CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCC*A`

  The asterisk is the printed terminal phosphorothioate bond. All three complete
  primer base sequences are checked directly against the extracted original table.
- 🟡 The ADT endpoint is derived from these primers and the printed bead scaffold;
  it is distinct from the GEX dual-index and Nextera ATAC endpoints. Commercial
  precursor duplexes, named primer mixtures and ligation-ready adapter geometry
  have not all been assembled in this model. Those are model limitations, not
  accepted claims of missing publication information.

- 🟢 2026-10-17 bounded re-read (article Methods lines 745–767; protocols.io
  TEA-seq V.2, DOI 10.17504/protocols.io.bqagmsbw, steps 45–49, 61–64): the GEM
  incubation is 37°C 45 min then 25°C 30 min followed by Quenching Agent; the
  preamplification (7 cycles) is followed by a 2.0× SPRIselect cleanup, now a
  separate graph step; the GEX fragmentation, end repair and A-tailing are one joint
  reaction (32°C 5 min, 65°C 30 min); ATAC indexing is 9 cycles (0.6×/1.6× SPRI),
  GEX TT indexing 14 cycles (0.6×/0.8× SPRI) and ADT indexing 15 cycles with KAPA HiFi.
- 🟢 Source discrepancy shown on the page: article Methods use 5 µL Pre-Amp
  Primers; protocols.io V.2 step 46 lists 4 µL. The model does not depend on the
  volume. The article's "dual-sided 1.6×" ADT cleanup names one ratio; protocols.io
  step 64 gives a single 1.6× cleanup.

## Sequencing

- 🟢 The raw joint run is R1/I1/I2/R2 = 50/10/16/90. The original article's masks
  retain GEX 28/10/10/90, ATAC 50/8/16/50 and ADT 28/8/no I2/90.
- 🟢 The author bench protocol step 80 explicitly identifies forward indexing.
  CG000373 Rev A establishes that forward ATAC I2 enters from P5 and reads the
  cell barcode directly. Reverse indexing traverses the eight-base spacer and
  requires a different cycle recipe; it must not be silently substituted.
- 🟡 Exact loaded ADT R2/I1 oligos have not been verified. The page shows the
  source-declared reads with that limitation, without inventing primer sites or
  feature-cycle coordinates. Selected TotalSeq-A and bead-barcode members remain
  unassigned; source follow-up continues.

## Apparent vendor-guide inconsistency

CG000338 Rev A p.80 prints the GEX prefix before i5 as
`AATGATACGGCGACCACCGAGATCT` (25 bases). CG000317 Rev E p.72 prints
`AATGATACGGCGACCACCGAGATCTACAC` (29 bases), in both the TT forward primer
and GEX product. Original page images confirm the difference. Both guides name
TT kit PN-1000215 and plate PN-3000431, so an established reagent variant cannot
be inferred. TEA's model retains its explicitly cited CG000338 Rev A layout;
SEC-seq retains its separately identified CG000317 layout. Both public pages
state the conflict. The physical prefix in the historical experiments remains
unverified.
