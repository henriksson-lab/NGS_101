# sci-Hi-C

## 1. What it is

Ramani et al.'s combinatorial-indexing Hi-C method, in which the intersection of two plate barcodes identifies a nucleus without physical isolation into one well at the start.

## 2. Sources and evidence

- 🟢 Ramani et al., *Nature Methods* 2017, DOI [10.1038/nmeth.4155](https://doi.org/10.1038/nmeth.4155); full text [PMC5330809](https://pmc.ncbi.nlm.nih.gov/articles/PMC5330809/).
- 🟢 DpnII-digested nuclei receive biotinylated indexed bridge adaptors in the first 96-well plate.
- 🟢 Nuclei are pooled, proximity-ligated, redistributed, and receive a second index through custom Illumina Y adaptors before streptavidin enrichment and PCR.
- 🟢 The detailed Supplementary Methods, steps 6–12 (`nmeth.4155_41592_2017_BFnmeth4155_MOESM254_ESM.pdf.txt`, lines 585–761), reverse-crosslink after second-well dilution, clean with AMPure, digest with **MseI + HaeIII** in the wells, end-repair, dA-tail, ligate the indexed Y adaptors, then pool and perform biotin pull-down.
- 🟡 The main-text overview (`nmeth.4155_PMC5330809.xml.txt`, line 62) instead lists pooling, streptavidin purification, then restriction digestion. That order conflicts with the detailed steps. The public graph follows the Supplementary Methods and its page states both source locations; whether the overview compresses or reverses the actual laboratory order is unresolved.
- 🟢 Supplementary Data print the bridge-adapter oligos and all 96 Y-adapter pairs. The
  model uses A1 representatives for both rounds: BC1 `GAGCTCGA` in the
  `GATC–BC1–GAATT` biotinylated bridge strand, and BC2 `TGACCTTG` in the SeqShort HI/LO
  Y adapter. The published internal biotin-dT, 5′ phosphate, terminal phosphorothioate,
  DpnII-compatible end, EcoRI overhang and dA/dT ligation end are retained.

## 3. Distinguishing chemistry

The first barcode is internal to the contact molecule. BC2 is the first eight bases of
both genomic reads (followed by three fixed bases), not an Illumina index read. Both
barcode rounds are required, so they are never collapsed into a generic dual-index library.

## 4. Bounded source pass — 2026-10-17

Read published article archived manuscript PMC5330809 JATS overview60–64,
Online Methods86–122 (fixation, sciHi-C, barcode design, processing). Read
published Springer supplement `MOESM254_ESM.pdf.txt` protocol sections2–5
(text435–590), amplification/cleanup and primer table755–846; inspected
rendered original PDF p22 to verify the reverse-primer text. Read published
Supplementary Spreadsheet1 representative bridge rows3/101 and Y-adapter
rows198/294, with worksheet headings. The manifest identifies the workbook
as the publisher/PMC duplicate. The superseded 065052 preprint was excluded.
Ran the bundle scraper; used hits were checked at original passages.

🟢 Supplementary protocol section4 requires PNK phosphorylation after pooling
and washing the first indexed nuclei and before in-situ ligation. The graph
now includes that step, names the T4 ligases and KAPA Hi-Fi PCR mix, and
records two 250-base reads from protocol13.10 and the main overview. BC2's
opposite-end feature no longer includes the adjacent fixed CGT-derived bases.

🟢 The final amplification uses Illumina F
`AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
and Barcoded Illumina R
`CAAGCAGAAGACGGCATACGAGAT[10 bp barcode]TGACTGGAGTTCAGACGTGTGCT`.
The ten-base PCR barcode is separate from BC1 and BC2; its placeholder is
now represented as a sample index. No acquired index read is inferred from
its presence. The graph gives Illumina F its source-derived exact input state.

🟡 A published junction difference remains unresolved: the reverse PCR primer
on rendered supplement p22 begins its template arm `TGACTGG...` after the
barcode, while the workbook SeqShort_BC1_LO end `...TCCAGTCAC` corresponds
to `GTGACTGG...`. The canonical complete Read-2 arm formerly drawn without
qualification is now inferred. The public page names both locations/values,
the model choice and remaining uncertainty; no source oligo is silently
repaired, and no exact PCR product is claimed from the printed pair. The
existing main-overview versus detailed-protocol digestion-order notice is
retained; steps6–12 were not newly re-read in this pass.

Deferred: protocol sections6–12 outside this pass; complete barcode families
and author code barcode-list originals; Protocol Exchange identity/version and
predecessor/vendor originals; internal biotin geometry and phosphorylation
propagation; carrier plasmid, wash and final recovery/QC states; right-arm
PCR junction and sample-specific ten-base index/acquisition. This is one
capped pass, not complete molecular validation or accepted source absence.
