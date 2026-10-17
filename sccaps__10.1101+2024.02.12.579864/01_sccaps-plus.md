# scCAPS+

## Sources read

- 🟢 Chen et al., bioRxiv preprint [doi:10.1101/2024.02.12.579864](https://doi.org/10.1101/2024.02.12.579864), the shared defining record.
- 🟢 Chen et al., version of record, *Genome Biology* 26:244 (2025), [doi:10.1186/s13059-025-03708-1](https://doi.org/10.1186/s13059-025-03708-1), especially Methods and Additional file 5.

## Molecular path

🟢 Single cells are FACS-sorted, lysed and tagmented with one of 96 preassembled i5/i7-barcoded Tn5 complexes. After SDS stripping, neutralization and 72 °C gap filling, 96 cell-barcoded reactions are pooled.

🟢 ACT+BF4− selectively oxidizes 5hmC to 5fC. Pinnick oxidation with sodium chlorite and 2-methyl-2-butene converts 5fC to 5caC. Pyridine–borane reduces 5caC to DHU, which PCR reads as T. 5mC and unmodified C remain C, making scCAPS+ selective for 5hmC.

🟢 Nextera XT index primers amplify the converted pool. The paper uses 120-bp paired-end sequencing with custom NovaSeq primers; Additional file 5 Tables S3–S4 define the transposome and sequencing oligos.

## Oligo source review

🟢 The shared paper's Additional file 5, Table S3 prints eight `/5Phos/` i5
and twelve `/5Phos/` i7 Tn5 transfer oligos and the shared reverse mosaic-end
strand `/5Phos/CTGTCTCTTATACACATCT`. For example, i5-1 is
`/5Phos/TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG`
and i7-1 is
`/5Phos/GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`.
The 8 × 12 pairings index 96 cells. These are transposome oligos, not PCR primers.

🟢 Table S4 prints custom Read 1
`GCGATCGAGGACGGCAGATGTGTATAAGAGACAG`, Read 2
`CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`, Index 1
`CTGTCTCTTATACACATCTGAGGCGGAGACGGTG`, and Index 2
`CTGTCTCTTATACACATCTGCCGTCCTCGATCGC` sequencing primers.
Table S1's three oligos are a 5hmC spike-in control, not genomic adapters.

🟡 Derived construct (model 2026-10-17): every i5 member is S5 `TCGTCGGCAGCGTC` +
spacer `TCCACGC` + 8-nt barcode + C15 `GCGATCGAGGACGGC` + ME; every i7 member is
S7 `GTCTCGTGGGCTCGG` + spacer `CTGTCCCTGTCC` + 8-nt barcode + D15
`CACCGTCTCCGCCTC` + ME. The Table S4 primers are exactly C15+ME (Read 1),
D15+ME (Read 2), ME_RC+revcomp(D15) (Index 1) and ME_RC+revcomp(C15) (Index 2);
`tools/selftest.py` checks the composed oligos against the Additional file 5
workbook. Read 1/Read 2 therefore start in the insert, and the custom index reads
cover the Tn5 i7/i5 cell barcodes rather than the Nextera XT sample indices further
out. 🟡 The XT P5/P7 ends and 8-nt sample indices are drawn INFERRED: the paper names
Nextera XT Index Kit FC-131-1001 but prints neither its primers nor the indices used.
🔴/open: index-read cycle counts are not stated (8 is inferred from barcode length);
the published text gives only 120-bp paired-end reads (Results, JATS line 90).

## Reaction details (published Methods, JATS lines 132, 144, 150–158)

🟢 Sort into 4 µl lysis buffer (16.67 mM TAPS-NaOH pH 8.5, 8.33 mM MgCl2, 0.67× NEB
buffer 4, 0.13% Triton X-100, 1 µg Qiagen protease); lyse 50 °C 3 h, 75 °C 30 min,
80 °C 15 min. Tagment 50 °C 15 min in 1× TAPS buffer, 6% PEG 3350, 1 µl 1×
preassembled Tn5 (10 µl). Strip with 1.11 µl 1.0% SDS, neutralize with 1.24 µl 5.0%
NP-40 or Triton X-100, add 12.35 µl 2× Phusion HF master mix, 72 °C 5 min. Pool 96
reactions with 200 ng 2-kb filler and 0.1% spike-ins, purify, convert, then amplify
with KAPA HiFi HotStart Uracil+ and Nextera XT index primers (98 °C 45 s; 13 × 98 °C
10 s / 60 °C 15 s / 72 °C 1 min; 72 °C 3 min) and select 500–700 bp with 0.5×–0.25×
AMPure XP. Transposomes: 8 × 12 Tn5MEDS-i5/i7 plus Tn5MEDS-Rev annealed in 1× TE
(95 °C 5 min, −1 °C per 12 s), assembled 1 h RT with in-house Tn5, diluted sevenfold.

## Bounded check 2026-10-17 (round breadth-first-16)

Published Genome Biology 2025 (PMC12359873) Methods and Additional file 5 are the
authority; the bioRxiv v1 preprint is superseded. Corrections: the schematic now
uses the Table S3 transposon arms and Table S4 custom primers instead of a generic
Nextera shell; the earlier note's i7-1 transcription had dropped one C (repository
error, not a published one); lysis/strip/purification/size-selection steps added.
The published scCAPS+ Methods (JATS line 154) contain no βGT step: 5hmC is
oxidized by ACT+BF4− (25 °C 16 h), then Pinnick NaClO2/2-methyl-2-butene (25 °C
16 h), then borane-reduced; the catalogue's “βGT protection” summary is a
repository error (reported for catalogue correction).
Deferred: Nextera XT Index Kit original; Amini et al. 2014 (ref. 28) and Picelli
et al. 2014 (ref. 27) transposome predecessors; TAPS/CAPS+ conversion predecessors
(refs. 15, 17); index-read cycle counts.

## Protocol boundary

This is a separate schematic from scTAPS despite the joint paper. The cell-barcoding front end is shared, but scCAPS+ has two selective chemical oxidations and a 5hmC-only readout; scTAPS instead uses mTET1 and reports 5mC+5hmC.
