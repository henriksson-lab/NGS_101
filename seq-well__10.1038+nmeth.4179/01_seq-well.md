# Seq-Well — research-note pointer

Seq-Well is documented from Gierahn *et al.*, 2017,
doi:[10.1038/nmeth.4179](https://doi.org/10.1038/nmeth.4179).

The detailed source record was originally written jointly with Drop-seq and remains
available as provenance in the
[combined Drop-seq / Seq-Well research note](../drop-seq__10.1016+j.cell.2015.05.002/01_drop-seq-seq-well.md#seq-well-from-the-master-protocol-moesm238-per-array).
The public schematic uses the Seq-Well paper’s Supplementary Table 1 and Master Protocol
as its authoritative molecular sources.

## Bounded source pass 2026-10-17 (round breadth-first-18)

Read: Supplementary Table 1 (MOESM237 text lines 260–300, five oligos), the Seq-Well Master
Protocol V1.0, 2/13/2017 (MOESM238 text lines 470–560 sealing/lysis/hybridization/bead removal,
560–760 RT/Exo I/WTA/SPRI, 760–860 Nextera XT tagmentation, library PCR, sequencing), and the
Supplementary Online Methods (NIHMS843531 supplement 1, lines 17–23, which defer to the protocol).

- 🟢 Oligos re-checked against Supplementary Table 1: bead oligo (linker–T7–SMART–AC–J12–N8–T30),
  TSO `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`, SMART PCR primer, New-P5-SMART hybrid with
  phosphorothioates before the final A and C, and the custom Read 1 primer. No transcription change.
- 🟢 Added to the graph: membrane sealing (37 °C, 30 min, clamped); guanidine thiocyanate lysis buffer
  with Sarkosyl and β-mercaptoethanol (20 min) and hybridization buffer (40 min); membrane removal and
  bead spin-off; Maxima H− RT with TSO (room temperature 30 min, then 52 °C 90 min); TE-TW/TE-SDS
  washes; Exo I (37 °C, 5 + 45 min); KAPA HiFi WTA with 1,500–2,000 beads per 50 µL, 4 cycles then
  9–12 cycles (13 total cycles for cell lines, 16 for primary cells); 0.6× SPRI; Nextera XT
  tagmentation of 750 pg (55 °C 5 min, NT neutralization); Nextera PCR mix with New-P5-SMART and
  N70X, 12 cycles; 0.6× SPRI (650–750 bp average; 420–800 bp sequenced).
- 🟢 Reads: Read 1 20 bp (custom primer), Index 1 8 bp (only when multiplexing), Read 2 50 bp.
- 🟡 Model correction: tagmented fragments without the bead end were previously "discarded" at
  tagmentation. No separation step removes them. Supplementary Table 1 note 4 says the P5-SMART primer
  selectively primes the bead-specific SMART sequence, so these fragments are not exponentially
  amplified. They are now an explicitly unresolved fate.

Open/deferred: the Macosko 2015 bead-synthesis details were not read. The Nextera XT guide (Rev E) and the NextSeq custom-primer
guide 15057456-01 are linked but unavailable (MANIFEST manual rows). The exact N70X and transposome
adapter sequences come from Illumina documents that were not read in this pass.

## Protocol Exchange protocol, first chemistry reading pass (2026-10-17)

Source: the Protocol Exchange protocol cited by the published article (reference 21 after the
publisher's erratum), doi:10.1038/protex.2017.006a, posted 14 March 2017. Two copies are on disk:
the protocols.io import (`protocolsio_…_v1.pdf.txt`, printed "Version 1", "Imported from Protocols
Exchange", 48 pp.) and the Exchange's own 30-page export, attached to that record
(`…_attachment_31_syjpcicn7.pdf.txt`, "Posted Date: March 14th, 2017"). One reader so far; an
independent second pass is still open.

- 🟢 Oligos (import lines 241–265; export lines 80–100): TSO `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`;
  "IS PCR Primer" `AAGCAGTGGTATCAACGCAGAGT` (the SMART PCR primer of Supplementary Table 1);
  P5-SMART PCR hybrid oligo `AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`;
  custom Read 1 primer `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`. The bases equal Supplementary
  Table 1. The bead oligo is not printed; beads are ordered as "mRNA Capture Beads (Chemgenes Cat No.
  MACOSKO-2011-10)" (import line 274). `tools/scrape_primers.py` finds only these oligos.
- 🟡 Modification difference: Supplementary Table 1 prints the P5-SMART primer as `…CAGAGT*A*C`
  (two phosphorothioates). The import prints no `*`. The export prints `…CAGAGT**A**C` in a reagent
  list that shows raw Markdown marks throughout (rendered page 3 checked), so the marks were probably
  lost in formatting. The model keeps Supplementary Table 1's modifications, and the page shows this
  difference.
- 🟢 Reaction order and conditions (import steps 80–131, lines 803–1145) agree with the Master
  Protocol already modelled. RT mix per array: Maxima H− RT, Ficoll PM-400, dNTPs, RNase inhibitor,
  5 µL 100 µM TSO; room temperature 30 min, then 52 °C 90 min. TE-TW/TE-SDS washes. Exo I 37 °C,
  5 + 45 min. WTA with KAPA HiFi and 0.4 µL 100 µM SMART primer, 1,500–2,000 beads per 50 µL:
  95 °C 3 min; 4 cycles of 98/65/72 °C; 9 cycles of 98/67/72 °C; 72 °C 5 min. The note gives 13 total
  cycles for cell lines and 16 for primary cells. 0.6× SPRI. Nextera XT tagmentation of 750 pg,
  55 °C 5 min, NT neutralization. Library PCR with 15 µL Nextera PCR mix, 1 µL 10 µM P5-SMART and
  1 µL 10 µM N70X: 95 °C 30 s; 12 cycles of 95/55/72 °C; 72 °C 5 min. 0.6× SPRI. Sequencing:
  Read 1 20 bp (custom primer), Read 2 50 bp, Index 8 bp only when multiplexing. It links the
  Illumina NextSeq (15057456-01) and MiSeq (15041638-01) custom-primer guides.
- No model change was needed. The page now shows the P5-SMART modification difference.

