# Stereo-seq — research notes

## 1. What it is

Stereo-seq combines a patterned array of clonally amplified DNA nanoballs with in-situ poly(A) RNA capture. Each DNB's coordinate identifier is decoded before tissue placement; captured cDNA retains that identifier and a molecular identifier.

## 2. Sources read

- 🟢 Chen et al., *Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays*, Cell 185, 1777–1792.e21 (2022), DOI [10.1016/j.cell.2022.04.003](https://doi.org/10.1016/j.cell.2022.04.003), [published PDF at Cambridge](https://www.repository.cam.ac.uk/bitstreams/1d4e7717-3e50-4448-8bda-a35871dfbfef/download).
- 🟢 The published STAR Methods reagent table (local text lines 1365–1472) prints DNB precursor, fixed oligo, splint, CID sequencing primer, DNB F/R, capture oligo, modified TSO, cDNA PCR primer, Stereo-seq library F/R and splint, read1/read2 and MDA primer. Methods lines 1676–1725 give the chip, reverse transcription, Tn5 and library amplification order; lines 1748–1755 assign CID/MID/read roles.
- 🟢 The 2021 author preprint v3, Methods lines 3900–3997, independently prints the chip and capture oligos. Its library section names later primers without printing their bases.
- 🟢 The paper-linked [MOSTA detailed protocol](https://db.cngb.org/stomics/assets/html/stereo.seq.html) includes a chip/library procedure and oligo table. Its full post-ligation capture-probe construct is distinct from the short phosphorylated capture oligo printed by Cell.

## 3. Construct-changing steps

1. Generate coordinate-identifier circles, rolling-circle amplify them and pattern the resulting DNBs on the chip.
2. Decode CID-to-position correspondence by sequencing.
3. Ligate molecular-identifier/oligo-dT capture sequence to the decoded DNB oligo.
4. Capture tissue RNA and perform on-chip reverse transcription with template switching.
5. Release and amplify CID–MID-bearing cDNA, construct a DNBSEQ library, circularize and form sequencing DNBs.

## 4. Evidence boundary

🟢 The published Cell table prints the TSO, cDNA-PCR primer, library F/R and splint, and read primers. The in-house Tn5 transposome loading adapter and exact terminal library junction are not specified in the inspected article. The actual historical MGI sequencer manual and separately cited Table S1 still need source review; absence is not accepted. Published Supplemental Figures S1–S7 are recovered in the author PDF (see below). The physical final library remains unconstructed, with an explicit validation warning.

🟢 There is a published-source discrepancy in the CID precursor: Cell STAR Methods reagent table (local lines 1365–1369) prints `...TTGAACT...`; bioRxiv v3 Methods precursor line 3900 and the linked MOSTA detailed protocol oligo table print `...TTGNAACT...`, adding one degenerate N before the fixed handle. The model uses `...TTGNAACT...` because both author documents show that synthesis design; the precursor actually used in the published experiment remains uncertain. The Cell Methods (local line 1684) call the capture tail “22 nt poly-T,” while its reagent table (lines 1412–1414) prints 21 T bases followed by a V anchor. The model follows the printed oligo.

## 5. Fresh original-source reconciliation

🟢 The recovered Cambridge published PDF has 38 pages, including **Supplemental Figures S1–S7 on pp. 25–38**. Its original figure captions (extracted text 1888 onward) describe chip dimensions, performance, expression/spatial comparisons, cell types and developmental results. These figures were visually surveyed and their full captions read. The former claim that the entire published supplement was inaccessible was incorrect. The separately cited **Table S1** remains an actual outstanding source: the library Methods refer to it for sample-specific read lengths. No run length is filled in from a modern kit default.

🟢 A further tail-length difference is visible in the actual MOSTA HTML: the short phosphorylated capture oligo in its procedure (local text 56) and the published Cell reagent table both contain **21 T bases followed by V**, while the MOSTA full capture-probe table (1276–1279) contains **23 T bases followed by V**. Both procedural descriptions call the tail 22-nt poly-T. The short ligation input and full CID-bearing capture construct are different molecular states; their literal tail lengths also differ. The model uses the published short ligation oligo, retains its 21T–V tail and records the unresolved full-probe difference on the public page. Neither source is silently normalized.

🟢 The selected RNA-library workflow follows the final Cell Methods: SuperScript II and the source-modified TSO; tissue removal; Exonuclease I; KAPA cDNA PCR; in-house Tn5 fragmentation; SDS stop; KAPA library F/R PCR; AMPure cleanup; reported DNB generation. The MOSTA procedure calls its RT reagent “RT Enzyme Mix”; the SuperScript II attribution comes from Cell. Its TSO tail is preserved literally as `/rG//rG//iXNA_G/`; the notation does not establish additional chemical identity for `iXNA_G`.

🟡 The page now removes the generic 18-position primer placeholder, unsupported final-circle construction and assigned RCA enzyme. DNB generation has an explicit validation warning explaining the unconstructed molecular states. Known library/splint/MDA/read oligos are retained separately from that model gap. The initial chip DNB preparation and the final RNA-library DNB generation are distinct stages.

The article actually names an **MGI DNBSEQ-Tx sequencer manual** for initial array RCA/loading and base calling; the historical revision has not been identified. This is distinct from the old hypothetical “adapter/component guide” checklist row. An unassigned Tn5 loading adapter or physical final junction remains a model question. No missing-information claim has passed the two-reader source gate.

## 6. Bounded pass 2026-10-17 (round breadth-first-19)

🟢 Re-read the published Cell PDF (Cambridge `download.pdf.txt`) Key Resources oligos lines 1365–1472 and STAR Methods lines 1673–1755. The chip-generation chain is now drawn from the graph instead of starting at a pre-decoded array: DNB library oligo1 (25-nt CID) and `/5phos/` oligo2 are ligated on splint oligo1 with T4 DNA ligase (37 °C, 2 h; lines 1676–1678), AMPure-purified, PCR-amplified for 12 cycles and purified again (lines 1678–1680). 🟡 Pairing of the splint halves with the oligo1 3′ end and oligo2 5′ end, and of DNB library F (`/5phos/TGTGAGCCAAGGAGT`) and R (`GAACGACATGGCTA`) with the oligo1 5′ end and oligo2 3′ end, is computed by the model; the Methods do not name the PCR polymerase. RCA, circularization and chip loading follow the MGI DNBSEQ-Tx manual (line 1680) and remain a warned boundary. The CID sequencing primer (`CTGCTGACGTACTGAGAGGCATGGCGACCTTATCAG`) ends where the 25-nt CID begins in oligo1, consistent with SE25 decoding (lines 1680–1681).

🟢 Tissue sections are methanol-fixed (−20 °C, 40 min) and permeabilized with 0.1% pepsin in 0.01 M HCl at 37 °C for 12 min before capture (lines 1692–1699); the capture step now carries that condition. The final AMPure purification uses 0.6× and 0.15× bead ratios (line 1724). Remaining open items are unchanged: Table S1 read lengths, the MGI manual revision, the in-house Tn5 loading adapter and the final DNB circle.

