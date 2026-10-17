# Yeast TeloTag nanopore sequencing

## 1. What it is

Sholes TeloTag marks the physical end of each *S. cerevisiae* telomere with poly(A) and a recognizable oligo before Oxford Nanopore whole-genome library construction.

## 2. Sources read

- 🟢 Sholes et al. (2022), Genome Res 32:616, [doi:10.1101/gr.275868.121](https://doi.org/10.1101/gr.275868.121), PMC8997346. Final Methods in `$CHEM_DATA/sources/yeast-telotag__10.1101+gr.275868.121/gr.275868.121_PMC8997346.xml.txt`, "DNA molecular end tagging" lines 185–193 and "Nanopore library preparation and sequencing" lines 201–207, and `gr.275868.121_supp_supp_gr.275868.121_Supplemental_Material.pdf.txt`, Supplemental Table S5, lines 676–709; the PDF page was checked visually. (Earlier passes read identical copies under `sources/doi_10.1101+gr.275868.121/`; the 2026-10-17 fetch placed the article in the protocol's own source directory.)
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparison with later long-read methods.

## 3. Construct and reaction order

1. 🟢 Terminal transferase (NEB M0315, 20 U) tails 2 µg HMW gDNA in 1× terminal transferase buffer, 1× CoCl2 and 5 mM dATP for 1 h at 37 °C; 10 min at 70 °C stops the reaction (line 187). 🟡 TdT can tail any accessible 3′ end, not only telomeres; no physical telomere enrichment follows.
2. 🟢 The TeloTag primer mix (2.5 mM as printed), 1× ThermoPol buffer, 2.5 mM dNTPs, 1 mM ATP and 4 U Sulfolobus DNA Polymerase IV (NEB M0327S) are added to the same reaction: 1 min 56 °C, 10 min 72 °C (line 191). One of five oligo(dT)-TeloTag primers anneals to the tail and the recessed strand is filled.
3. 🟢 400 U T4 ligase (NEB M0202) is added, 20 min at 12 °C, sealing the nick; tagged gDNA is purified with "Ampere XP" (AMPure XP) beads (line 191).
4. 🟢 Short Read Eliminator XS (Circulomics SS-100-121-01) retains molecules >10 kb (line 193).
5. 🟢 1 µg tagged DNA is barcoded and adapted with the ONT native barcoding kit (EXP-NBD104 with SQK-LSK109), run on R9.4.1 MinION or Flongle flow cells (line 203). The kit's end-preparation and barcode-ligation chemistry are not printed; the model draws that as a warned kit boundary.

🟢 Supplemental Table S5 prints the five TeloTag primers (5′→3′, with no terminal modifications stated):

| Member | Oligo as printed |
|---|---|
| 1 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT` |
| 2 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCAC` |
| 3 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACA` |
| 4 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCA` |
| 5 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACC` |

The shared 20-base core is followed by 35 T bases in member 1 and 36 T bases plus a 3-base terminal variant in members 2–5 🟡 (counted from the table). The Methods text points to Supplemental Table S6, but that table lists plasmids; the primer mix appears in S5. The variable length of the terminal-transferase poly(A) addition remains unknown. The model shows a representative 35-base pairing footprint, not a measured tail length 🟡. The EXP-NBD104 barcodes and SQK-LSK109 motor-adapter oligos remain kit-specific 🔴.

Version history (not a published discrepancy): the superseded author bioRxiv v1 Methods (lines 872–874) gave a >5-kb Short Read Eliminator threshold; the final article (line 193) gives >10 kb, which the model uses. Under the published-version policy this ordinary revision is no longer displayed as a source conflict on the page; the S6/S5 table-number error in the final article remains displayed.

## 4. Final library and sequencing

The TeloTag is a terminal landmark in an otherwise long genomic molecule. Nanopore sequencing uses a motor-loaded adapter, not a sequencing primer. 🟢 Only reads with barcodes on both ends pass; telomere-containing reads (TideHunter) and TeloTag-containing reads (SeqKit) are then selected and trimmed (Porechop) computationally (lines 203–207).

## 5. Open questions

- EXP-NBD104 / SQK-LSK109 historical kit protocols (end-prep, barcode ligation, barcode and motor-adapter sequences) not inspected.
- Poly(A) tail length and which TeloTag member pairs at each telomere phase are not measured.
- The printed "2.5 mM TeloTag primer mix" concentration is transcribed as written.
