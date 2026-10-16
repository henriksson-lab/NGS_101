# Yeast TeloTag nanopore sequencing

## 1. What it is

Sholes TeloTag marks the physical end of each *S. cerevisiae* telomere with poly(A) and a recognizable oligo before Oxford Nanopore whole-genome library construction.

## 2. Sources read

- 🟢 Sholes et al. (2022), [doi:10.1101/gr.275868.121](https://doi.org/10.1101/gr.275868.121), full text and `_data/sources/doi_10.1101+gr.275868.121/gr.275868.121_supp_supp_gr.275868.121_Supplemental_Material.pdf.txt`, Supplemental Table S5, lines 679–707; the PDF page was checked visually.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), comparison with later long-read methods.

## 3. Construct and reaction order

1. 🟢 Terminal transferase adds poly(A) to native telomeric 3′ ends.
2. 🟢 One of five oligo(dT)-TeloTag primers anneals to that tail.
3. 🟢 Sulfolobus DNA polymerase IV fills the recessed strand without strand displacement; T4 ligase seals the nick.
4. 🟢 Oxford Nanopore adapters are added with the manufacturer workflow.

🟢 Supplemental Table S5 prints the five TeloTag primers (5′→3′, with no terminal modifications stated):

| Member | Oligo as printed |
|---|---|
| 1 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT` |
| 2 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCAC` |
| 3 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACA` |
| 4 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCA` |
| 5 | `CAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTACC` |

The shared 20-base core is followed by 35 T bases in member 1 and 36 T bases plus a 3-base terminal variant in members 2–5 🟡 (counted from the table). The Methods text points to Supplemental Table S6, but that table lists plasmids; the primer mix appears in S5. The variable length of the terminal-transferase poly(A) addition remains unknown. The model shows a representative 35-base pairing footprint, not a measured tail length 🟡. The EXP-NBD104 barcodes and SQK-LSK109 motor-adapter oligos remain kit-specific 🔴.

🟢 The author bioRxiv v1 Methods (`2021.06.07.447263_v1.full.pdf.txt`, lines 872–874) describe Short Read Eliminator XS as retaining molecules **above 5 kb**. The final Genome Research Methods (`gr.275868.121_PMC8997346.xml.txt`, line 193) say the same kit retains molecules **above 10 kb**. The page follows the final publication for the size-selection threshold; why the wording changed between versions is unknown. This source-version discrepancy is separate from the S6/S5 primer-table error.

## 4. Final library and sequencing

The TeloTag is a terminal landmark in an otherwise long genomic molecule. Nanopore sequencing uses a motor-loaded adapter, not a sequencing primer.
