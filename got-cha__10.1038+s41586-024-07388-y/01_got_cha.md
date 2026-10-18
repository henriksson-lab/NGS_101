# GoT-ChA

**Defining source (published version):** Izzo F, Myers RM *et al.*, “Mapping genotypes to chromatin accessibility profiles in single cells”, *Nature* 629, 1149–1157 (2024), doi:[10.1038/s41586-024-07388-y](https://doi.org/10.1038/s41586-024-07388-y), author manuscript PMC11139586. It supersedes the bioRxiv preprint doi:10.1101/2022.05.11.491515 (Myers, Izzo *et al.*, 2022) (which named this directory until 2026-10-17); preprint values are not used. Earlier versions of this note cited “Namikawa et al.” and called PMC11139586 a detailed protocol; both were repository errors.

## Sources read (2026-10-17 bounded pass)

- 🟢 PMC11139586 JATS (`$CHEM_DATA/sources/got-cha__10.1038+s41586-024-07388-y/s41586-024-07388-y_PMC11139586.xml.txt`) lines 205–217 (“scATAC-seq with GoT-ChA”), 219–229 (“Multiplex-adapted GoT-ChA”), 243–245 (“GoT-ChA primer design”), 347 (Extended Data Fig. 1 legend).
- 🟢 Springer ESM `41586_2024_7388_MOESM3_ESM.xlsx.txt` (Supplementary Table 1, “GoT-ChA primer sequences”), rows 1–34.
- The published files are in this protocol's source bundle (md5-identical copies also remain in the older `doi_10.1038+s41586-024-07388-y` directory).

## Molecular boundary

1. 🟢 Nuclei are prepared per 10x CG000169 Rev D and processed per Chromium Next GEM Single Cell ATAC CG000209 Rev F. At GEM generation (Step 2.1), 1 µL of 22.5 µM GoT-ChA primer mix (GoT-ChA_R1N + GoT-ChA_Rev) is added to the barcoding mix. The target is amplified exponentially while ATAC fragments are amplified linearly.
2. 🟢 At Post GEM Incubation Cleanup (Step 3.2) the beads are eluted in 45.5 µL: 5 µL goes to GoT-ChA, 40 µL to the standard ATAC library.
3. 🟢 Hemi-nested PCR with P5 and biotinylated GoT-ChA_nested: 95 °C 3 min; 15 × (95 °C 20 s, 65 °C 30 s, 72 °C 20 s); 72 °C 5 min. Then 1.2X SPRIselect.
4. 🟢 Dynabeads M-280 Streptavidin, washed 3× in 1X SSPE, bind for 15 min at room temperature. They are washed 2× in SSPE and 1× in 10 mM Tris-HCl pH 8.0, then resuspended in water.
5. 🟢 On-bead indexing PCR with P5 + RPI-X, 6–10 cycles, same program.
6. 🟢 NovaSeq 6000: Read 1N 50, i7 8, i5 16, Read 2N 50 cycles. ATAC is sequenced to 25,000 and GoT-ChA to 5,000 read pairs per nucleus.

## Oligos (verbatim)

| Oligo | Source | 5′→3′ |
|---|---|---|
| GoT-ChA_R1N | Methods 211 | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG - [20–22bp locus specific]` |
| GoT-ChA_Rev | Methods 211 (single target) | `20–22bp locus specific primer` (no handle stated) |
| GoT-ChA_Rev | Methods 223 (multiplex); Supp. Table 1 | `AGCAAGTGAGAAGCATCGTGTC - [20–22bp locus specific]` |
| GoT-ChA_nested | Methods 215 | `/5BiosG/CCTTGGCACCCGAGAATTCCA-[20–22bp locus specific]` |
| P5 | Methods 215 | `AATGATACGGCGACCACCGAGATCTACAC` (= `illumina.P5`, {{= len(illumina.P5) }} nt) |
| p5 | Supp. Table 1 row 9 | `AATGATACGGCGACCACCGAGA` |
| RPI-X | Methods 215; Supp. Table 1 row 33 | `CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` |

🟢 Supplementary Table 1 prints full per-locus sets for JAK2-V617, TP53-M133, TP53-R284, NRAS-Q61 and FOXO1_S22. The IDT order format has three 3′ phosphorothioates on R1 and Rev and `/5BiosG/` on nested. For example, JAK2 R1 is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTGATGGCAGTTGCAGGTCC*A*T*A`. The locus portions are examples, not protocol constants.

## Published discrepancies (shown on the page)

- 🟢 GoT-ChA_Rev handle: the single-target Methods (line 211) describe a plain locus-specific primer. The multiplex section (line 223) says the IS2 handle “was added”. Supplementary Table 1 lists IS2-handled Rev primers for every listed locus. The model leaves the single-target 5′ handle unassigned.
- 🟢 P5: Methods give 29 nt, Supplementary Table 1 “Common Primers” give 22 nt. The model uses the Methods sequence, which is the one tied to the reactions.
- 🟡 Observation, not modelled: Supp. Table 1 NRAS-Q61 GoT-ChA_nested shares its locus portion (`GGGACAAACCAGATAGGCAGAA`) with NRAS GoT-ChA_Rev, unlike the design rule (line 245) that the nested primer lies ≥100 bp downstream of R1N. Its intended geometry was not checked.

## Open questions

- The 10x CG000209 Rev F original was not read in this pass. Its gel-bead primer, barcoding polymerase and the SI-NA-X ATAC indexing structure (Supp. Table 1 row 34) remain role-level.
- The polymerase/master mix for the hemi-nested and indexing PCRs is not named in the inspected passage.
- The multiplex-adapted route (P5/IS2 pre-amplification for 8 cycles, nested PCR reduced to 8 cycles, and a longer annealing time for the first 6 in-droplet cycles) is not drawn. Neither are GoT-ChA-ASAP (ADT spike-in) and the DOGMA-seq extension.

## Local validation migration (2026-10-18)

Re-read published PMC11139586 Methods, `scATAC-seq with GoT-ChA`, text lines
205–217; ran the bundle oligo scraper and checked the P5/RPI-X passages in that
published original. The superseded preprint was not used for active values.

🟢 The cell barcode is the 16-base i5 identity, distinct from the eight-base RPI-X
sample index. The local genotype and parent ATAC plans now declare 50/8/16/50,
with structured read roles. Final endpoint roles and identifier features are owned
here. The parent ATAC branch explicitly delegates library construction to CG000209
Rev F rather than inheriting a generic final PCR with invented reagent stubs.

🟡 The parent kit revision is identified but was not inspected in this bounded
pass. Target-panel and bead junctions remain unassembled, and unavailable PCR
geometry now uses the actual incoming target/captured state. The previously recorded
published P5/GoT-ChA_Rev discrepancies remain visible. No full chemistry or
source-absence review is claimed.
