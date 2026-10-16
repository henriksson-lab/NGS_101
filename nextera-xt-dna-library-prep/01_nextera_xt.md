# Nextera XT DNA Library Prep

## Sources read

🟢 Illumina, *Nextera XT DNA Library Prep Product Documentation*, document 15031942 v07, pp. 1–15. It defines one-step enzymatic fragmentation/partial-adapter transfer, neutralization, limited-cycle index PCR, cleanup and normalization.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, supplies the canonical Nextera mosaic end, S5/S7 entry points and sequencing primers through `lib/nextera.py`.

🟢 Illumina, *Nextera XT Library Prep: Tips and Troubleshooting*, publication 770-2015-015, pp. 1 and 5–6, confirms that transposomes carry partial adapters, limited-cycle PCR completes the sequencing arms, and bead normalization yields single-stranded library. The schematic endpoint is the double-stranded library after PCR cleanup, before optional bead normalization.

## Molecular workflow

1. An engineered solution transposome tagmentates 1 ng double-stranded DNA in 15 minutes.
2. Neutralize Tagment Buffer stops the transposase.
3. Limited-cycle PCR amplifies the insert and adds both index-adapter sequences.
4. Bead cleanup removes short products; libraries are normalized or quantified and pooled.

## Scope and endpoint

🟢 The guide recommends Nextera XT for genomes below 5 Mb, plasmids, double-stranded cDNA and PCR amplicons longer than 300 bp; Illumina DNA Prep is the alternative for large or complex genomes.

🟡 The inspected product guide and technical note describe partial-adapter transfer but do not print a base-by-base inventory of the loaded solution transposome. Its reagent composition remains an open source-review question. The final library has the public Nextera S5/ME–insert–ME/S7 architecture. The model displays the heterologous S5/S7 product because it can acquire both P5 and P7 and amplify.

🟢 The guide offers Nextera XT v2 and other compatible dual-index kits. The index-reference document prints the N7xx and S5xx families, while the pooling guide gives valid combinations. A particular pair is selected for each sample; it is not a fixed sequence of the prep chemistry.
