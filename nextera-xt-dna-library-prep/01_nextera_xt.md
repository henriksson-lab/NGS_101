# Nextera XT DNA Library Prep

## Sources read

🟢 Illumina, *Nextera XT DNA Library Prep Product Documentation*, document 15031942 v07, pp. 1–15. It defines one-step enzymatic fragmentation/partial-adapter transfer, neutralization, limited-cycle index PCR, cleanup and normalization.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, supplies the canonical Nextera mosaic end, S5/S7 entry points and sequencing primers through `lib/nextera.py`.

## Molecular workflow

1. An engineered solution transposome tagmentates 1 ng double-stranded DNA in 15 minutes.
2. Neutralize Tagment Buffer stops the transposase.
3. Limited-cycle PCR amplifies the insert and adds both index-adapter sequences.
4. Bead cleanup removes short products; libraries are normalized or quantified and pooled.

## Scope and endpoint

🟢 The guide recommends Nextera XT for genomes below 5 Mb, plasmids, double-stranded cDNA and PCR amplicons longer than 300 bp; Illumina DNA Prep is the alternative for large or complex genomes.

🟡 The commercial transposome mixture is not disclosed base by base. The final library nevertheless has the public Nextera S5/ME–insert–ME/S7 architecture. The model displays the heterologous S5/S7 product because it is the product that can acquire both P5 and P7 and amplify.
