# TruSeq Stranded mRNA

## Sources read

🟢 Illumina, *TruSeq Stranded mRNA Reference Guide*, document 1000000040498 v00, pp. 1–5 and 9–23. It defines the two-round poly(A) capture, heat/divalent-cation fragmentation, random-prime first strand, dUTP-marked second strand, A-tailing, adapter ligation and PCR enrichment.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, supplies the canonical TruSeq adapter and sequencing-primer sequences through `lib/illumina.py`.

## Molecular workflow

1. Poly-T magnetic beads select polyadenylated RNA twice. Divalent cations at elevated temperature fragment the selected RNA.
2. Random primers and reverse transcriptase synthesize first-strand cDNA. Actinomycin D suppresses DNA-dependent synthesis.
3. RNase H and DNA polymerase I make the complementary strand while dUTP replaces dTTP.
4. The blunt cDNA is dA-tailed and receives T-overhang TruSeq adapters.
5. PCR enriches adapter-bearing molecules. The enrichment polymerase does not copy past dUTP, so the marked second strand is quenched.

## Direction and indexes

🟢 The guide states that Read 1 maps to the antisense strand and Read 2 maps to the sense strand. That orientation is encoded in the insert label rather than left implicit.

🟢 The guide supports single-index LS and dual-index HS libraries. The schematic shows the dual-index form. 🟡 Specific index bases depend on the purchased adapter plate and remain typed sample-index placeholders.
