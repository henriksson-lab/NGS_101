# spatial-ATAC–RNA-seq

**Defining paper:** Zhang et al., *Nature* (2023), doi:10.1038/s41586-023-05795-1 (spatial ATAC–RNA-seq branch).

## Sources read (bounded pass 2026-10-17, round breadth-first-19)

- 🟢 Defining article, published version: Zhang D. et al., *Nature* 616, 113–122 (2023), PMC10076218 JATS (`s41586-023-05795-1_PMC10076218.xml.txt`) Methods "Preparation of transposome" lines 245–252, "DNA barcode sequences…" 254–256, "Spatial ATAC–RNA-seq" 262–284, "Spatial CUT&Tag–RNA-seq" 286–292, "Data preprocessing" 294–298; Results line 136.
- 🟢 Its Supplementary Information (`…MOESM1_ESM.pdf.txt`) Table 4 lines 594–657 (RT primer, linkers 1/2, PCR Primers 1/2, N501, N701–N707), Table 5 lines 658–1066 (/5Phos/ barcode A-1…A-100), Table 6 lines 1067–1475 (barcode B-1…B-100).
- 🟡 Later detailed protocol (Nature Protocols 2025, doi:10.1038/s41596-025-01145-9): only its Supplementary Table 1 workbook was available; its oligos are letter-identical to the 2023 Table 4 entries used here, but the 2023 article and supplement are the authority. The 2025 main text was not available.
- Oligo scraper run on the source bundle; every sequence used was checked in the 2023 passages above.

## Step by step (2023 Methods lines 262–284)

1. 🟢 Thaw, fix (0.2% formaldehyde 5 min, glycine quench), permeabilize (lysis buffer with digitonin/NP40, 15 min).
2. 🟢 Tagment with home-made Tn5 transposome (Diagenode unloaded Tn5 C01070010 loaded with Tn5ME-A and /5Phos/Tn5ME-B, each with /5Phos/Tn5MErev), 37 °C 30 min; stop with 40 mM EDTA.
3. 🟢 In situ RT: Maxima H Minus, biotinylated RT primer, room temperature 30 min then 42 °C 90 min.
4. 🟢 Barcode A (pre-annealed with linker 1) ligated with T4 DNA ligase in the first channel set, 37 °C 30 min; then barcode B (pre-annealed with linker 2) in the perpendicular set.
5. 🟢 Reverse crosslinking: proteinase K/SDS 58 °C 2 h, 65 °C overnight; Zymo DCC-5.
6. 🟢 Streptavidin C1 beads capture biotinylated cDNA; the supernatant is the ATAC branch.
7. 🟢 ATAC library: NEBNext with indexed i7 and P5 primer; 72 °C 5 min, then 5 cycles plus qPCR-determined cycles; 1× AMPure.
8. 🟢 cDNA: on-bead template switching (Maxima H Minus + TSO), KAPA HiFi pre-amplification with PCR Primers 1 and 2 (5 cycles + qPCR), 0.8× AMPure; Nextera XT tagmentation of 1 ng, 12-cycle PCR with P5 and indexed P7, 0.7× AMPure.

## Oligos (verbatim, 2023 Methods and Supplementary Table 4–6)

- Tn5ME-A `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`; Tn5MErev `/5Phos/CTGTCTCTTATACACATCT`; Tn5ME-B `/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG`.
- RT primer `/5Phos/CATCGGCGTACGACTNNNNNNNNNN/iBiodT/TTTTTTTTTTTTTTTVN`.
- Barcode A-1 `/5Phos/AGGCCAGAGCATTCGAACGTGATGTGGCCGATGTTTCG` (A-1…A-100); barcode B-1 `CAAGCGTTGGCTTCTCGCATCTAACGTGATATCCACGTGCTTGAG` (B-1…B-100, no 5′ phosphate printed).
- Ligation linker 1 `AGTCGTACGCCGATGCGAAACATCGGCCAC`; linker 2 `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`.
- Template switch primer `AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G`; PCR Primer 1 `CAAGCGTTGGCTTCTCGCATCT`; PCR Primer 2 `AAGCAGTGGTATCAACGCAGAGT`.
- N501 and N701–N707 as printed in Table 4; the model draws N501/N701.

🟡 Derived: Tn5ME-B and the RT primer share the 15-nt handle `CATCGGCGTACGACT`; linker 1 is its reverse complement followed by the reverse complement of the barcode-A ligating arm, and linker 2 bridges the barcode-A outgoing arm to the barcode-B ligating arm. The schematic derives both splinted nicks from these sequences. N501's 3′ arm lands on Tn5ME-A and N701's 3′ arm on the barcode-B PCR arm (= PCR Primer 1).

## Read layout

🟢 NovaSeq 6000, paired-end 150. For DNA reads, linkers 1/2 filter Read 2, which carries barcodes A and B; Read 1 is genomic. For RNA, Read 2 is parsed for barcode A, barcode B and the UMI; Read 1 is mapped. Index-read lengths are not printed.

## Open questions

- 🔴 Which "P5 PCR primer" (and "new P5 PCR primer" in the qPCR) the DNA PCR used is named but not sequence-identified beyond the Table 4 N501 entry; the model uses N501.
- 🔴 Which P5 and indexed P7 primers the RNA Nextera XT 12-cycle PCR used is not printed, so the RNA final library stays role-level and its tagmentation off-product fate unresolved.
- 🟡 Candidate A–A / non-B-tagged transposition products stay with the DNA branch (no physical separation); fate not observed in the source.
- Deferred: 2025 Nature Protocols main text; DBiT-seq predecessor (doi:10.1016/j.cell.2020.10.026); Diagenode Tn5/pA-Tn5 assembly guidelines; Nextera XT kit primers.
