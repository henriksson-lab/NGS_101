# spatial-CUT&Tag–RNA-seq

**Defining paper:** Zhang et al., *Nature* (2023), doi:10.1038/s41586-023-05795-1 (spatial CUT&Tag–RNA-seq branch). Kept separate from spatial-ATAC–RNA-seq because antibody-targeted pA–Tn5 replaces untargeted accessibility tagmentation; RT, ligations, separation and the cDNA library are shared (Methods line 288).

## Sources read (bounded pass 2026-10-17, round breadth-first-19)

- 🟢 Defining article, published version: Zhang D. et al., *Nature* 616, 113–122 (2023), PMC10076218 JATS (`s41586-023-05795-1_PMC10076218.xml.txt`) Methods "Preparation of transposome" lines 245–252, "DNA barcode sequences…" 254–256, "Spatial ATAC–RNA-seq" 262–284, "Spatial CUT&Tag–RNA-seq" 286–292, "Data preprocessing" 294–298; Results line 136.
- 🟢 Its Supplementary Information (`…MOESM1_ESM.pdf.txt`) Table 4 lines 594–657 (RT primer, linkers 1/2, PCR Primers 1/2, N501, N701–N707), Table 5 lines 658–1066 (/5Phos/ barcode A-1…A-100), Table 6 lines 1067–1475 (barcode B-1…B-100).
- 🟡 Later detailed protocol (Nature Protocols 2025, doi:10.1038/s41596-025-01145-9): only its Supplementary Table 1 workbook was available; its oligos are letter-identical to the 2023 Table 4 entries used here, but the 2023 article and supplement are the authority. The 2025 main text was not available.
- Oligo scraper run on the source bundle; every sequence used was checked in the 2023 passages above.

## Step by step (2023 Methods lines 286–290, shared steps 262–284)

1. 🟢 Thaw, fix (0.2% formaldehyde 5 min, glycine quench), wash, permeabilize with NP40-digitonin wash buffer 5 min.
2. 🟢 Primary antibody against the histone mark (1:50, 4 °C overnight); guinea pig anti-rabbit secondary antibody (1:50, 30 min).
3. 🟢 pA–Tn5 adaptor complex (Diagenode C01070002 loaded with the same Tn5ME-A/B oligos) 1:100 in 300-wash buffer, 1 h; tagmentation with 10 mM MgCl2 at 37 °C 1 h; stop with 40 mM EDTA.
4. 🟢 RT, barcode A and B ligations, reverse crosslinking and streptavidin separation as in spatial ATAC–RNA-seq.
5. 🟢 CUT&Tag library from the supernatant: NEBNext with P5 PCR primer and indexed i7 primer; 58 °C 5 min, 72 °C 5 min, then 12 cycles; 1.3× AMPure.
6. 🟢 cDNA library as in spatial ATAC–RNA-seq (on-bead TSO, PCR Primers 1/2, Nextera XT, 12 cycles).

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


## Native graph migration (2026-10-18, batch10)

The source-owned graph now follows the original 2023 assay through tissue processing,
transposition, RT, both spatial ligations, reverse crosslinking and streptavidin
partition. DNA and RNA are independent physical populations. The source RT primer
contains a physical internal biotin-dT followed by T15VN; its ordinary DNA copies
retain that nucleotide without inheriting biotin. Modified TSO originals and their
DNA copies likewise have distinct chemistry.

The DNA PCR receives the actual native gapped material after the published initial
extension. N501/N701 and compatible standard Nextera sequencing primers form an
**inferred reference**, with sample assignment and first-cycle completion unresolved.
The RNA branch retains on-bead template switching, five PCR cycles, physical bead
removal, free-solution amplification, 0.8× recovery, a 1ng aliquot, XT tagmentation,
neutralization, twelve PCR cycles and 0.7× recovery. Possible free RNA, fragments,
primers and non-target DNA remain explicit; bead removal alone does not prove their
complete removal. Native final RNA primer ends and index acquisition stay unassigned.
Both preparations retain the source paired-end150 acquisition; physical identifier
widths do not supply index-read cycle counts.

This pass read the published Methods244–298 and original 2023 SupplementaryInformation
Table4 (594–657), first A1/B1 rows and their barcode-family context in Tables5–6
(658–681,1067–1088). It inspected two originals; it did not reread the later 2025
protocol supplement or rendered figures. The earlier source reviews remain history.
Actual study primer/index choices, native gaps and the detailed XT oligos remain
follow-up questions, without an accepted publication-absence claim.
