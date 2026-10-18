# TARGET-seq — research notes

## Active native migration, 2026-10-18

Selected full-length96well HSPC route. Native oligo(dT)-ISPCR, sample-specific mRNA primers and SMARTScribe produce actual ss cDNA; the physical TSO ends rGrG+G (RNA/RNA/LNA), while its copied cDNA is normal DNA. Native genomic DNA, endogenous/ERCC RNA and separate target cDNA populations remain explicit. SeqAmp22cycle preamp receives actual mixed originals/copies with ISPCR and target primer pairs; its product is not a final Illumina construct.

Mixed physical allocations feed transcriptome and genotype preparation; they are not molecule-type purification. The detailed protocol uses14µl plus11µlwater and16µlbeads labelled0.6X, whereas STARMethods gives15µl. The graph selects the detailed14µl route and leaves the published difference visible. Two EtOH washes,8µlEB elution/7.5µl transfer precede700pg ATM55°C6min, NT and explicit initial72°C3min extension. Actual partial native strands enter NPM14cycle PCR; no9base gap size is asserted. Source postPCR dilution/beads/twoEtOH/elution/transfer and equimolar pooling are explicit. The KeyResources v2SetA/C versus detailedStep25FC-131-1001 index-kit difference stays public; no canonical8base width or executed transcriptome index length is assigned. The source biological run is75nt single-end.

The independent1.5µl unpurified genotype aliquot undergoes dilution and CS1/CS2-tailed targetPCR1 (20cycles), an actual2.5µl aliquot enters AccessArray PCR2 (10cycles), and unallocatedPCR1 remains unresolved. Only after the source10base single-direction cell barcode is attached are other cell products pooled and0.8X purified. Ordered AccessArray strands/barcode side are unassigned. MiSeq acquired151/10/151 and four source LNA-modified custom sequencing primers are distinct from unmodified CS1/CS2 PCRtags. Target selection and computational QC do not destroy original/RNA/nonproductive products; their actual physical recovery is unresolved.

The exact57nt oligo(dT)-ISPCR is independently source-tested against the original lowercase KeyResources string. The bundle oligo guard cannot match that extraction and is narrowly exempted; no absence conclusion follows. The separate3′-TARGET/K562/cell-line recipes remain outside this selected native graph.

This bounded pass read 3 originals and reran the mandatory bundle primer scraper. Original/source identity was checked in the external manifest. Historical review status and the two-reader missing-information gate are unchanged.

## Historical observations

The earlier observations below retain their original review scope; current source-owned declarations above replace earlier canonical endpoint and recipe defaults.


Defining source: Rodriguez-Meira *et al.*, *Molecular Cell* (2019), [10.1016/j.molcel.2019.01.009](https://doi.org/10.1016/j.molcel.2019.01.009) ([open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6436961/)). Detailed implementation: Rodriguez-Meira *et al.*, *STAR Protocols* (2020), [10.1016/j.xpro.2020.100125](https://doi.org/10.1016/j.xpro.2020.100125).

## 1. What it is

TARGET-seq combines targeted, biallelic mutation detection from genomic DNA and coding DNA with whole-transcriptome sequencing and index-FACS phenotypes from the same plate-sorted cell.

## Sources read

- 🟢 The defining paper's Key Resources, cDNA synthesis, targeted genotyping and whole-transcriptome library sections were read, including the published ordering forms on p. 18 of the open author PDF.
- 🟢 The detailed protocol was used to check primer nesting and branch order, without changing the defining method.
- 🟢 Bounded pass 2026-10-17 (round breadth-first-20): PMC6436961 JATS Key Resources lines 433–500, STAR Methods cDNA Synthesis/Targeted NGS Single-Cell Genotyping/Nextera XT sections lines 593–611; Table S1b (mmc2 workbook, lysis/RT/PCR reagents and programs); Supplemental Experimental Procedures (mmc1) full-length protocol steps 24–46 (lines 7140–7434) and materials list lines 6800–6828. The Fluidigm Access Array guide and the STAR Protocols paper were not read in this pass.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1016/j.molcel.2019.01.009`.

## Key oligos and how they interlock

- 🟢 Full-length transcriptome cDNA uses the paper's oligo(dT)-ISPCR primer, an LNA-ended template-switch oligo 5′-`AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`-3′, and ISPCR 5′-`AAGCAGTGGTATCAACGCAGAGT`-3′. The template-switch oligo's terminal `rG` bases and `+G` LNA are retained in the ordering form.
- 🟢 Targeted pre-amplification primers for gDNA and optionally cDNA are added during RT/PCR. gDNA-specific pairs include intronic sequence; cDNA-specific pairs target coding sequence.
- 🟢 Genotyping PCR1 adds universal CS1 `ACACTGACGACATGGTTCTACA` and CS2 `TACGGTAGCAGAGACTTGGTCT` tails.
- 🟢 PCR2 attaches Illumina-compatible PE1/PE2 arms and a 10-base cell index through CS1/CS2.
- 🟢 Published CS1-seq 5′-`A+CA+CTG+ACGACATGGTTCTACA`-3′ and CS2-seq 5′-`T+AC+GGT+AGCAGAGACTTGGTCT`-3′ custom sequencing primers use those binding regions. `+` denotes an LNA base in the ordered primer; the unmodified A/C/G/T strings in the construct model are their binding sequences. The same Key Resources list also gives CS1rc-seq and CS2rc-seq alternatives.
- 🟢 The paper's separate 3′ TARGET-seq implementation uses barcoded oligo(dT)-ISPCR primers listed in Supplementary Table S2C for 384-well processing. Those barcodes belong to that branch, not to the full-length genotype construct drawn here.

## Construct-changing steps

1. 🟢 Sort one immunophenotyped cell per well into 4.2 µL lysis buffer (Triton X-100, RNase inhibitor, dNTPs, oligo(dT)-ISPCR, Qiagen protease, ERCC). The protease is heat-inactivated at 72 °C for 15 min (Table S1b, TARGET-seq column).
2. 🟢 RT with SMARTScribe, TSO-LNA and mRNA target-specific primers: 42 °C 90 min, then 10 cycles of 50 °C/42 °C 2 min, then 72 °C 15 min (Table S1b).
3. 🟢 SeqAmp PCR with ISPCR plus gDNA and cDNA target-specific primers: 98 °C 3 min, then 20 (cell lines) or 22 (HSPCs) cycles of 98 °C 15 s, 67 °C 20 s, 72 °C 6 min (Table S1b; Methods line 597).
4. 🟢 Split the unpurified product. 15 of 25 µL are diluted with 11 µL water and purified with 16 µL AMPure XP (0.6:1) into 8 µL EB (Methods line 597). A 1.5 µL aliquot is the genotyping input; the detailed protocol dilutes it 1:2 first (mmc1 step 38).
5. 🟢 Transcriptome: Nextera XT at one-quarter volume, 700 pg cDNA with ATM at 55 °C for 6 min, then NT neutralization. NPM PCR is 72 °C 3 min, then 14 cycles with 55 °C annealing. A 0.6:1 AMPure XP cleanup follows (mmc1 steps 25–35).
6. 🟢 Genotype PCR1: KAPA 2G Robust HS ReadyMix with CS1/CS2-tailed target primers (Table S2B), 20 cycles, 60 °C annealing (mmc1 steps 38–40).
7. 🟢 Genotype PCR2: FastStart High Fidelity with 2 µM Fluidigm Access Array single-direction barcode primers, 10 cycles, 60 °C annealing. Barcodes and adaptors add 103 bp (mmc1 steps 41–45).
8. 🟢 Barcoded cells are pooled with a Mosquito HTS and the pool is purified with 0.8:1 AMPure XP (Methods line 601).
9. 🟢 Sequence both branches and match them by well/cell index and FACS record.

## Final libraries and read layout

- 🟢 The targeted genotype library used 151-cycle Read 1, a 10-cycle index read and 151-cycle Read 2 on a MiSeq. Read 1 and Read 2 use CS1-seq and CS2-seq. The index read uses a mix of CS1rc-seq and CS2rc-seq, 500 nM each (Methods line 603).
- 🔴 The inspected sources do not say which PE arm carries the 10-nt Access Array barcode, so the CS1rc/CS2rc primer that reads it is also unstated. They do not print the PE1/PE2 primer sequences either. Before a 2026-10-17 correction, the model drew the cell index on the P5 side. Both arms are now neutral tokens, and the index read is not declared as a primer site.
- 🟢 The full-length transcriptome branch used Nextera XT and 75-base single-end NextSeq 500/550 High Output v2 reads (Methods line 607). The model declares Read 1 = 75. It also lists both index reads as acquired, but their cycle counts are not stated.
- ⚠️ The published sources disagree on the index kit. The Key Resources (lines 443–445) list Nextera XT Index Kit v2 Sets A/C, while supplemental step 25 (mmc1 line 7150) gives the Nextera XT Index Kit FC-131-1001. The page shows both and the model assigns no index set.
- The public final-library panel focuses on the distinctive CS1/CS2 genotype branch; the matched Nextera branch remains visible in the workflow.

## Open questions

- Which arm carries the Access Array barcode, and what the PE1/PE2 primer sequences are. This needs the Fluidigm Access Array Barcode Library original (deferred).
- Which Nextera XT index kit was used for each plate.
- Table S2A/S2B target primers were not re-transcribed in this pass, and the STAR Protocols (2020) detailed method was not re-read.

## Boundary

TARGET-seq+ is a later named optimization incorporating Smart-seq3 features. It should not be silently represented as original TARGET-seq.
