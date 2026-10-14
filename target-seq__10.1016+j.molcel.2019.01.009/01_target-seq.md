# TARGET-seq — research notes

Defining source: Rodriguez-Meira *et al.*, *Molecular Cell* (2019), [10.1016/j.molcel.2019.01.009](https://doi.org/10.1016/j.molcel.2019.01.009) ([open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6436961/)). Detailed implementation: Rodriguez-Meira *et al.*, *STAR Protocols* (2020), [10.1016/j.xpro.2020.100125](https://doi.org/10.1016/j.xpro.2020.100125).

## 1. What it is

TARGET-seq combines targeted, biallelic mutation detection from genomic DNA and coding DNA with whole-transcriptome sequencing and index-FACS phenotypes from the same plate-sorted cell.

## Sources read

- 🟢 The defining paper's Key Resources, cDNA synthesis, targeted genotyping and whole-transcriptome library sections were read, including the published ordering forms on p. 18 of the open author PDF.
- 🟢 The detailed protocol was used to check primer nesting and branch order, without changing the defining method.
- The source fetch command is `python3 tools/get_sources.py --doi 10.1016/j.molcel.2019.01.009`.

## Key oligos and how they interlock

- 🟢 Full-length transcriptome cDNA uses the paper's oligo(dT)-ISPCR primer, an LNA-ended template-switch oligo 5′-`AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`-3′, and ISPCR 5′-`AAGCAGTGGTATCAACGCAGAGT`-3′. The template-switch oligo's terminal `rG` bases and `+G` LNA are retained in the ordering form.
- 🟢 Targeted pre-amplification primers for gDNA and optionally cDNA are added during RT/PCR. gDNA-specific pairs include intronic sequence; cDNA-specific pairs target coding sequence.
- 🟢 Genotyping PCR1 adds universal CS1 `ACACTGACGACATGGTTCTACA` and CS2 `TACGGTAGCAGAGACTTGGTCT` tails.
- 🟢 PCR2 attaches Illumina-compatible PE1/PE2 arms and a 10-base cell index through CS1/CS2.
- 🟢 Published CS1-seq 5′-`A+CA+CTG+ACGACATGGTTCTACA`-3′ and CS2-seq 5′-`T+AC+GGT+AGCAGAGACTTGGTCT`-3′ custom sequencing primers use those binding regions. `+` denotes an LNA base in the ordered primer; the unmodified A/C/G/T strings in the construct model are their binding sequences. The same Key Resources list also gives CS1rc-seq and CS2rc-seq alternatives.
- 🟢 The paper's separate 3′ TARGET-seq implementation uses barcoded oligo(dT)-ISPCR primers listed in Supplementary Table S2C for 384-well processing. Those barcodes belong to that branch, not to the full-length genotype construct drawn here.

## Construct-changing steps

1. 🟢 Sort one immunophenotyped cell per well into mild protease-containing lysis buffer.
2. 🟢 Reverse-transcribe RNA and co-amplify whole-transcriptome cDNA, target cDNA and target gDNA in one reaction.
3. 🟢 Split the common product. Tagment purified Smart-seq+ cDNA with Nextera XT for expression; perform two nested PCRs on the target aliquot for genotyping.
4. 🟢 Sequence and match both branches by the well/cell index and FACS record.

## Final libraries and read layout

- 🟢 The targeted genotype library used 151-cycle Read 1, a 10-cycle index read and 151-cycle Read 2 with CS1/CS2 custom primers.
- 🟢 The full-length transcriptome branch used Nextera XT and 75-base single-end sequencing in the defining study.
- The public final-library panel focuses on the distinctive CS1/CS2 genotype branch; the matched Nextera branch remains visible in the workflow.

## Boundary

TARGET-seq+ is a later named optimization incorporating Smart-seq3 features. It should not be silently represented as original TARGET-seq.
