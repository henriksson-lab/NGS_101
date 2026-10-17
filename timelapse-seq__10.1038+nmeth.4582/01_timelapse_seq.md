# TimeLapse-seq — K562 total-RNA route

[Schofield et al., Nature Methods (2018)](https://doi.org/10.1038/nmeth.4582) uses chemical recoding to add time information to RNA-seq. Cellular RNA synthesized during a pulse can incorporate 4-thiouridine. Periodate/TFEA chemistry converts labelled sites into a cytidine analogue **C***, not ordinary RNA cytidine; its altered pairing produces apparent U-to-C transitions in sequencing reads.

## Represented source branch

The graph represents the original **K562 four-hour total-RNA experiment**, without affinity enrichment. It distinguishes RNA classes by pulse incorporation. A transcript containing a pulse-labelled tract may also contain a pre-existing 5′ tract; the classification does not date every base of the molecule. A representative converted site does not mean every labelled molecule or site converted; possible unconverted pulse-labelled RNA remains explicit.

The original Supplementary Protocol gives this order: modified RNeasy extraction; NaIO4/TFEA recoding; RNAClean purification; DTT reducing treatment; a second RNAClean purification; named SMARTer Pico kit entry. The graph now preserves those stages. The main article's K562 mRNA paragraph uses 10 ng of total RNA per library and paired-end 150 bp HiSeq 4000 reads; the Supplementary Protocol recoding example (step 13) starts from 2 µg RNA.

The targeted Actb/Gapdh/MYC assays, transcriptome-wide MEF experiments and enriched TT-TimeLapse experiments are distinct source branches. Targeted primers and their TruSeq adapter completion are not imported into this total-RNA graph.

## Published pulse discrepancy

The original Supplementary Protocol, MOESM3 step 2, prints **100 mM** final 4sU for a four-hour K562 pulse. The original main article's “Samples for TimeLapse-seq analysis of K562 mRNA” paragraph explicitly reports **100 µM**. The original Supplementary Note, “Important parameters,” agrees with 100 µM and gives a successful range of 10 µM–1 mM.

The page follows the reported 100 µM experiment for its source description and leaves the graph's pulse concentration unspecified. The intended value of the Protocol's 100 mM entry remains uncertain. This discrepancy is also shown on the public schematic.

## Targeted-assay oligos retained as a separate branch

🟢 Original Supplementary Table S1b (MOESM4, Table S1, rows headed “Targeted TimeLapse-seq primer sequences”) prints twelve targeted primer pairs: four mouse Actb, four mouse Gapdh and four human MYC pairs. Each forward oligo has the common prefix `CTACACGACGCTCTTCCGATCTNNNNN` before its gene-specific arm; each reverse oligo has `CAGACGTGTGCTCTTCCGATCT`. The five N positions are retained as source placeholders; no identifier meaning is inferred here.

The complete source mActB1 example, with spaces separating common and gene-specific parts, is:

- Forward: `CTACACGACGCTCTTCCGATCTNNNNN ACTGAGCTGCGTTTTACACCC`.
- Reverse: `CAGACGTGTGCTCTTCCGATCT TCCTGAGTCAAAAGCGCCAAAAC`.

The same table's “PCR Amplification primers” row prints `CAGACGTGTGCTCTTCCGATC` and `CTACACGACGCTCTTCCGATCT`. The shorter amplification primer must not be silently replaced by the longer common reverse-oligo prefix.

These positive source facts describe the targeted RT-PCR/TruSeq-index branch. They are retained for comparison and do not enter the K562 total-RNA graph or define the SMARTer kit's adapters.

## Conditional kit implementation

The paper names SMARTer Stranded Total RNA-Seq Pico Input with ribosomal cDNA depletion. The [original-generation manual 112216](https://manuals.plus/m/fa0d970016b83cb070dba8c79d0f92b6efcf3bab7ff2ce4ba5361703a491b5dd.pdf), cat. 635005–635007, provides a dated, source-positive implementation:

1. Input-dependent fragmentation or a no-fragmentation option.
2. N6/random priming and SMARTScribe RT with LNA template switching.
3. PCR1 adds full indexed Illumina adapters.
4. **Two consecutive AMPure cleanups** before ZapR/R-Probe ribosomal cDNA cleavage.
5. PCR2 amplifies remaining intact library templates, followed by final cleanup.

These kit stages are conditional in the graph because the study's exact revision and input-dependent choice are unassigned. The unresolved input preparation is a warning-producing step with its reason attached. Exact RT/TSO/primer oligos, selected indexes and index widths are not assigned from a generic shared recipe. The manual's proprietary-mix labels are source statements, not a publication-wide absence claim.

The cleaved ribosomal-cDNA products are shown as an inferred logical class with unresolved physical fate. They are not depicted as a separate fraction removed before PCR2. The graph does not establish complete depletion, quantitative amplification or conversion yields.

## Sequencing and other source obligations

The represented K562 experiment used paired-end 150-base sequencing on HiSeq 4000. The MEF transcriptome-wide branch used paired-end 100, while targeted MEF libraries used paired-end 75 on HiSeq 2500. Exact physically loaded run primers remain unassigned; a PCR or kit-name match does not establish those oligos.

The TT-TimeLapse K562 branch additionally follows [Duffy 2015 MTS-biotin enrichment](https://doi.org/10.1016/j.molcel.2015.07.023), with TimeLapse-specific wash/purification changes. It is outside this non-enriched graph. The original Duffy article, figure supplement and three previously unretrieved workbooks have now been recovered. Their tables contain alignment statistics, miRNA turnover results, synthetic RNA spikes and qPCR primers. They do not contain the separately cited Detailed Protocol; that named original remains to be located. This finding is restricted to those inspected files.

The source audit remains in progress. Unassigned experimental kit choices are model questions, not invented unavailable specification documents. The remaining named original and linked-source work stays explicit; no missing-information claim is accepted.
