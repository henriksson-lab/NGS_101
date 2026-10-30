# Subassembly (Hiatt) — native end tags and nested sublibraries

Hiatt et al. (2010), *Parallel, tag-directed assembly of locally derived short sequence reads*, [10.1038/nmeth.1416](https://doi.org/10.1038/nmeth.1416), Nature Methods 7:119–122.

## Sources and version boundary

The authoritative preparation source in this pass is the **published publisher supplementary PDF**, `nmeth.1416_41592_2010_BFnmeth1416_MOESM439_ESM.pdf`, manifest MD5 `26bb14af6e4c525034e9b54afa9cb32a`. Inspected: Supplementary Protocols 1–3, the metagenomic route in detailed Protocol 2, and Table 4. The rendered original Table 4, printed p. 11 / PDF p. 12, was inspected to resolve row assignment, line wraps, 5′ phosphates and the terminal phosphorothioate bond. The standard source fetcher and full-bundle oligo scraper were run; Table 4 establishes the sequences and modifications.

PMC2848820 is an **accepted author manuscript**, not the final published article. Its Introduction/native-tag explanation and Methods overview were read as version-qualified orientation. The published SI supplies the actual modeled operations and oligos; no complete equivalence of the two article versions is claimed. The final article's full Methods remain uninspected. The PMC manuscript supplement download is gated, while the publisher SI is available; their complete equivalence has not been established from a filename.

This is a capped two-original pass, not a full source/linked-material or missing-information review. The private graph follows the **published SI metagenomic route**, with its separate merging-library output. The P. aeruginosa branch has different initial/final size ranges and omits metagenomic-specific resizing; it is not silently mixed into this instance.

## What identifies a molecule

🟢 Subassembly uses short reads from the native ends of source fragments as tags. The author-manuscript explanation makes this explicit; the published SI's fixed tag-adjacent adapters and 20-base tag read are consistent with that native-end strategy. The model does **not** turn this into an added random 20-base UMI. Tag collisions and repeated genomic sequences remain possible. Genomic end sequences are unassigned; the displayed N20 component stands for the source tag-read span, not a synthesized oligo.

🟢 Complexity is restricted before PCR. Many copies of the retained fragment population are blunt-ligated into concatemers and sheared. An original fragment end remains accessible through its tag-adjacent adapter while a second read samples different shearing breakpoints. Computational grouping and local assembly generate a synthetic long sequence; they do not produce a physically joined long molecule.

## Published metagenomic operations

1. 🟢 Fragment metagenomic DNA with Bioruptor, purify by QIAquick, End-It end-repair and purify again. Select 400–500 bp by 6% TBE electrophoresis and ethanol recovery. Add terminal 3′ A with the source Taq/dATP reaction at 70°C for 20 min and purify.
2. 🟢 Anneal Ad1/Ad1_rc and Ad2/Ad2_rc separately; ligate the source T-tailed adapters with Quick Ligase. Select 400–800 bp to remove free adapters. A Scene reconstructs each source-written duplex, its terminal 3′ T and phosphorylated bottom strand. Full insert/adaptor product geometry remains unassembled.
3. 🟢 Serial dilution followed by real-time PCR imposes a complexity bottleneck. Both amplification primers are physically 5′ phosphorylated. Stop amplification before the end of log phase. Purify; the metagenomic route additionally selects 500–600 bp and repeats amplification. Step 8c divides material among eight reactions, then pools and purifies products. Neither the dilution choice nor the amplification kinetics is replaced with a universal fixed cycle count.
4. 🟢 A separate metagenomic aliquot enters **Step 8d merging-library PCR** with Illum_amp_f_Ad1 and Illum_amp_r_Ad2, using iProof HF. Published SI specifies independent 36×36 paired-end sequencing with Ad1_seq and Ad2_seq. This library pairs native tags from opposite ends computationally. It is not the product of the main breakpoint-adapter workflow and is not physically pooled into that workflow by assumption. The author-manuscript Results reports a different trimmed/read description; the model uses published SI's specified acquisition without treating version differences as a publication error.
5. 🟢 The main aliquot undergoes blunt intermolecular ligation of phosphorylated PCR products into high-molecular-weight concatemers. The **existing molecular population supplies both partners**; no new Ad1 amplification primer is pipetted into this ligation. Bioruptor re-fragmentation, QIAquick cleanup, End-It repair/cleanup, Taq A-tailing/cleanup follow.
6. 🟢 Anneal and ligate Illum_rev/Illum_rev_rc as the breakpoint-adjacent adapter; purify. Its source top strand carries a phosphorothioate bond immediately before the terminal T and its bottom strand a 5′ phosphate. These modifications belong to the physical input oligos, not every polymerase-made copy.
7. 🟢 Physically separate aliquots enter Ad1-specific and Ad2-specific PCRs, each paired with Illum_amp_r. Purify, select the metagenomic 300–450 bp range, perform the final material-generation PCR and QIAquick cleanup. Published SI specifies **Read 1 = 20-base tag, Read 2 = 76-base breakpoint**, with Ad1_seq or Ad2_seq and Illum_seq_r. Complete final product/run-primer geometry remains unassembled.

The model preserves physical fractions and untransferred aliquots. Soluble reaction products persist until physical purification. It does not claim perfect ligation, clonal purity, unbiased size recovery, collision-free tags or a completely assembled final library.

## Source oligos

| Source input | Verbatim or model-derived source ordering string |
|---|---|
| Ad1 | `{{= chem.AD1 }}` |
| Ad1_rc | `/5Phos/{{= chem.revcomp(chem.AD1[:-1]) }}` |
| Ad2 | `{{= chem.AD2 }}` |
| Ad2_rc | `/5Phos/{{= chem.revcomp(chem.AD2[:-1]) }}` |
| Ad1_amp | `{{= chem.OLIGOS['ad1amp'].ordering() }}` |
| Ad2_amp | `{{= chem.OLIGOS['ad2amp'].ordering() }}` |
| Illum_amp_f_Ad1 | `{{= chem.F1 }}` |
| Illum_amp_f_Ad2 | `{{= chem.F2 }}` |
| Illum_amp_r | `{{= chem.RBP }}` |
| Merging reverse | `{{= chem.MERGE_R }}` |
| Illum_rev | `{{= chem.BP[:-1] }}*{{= chem.BP[-1] }}` |
| Illum_rev_rc | `/5Phos/{{= chem.revcomp(chem.BP[:-1]) }}` |
| Ad1_seq / Ad2_seq | `{{= chem.SEQ1 }}` / `{{= chem.SEQ2 }}` |
| Illum_seq_r | `{{= chem.SEQBP }}` |

The source `*` denotes a phosphorothioate **bond** and `/5Phos/` a non-nucleotide terminal phosphate. Canonical matching P5/P7 regions are imported; the other fixed flanks are not assigned undocumented sample-index functions.

## Lineage and remaining scope

BAsE-Seq's Introduction cites and describes Hiatt's earlier preparation, but that paragraph alone does not assert direct adaptation of a molecular component. Similar tag-directed assembly goals do not establish shared reaction chemistry or justify importing BAsE/Loop constructors.

Deferred: final published main article/Methods, the gated manuscript SI's exact identity relative to the publisher SI, uninspected SI figures/results, complete native-DNA/concatemer/breakpoint/read geometry, exact physically retained primer modifications through copying, historical source kit documentation where relevant, remaining linked-source review and independent review. No missing-information claim is accepted; the preview remains private.
