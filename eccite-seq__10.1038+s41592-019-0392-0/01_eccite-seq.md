# ECCITE-seq — research notes

## 1. What it is

ECCITE-seq extends a 10x 5′ workflow to recover transcriptome, clonotypes, surface proteins, sample hashtags and directly captured sgRNAs with a shared single-cell barcode.

## 2. Sources read

- 🟢 Mimitou et al., *Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells*, Nature Methods 16, 409–412 (2019), DOI [10.1038/s41592-019-0392-0](https://doi.org/10.1038/s41592-019-0392-0), open manuscript [PMC6557128](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557128/).
- 🟢 Methods and Figure 1 establish the bead-TSO-complementary antibody tags, separate protein/hashtag handles, scaffold RT primer, guide copying and template switching.

## 3. Construct-changing steps

1. Label cells with antibody-derived tags and sample hashtags carrying distinct amplification handles.
2. Add a scaffold-complementary RT primer to the 10x 5′ RT mix.
3. Copy sgRNA through its variable protospacer and template-switch onto the bead-derived handle.
4. Transfer cell barcode and UMI to transcripts and all short-tag modalities.
5. Enrich gene-expression/V(D)J, ADT, hashtag and guide-derived oligo libraries separately.

## 4. Evidence boundary

🟢 The paper's supplement supplies oligo sequences and PCR recipes. 🟡 The page emphasizes library-changing branches and represents experiment-specific barcodes as placeholders. It does not collapse ECCITE-seq into the later 5′ direct-capture Perturb-seq protocol already documented separately.
