# ECCITE-seq — research notes

## 1. What it is

ECCITE-seq extends a 10x 5′ workflow to recover transcriptome, clonotypes, surface proteins, sample hashtags and directly captured sgRNAs with a shared single-cell barcode.

## 2. Sources read

- 🟢 Mimitou et al., *Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells*, Nature Methods 16, 409–412 (2019), DOI [10.1038/s41592-019-0392-0](https://doi.org/10.1038/s41592-019-0392-0), open manuscript [PMC6557128](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557128/).
- 🟢 Methods and Figure 1 establish the bead-TSO-complementary antibody tags, separate protein/hashtag handles, scaffold RT primer, guide copying and template switching.
- 🟢 [Author supplementary protocol, pp. 5–6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0392-0/MediaObjects/41592_2019_392_MOESM3_ESM.pdf), **Oligonucleotide sequences** table, prints tag examples, amplification primers, guide RT primer and their terminal phosphorothioates (`*`). Supplementary Tables 3–5 separately list the experiment's antibody and hashtag barcode assignments.

## 3. Ordered oligos

🟢 The author protocol prints these 5′→3′. Lowercase `x` is the chosen sample index, and `*` marks a phosphorothioate bond; both are retained as written.

| Role | Sequence |
|---|---|
| Protein-tag example | `/5AmMC12/CCTTGGCACCCGAGAATTCCAxxxxxxxxCCCATATAAGAAA` |
| Hashtag example | `/5AmMC12/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTxxxxxxxxCCCATATAAGAAA` |
| Guide RT primer, sg_RT_v4 | `AGCAAGTGAGAAGCATCGTGTCAAAGCACCGACTCGGTGCCAC` |
| 10x SI-PCR | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C` |
| Protein-tag additive | `CCTTGGCACCCGAGAATT*C*C` |
| Hashtag additive | `GTGACTGGAGTTCAGACGTGTGC*T*C` |
| Guide-tag additive | `AGCAAGTGAGAAGCATCGTG*T*C` |
| Small-RNA RPI1 example | `CAAGCAGAAGACGGCATACGAGxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTC*C*A` |
| TruSeq D701_s example | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCAGACGTGT*G*C` |
| Guide-index Next_nst_x | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTATTTCTAGCTCTAAA*A*C` |

The protocol captions call the example tag barcodes 12 nt, but each printed example contains eight `x` characters. That mismatch remains unresolved; the model does not turn either into an exact barcode length. The same table also prints four TCR γ/δ enrichment primers. They belong to the optional receptor branch and are `R1_hTRDC` `AGCTTGACAGCATTGTACTTCC`, `R1_hTRGC` `TGTGTCGTTAGTCTTCATGGTGTTCC`, `R2_hTRDC` `TCCTTCACCAGACAAGCGAC`, and `R2_hTRGC` `GATCCCAGAATCGTGTTGCTC`.

## 4. Construct-changing steps

1. Label cells with antibody-derived tags and sample hashtags carrying distinct amplification handles.
2. Add a scaffold-complementary RT primer to the 10x 5′ RT mix.
3. Copy sgRNA through its variable protospacer and template-switch onto the bead-derived handle.
4. Transfer cell barcode and UMI to transcripts and all short-tag modalities.
5. Enrich gene-expression/V(D)J, ADT, hashtag and guide-derived oligo libraries separately.

## 5. Evidence boundary

🟢 The paper's author protocol supplies the oligo sequences and PCR recipes above. 🟡 The page emphasizes library-changing branches and represents experiment-specific barcodes as placeholders. The commercial 10x 5′ bead oligo is not given base for base in the ECCITE-seq paper; no exact bead sequence is inferred from these tag primers. It does not collapse ECCITE-seq into the later 5′ direct-capture Perturb-seq protocol already documented separately.
