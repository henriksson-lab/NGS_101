# ECCITE-seq — research notes

## 1. What it is

ECCITE-seq extends a 10x 5′ workflow to recover transcriptome, clonotypes, surface proteins, sample hashtags and directly captured sgRNAs with a shared single-cell barcode.

## 2. Sources read

- 🟢 Mimitou et al., *Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells*, Nature Methods 16, 409–412 (2019), DOI [10.1038/s41592-019-0392-0](https://doi.org/10.1038/s41592-019-0392-0), open manuscript [PMC6557128](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557128/).
- 🟢 Methods and Figure 1 establish the bead-TSO-complementary antibody tags, separate protein/hashtag handles, scaffold RT primer, guide copying and template switching.
- 🟢 [Author supplementary protocol, pp. 5–6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0392-0/MediaObjects/41592_2019_392_MOESM3_ESM.pdf), **Oligonucleotide sequences** table, prints tag examples, amplification primers, guide RT primer and their terminal phosphorothioates (`*`). Supplementary Tables 3–5 separately list the experiment's antibody and hashtag barcode assignments.
- 🟢 The linked [author Protocol Exchange version](https://doi.org/10.1038/protex.2019.025), posted in May 2019, repeats the tag examples and details the 0.6× size split and separate ADT, hashtag, and guide PCR branches. It does not resolve the eight-`x` versus “12nt” caption discrepancy.
- 🟢 The [CITE-seq predecessor](https://doi.org/10.1038/nmeth.4380) uses a 3′ poly(A)-tailed ADT and a streptavidin–biotin disulfide conjugate. The [Cell Hashing predecessor](https://doi.org/10.1186/s13059-018-1603-1) instead uses covalent iEDDA conjugation, 12-base HTO barcodes, and separate ADT/HTO amplification handles. ECCITE adapts the latter's tag strategy to 10x 5′ capture; the CITE 3′ construct is not imported into the model.
- 🟢 The earlier 10x [CG000109 Rev D technical note](https://raw.githubusercontent.com/Teichlab/scg_lib_structs/master/data/10X-Genomics/CG000109_AssayConfiguration_VDJ_RevD.pdf), an original vendor PDF mirrored by scg_lib_structs, gives a 22-nt partial R1, 16-nt cell barcode, 10-nt UMI, 13-nt template-switch oligo and separate GEX/V(D)J layouts. The full CG000086 kit guide cited by that note and the exact revision used for ECCITE remain open.

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

The protocol captions call both example tag barcodes 12 nt, but each printed example contains eight `x` characters. That mismatch remains unresolved. The original Supplementary Table 3 (MOESM1 lines 219–387) prints 49 antibody barcodes of 15 bases and three of 12 bases; Tables 4–5 (lines 389–460) print 12-base mouse and human hashtags. The model draws representative exact barcodes from all three groups: Table 3's first `TAGAACCGTATCCGT`, one short Table 3 entry `AATAGCGGAGCC`, and mouse Hashtag_A `AGGACCATCCAA`. It does not infer the eight placeholder positions in the protocol examples to be the ordered conjugate barcodes. The same table also prints four TCR γ/δ enrichment primers. They belong to the optional receptor branch and are `R1_hTRDC` `AGCTTGACAGCATTGTACTTCC`, `R1_hTRGC` `TGTGTCGTTAGTCTTCATGGTGTTCC`, `R2_hTRDC` `TCCTTCACCAGACAAGCGAC`, and `R2_hTRGC` `GATCCCAGAATCGTGTTGCTC`.

## 4. Construct-changing steps

1. Label cells with antibody-derived tags and sample hashtags carrying distinct amplification handles.
2. Add a scaffold-complementary RT primer to the 10x 5′ RT mix.
3. Copy sgRNA through its variable protospacer and template-switch onto the bead-derived handle.
4. Transfer cell barcode and UMI to transcripts and all short-tag modalities.
5. Enrich gene-expression/V(D)J, ADT, hashtag and guide-derived oligo libraries separately.

## 5. Evidence boundary

🟢 The paper's author protocol supplies the oligo sequences and PCR recipes above, while its Supplementary Tables 3–5 supply ordered conjugate barcode examples. 🟡 The page emphasizes library-changing branches and uses three exact example barcodes while leaving final-library handle regions schematic; it does not assert their exact end-to-end orientation. The commercial 10x 5′ bead oligo is not given base for base in the ECCITE-seq paper. CG000109 Rev D prints an early 10x 5′ bead oligo, but the ECCITE kit revision has not been established, so its sequence is not imported as an ECCITE-specific fact. This page does not collapse ECCITE-seq into the later 5′ direct-capture Perturb-seq protocol documented separately.
