# ECCITE-seq — research notes

## 1. What it is

ECCITE-seq extends a 10x 5′ workflow to recover transcriptome, clonotypes, surface proteins, sample hashtags and directly captured sgRNAs with a shared single-cell barcode.

## 2. Sources read

- 🟢 Mimitou et al., *Multiplexed detection of proteins, transcriptomes, clonotypes and CRISPR perturbations in single cells*, Nature Methods 16, 409–412 (2019), DOI [10.1038/s41592-019-0392-0](https://doi.org/10.1038/s41592-019-0392-0), open manuscript [PMC6557128](https://pmc.ncbi.nlm.nih.gov/articles/PMC6557128/).
- 🟢 Methods and Figure 1 establish the bead-TSO-complementary antibody tags, separate protein/hashtag handles, scaffold RT primer, guide copying and template switching.
- 🟢 [Author supplementary protocol, pp. 5–6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0392-0/MediaObjects/41592_2019_392_MOESM3_ESM.pdf), **Oligonucleotide sequences** table, prints tag examples, amplification primers, guide RT primer and their terminal phosphorothioates (`*`). Supplementary Tables 3–5 separately list the experiment's antibody and hashtag barcode assignments.
- 🟢 The linked [author Protocol Exchange version](https://doi.org/10.1038/protex.2019.025), posted in May 2019, repeats the tag examples and details the 0.6× size split and separate ADT, hashtag, and guide PCR branches. It does not resolve the eight-`x` versus “12nt” caption discrepancy.
- 🟢 The [CITE-seq predecessor](https://doi.org/10.1038/nmeth.4380) uses a 3′ poly(A)-tailed ADT and a streptavidin–biotin disulfide conjugate. The [Cell Hashing predecessor](https://doi.org/10.1186/s13059-018-1603-1) instead uses covalent iEDDA conjugation, 12-base HTO barcodes, and separate ADT/HTO amplification handles. ECCITE adapts the latter's tag strategy to 10x 5′ capture; the CITE 3′ construct is not imported into the model.
- 2026-10-17 bounded pass: PMC6557128 JATS text (Results para. 1–2, Methods lines 96–114), MOESM3 supplementary protocol (RT line 104, PCR recipes lines 205–310, oligo list lines 383–429), MOESM1 Supplementary Tables 3–5 via the self-test. Not re-read: Protocol Exchange PDF, 10x CG000086 user guide (not on disk), CITE-seq and Cell Hashing predecessors.
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

🟢 Methods, "ECCITE-seq on 10× Genomics instrument" (PMC6557128 XML text lines 102–114):

1. Stain cells with antibody-derived tags and sample hashtags carrying distinct amplification handles.
2. Load the 10x Single Cell V(D)J workflow; spike 12 pmol scaffold RT primer (gd_RT_v4; MOESM3 line 399 calls it sg_RT_v4) into RT when guides are wanted.
3. In the GEM, tag 3′ ends (`CCCATATAAGAAA`) anneal to the bead TSO and are extended; transcripts and sgRNAs are reverse-transcribed and template-switch onto the bead oligo.
4. Amplify cDNA with spiked hashtag (1 pmol), guide-tag (1 pmol) and protein-tag (2 pmol) additive primers.
5. 0.6X SPRI: mRNA cDNA stays on beads (10x V(D)J transcript and TCR α/β libraries; TCR γ/δ by two nested PCRs with the printed R1_/R2_hTRDC/hTRGC primers); tag/guide products stay in the supernatant.
6. Bring the supernatant to 2.0X SPRI, wash, then a second 2.0X SPRI to remove excess single-stranded oligos.
7. Separate KAPA HiFi PCRs: SI-PCR + RPI-x (protein tag, ~180 bp), SI-PCR + D7xx_s (hashtag), SI-PCR + Next_nst_x (v4) or smRNA_nst_x (v3) (guide). "x" is the i7 sample index. MOESM3 (Version March 2019) step 7 first enriches guides with SI-PCR + gd_add_v4 (8–10 cycles, ~220 bp; gd_add_v4 is not in the oligo list and is presumed to be the guide-tag additive), then indexes with Next_nst_x (~210 bp); tag libraries are purified at 1.6X SPRI. Methods call the extra guide enrichment optional.

## 4a. Final libraries (derived)

🟡 Assembled from the printed oligos (computed in `tools/eccite_seq.py`; bases checked against MOESM3 by the self-test):

- Common P5 side: SI-PCR (`P5` + partial TruSeq Read 1), then the bead oligo's last nine Read 1 bases, 16-nt cell barcode, 10-nt UMI and the 13 TSO bases complementary to the tag 3′ end. The bead bases and barcode/UMI lengths are not printed by ECCITE sources (CG000109 Rev D, revision not established) and are drawn inferred/placeholder.
- Protein tag: barcode, then `revcomp` of the RPI1 handle (`TGGAATTCTCGGGTGCCAAGG…`, small-RNA RA3), i7, P7 without the terminal `AT` (RPI1 lacks it). Index 1 and Read 2 primers are derived from that handle; the run's primer mix is not stated.
- Hashtag: barcode, then the TruSeq Read 2 adapter printed in the hashtag oligo, i7, P7 (D701_s).
- Guide: protospacer, the scaffold 5′ bases annealed by Next_nst_x (`GTTTTAGAGCTAGAAATA`), Nextera Read 2 adapter, i7, P7.
- Transcript: 10x V(D)J-protocol library, not inspected here; drawn as an inferred single-index TruSeq shell after the TSO bases.

The model previously placed tag handles on the P5 side and used TruSeq Read 2 arms for all tag and guide libraries; that was a repository modelling error, corrected 2026-10-17.

🟢 Reads (Methods): 26 cycles Read 1, 8 cycles index, Read 2 39 cycles (HiSeq 2500 rapid run) or 91 cycles (NovaSeq 6000 S2). The model declares Read 1 and Index 1 lengths and leaves Read 2 unfixed.

## 5. Evidence boundary

🟢 The paper's author protocol supplies the oligo sequences and PCR recipes above, while its Supplementary Tables 3–5 supply ordered conjugate barcode examples. 🟡 The page assembles tag and guide final libraries end to end from the printed oligos (section 4a) with three exact example barcodes; bead-oligo bases, barcode/UMI lengths, the transcript library adapters and the protein-tag sequencing primers remain inferred or derived. The commercial 10x 5′ bead oligo is not given base for base in the ECCITE-seq paper. CG000109 Rev D prints an early 10x 5′ bead oligo, but the ECCITE kit revision has not been established, so its sequence is not imported as an ECCITE-specific fact. This page does not collapse ECCITE-seq into the later 5′ direct-capture Perturb-seq protocol documented separately.
