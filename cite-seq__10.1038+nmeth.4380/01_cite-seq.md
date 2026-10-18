# CITE-seq — research notes

## 1. What it is

CITE-seq converts surface-protein abundance into sequenceable antibody-derived tags while measuring the transcriptome from the same single cell.

## 2. Sources read

- 🟢 Stoeckius et al., *Simultaneous epitope and transcriptome measurement in single cells*, Nature Methods 14, 865–868 (2017), DOI [10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380), open manuscript [PMC5669064](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669064/).
- 🟢 Main Methods, especially the **Antibody-oligo sequences** subsection, and [author supplementary protocol, pp. 6–7](https://static-content.springer.com/esm/art%3A10.1038%2Fnmeth.4380/MediaObjects/41592_2017_BFnmeth4380_MOESM3_ESM.pdf), establish the poly(A)-tailed antibody tag, DNA-dependent extension during RT, size separation and independent ADT library PCR. The protocol has an additional v2_BC7 oligo absent from the paper's listed CBMC subset.

## 3. Published antibody oligos

🟢 For the 10x Single Cell 3′ **v2** CBMC experiment, the author protocol prints 14 ordered oligos with `/5AmMC12/` at the 5′ end, the common `CCTTGGCACCCGAGAATTCC` TruSeq small-RNA Read-2 handle, a shared `A` followed by a **six-base** antibody identifier, one degenerate `B`, and a poly(A) tail. The seven-base `A` + identifier regions in BC1–BC14 order are `AATCACG`, `ACGATGT`, `ATTAGGC`, `ATGACCA`, `AACAGTG`, `AGCCAAT`, `ACAGATC`, `AACTTGA`, `AGATCAG`, `ATAGCTT`, `AGGCTAC`, `ACTTGTA`, `AAGTCAA`, and `AAGTTCC`. All 14 lengths were checked against the printed sequence rows; the defining paper prints 13 of these.

🟢 A complete example (v2_BC1) is `/5AmMC12/CCTTGGCACCCGAGAATTCCAATCACGBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA`. The antibody identifier is the six bases `ATCACG` in Supplementary Table 2 (PMC supplement-1, lines 45–55); the preceding `A` is shared by the v2_BC1–14 family (Springer author protocol, lines 403–434). The terminal `B` and poly(A) tract are retained as ordered; they are not antibody identifiers. The paper separately prints Drop-seq and 10x 3′ v1 species-mixing designs with Nextera Read-2 and Read-1 handles, respectively. Those are different branches of CITE-seq, not replacements for the v2 CBMC construct.

🟢 The author protocol also prints both ADT library PCR oligos for the 10x v2 branch: PE 1.0 P5 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T` and Small RNA RPI1 `CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. The two `*` positions in P5 are terminal phosphorothioate bonds. It separately gives a different P5-SMART-PCR hybrid oligo for the Drop-seq branch.

🟢 The paper's Online Methods (lines 94–132 of the PMC JATS text) name the two species-mixing oligos as BC6 and BC12 for both Drop-seq and 10x 3′ v1, with different platform-specific PCR handles. Its Supplementary Table 2 (PMC supplement-1, lines 178–196) assigns `GCCAAT`/BC6 to human CD29 and `CTTGTA`/BC12 to mouse CD29. These assignments apply to the species-mixing experiments, not to the 10x v2 CBMC panel.

🟢 The paper calls its 10x v2 CBMC assay a **13-antibody panel** (JATS line 70), prints the v2 family under “CBMC profiling” (lines 106–132), and maps the antibodies to identifiers in Supplementary Table 2 (lines 45–176). Crosswalking the table to the author's full v2 sequence list (protocol lines 403–434) resolves the panel as CD3e→v2_BC1, CD19→BC2, CD16→BC3, CD4→BC4, CD45RA→BC5, CD8a→BC6, CCR5→BC8, CCR7→BC9, CD11c→BC10, CD34→BC11, CD56→BC12, CD14→BC13, and CD10→BC14. The author's v2_BC7 (`ACAGATC`) is an additional printed oligo; Table 2 assigns no panel antibody to this identifier. The paper's 13 listed v2 oligos likewise omit BC7. This crosswalk identifies the CBMC panel without assuming that every oligo in the author protocol was used.

🟡 The four CD8a FACS-sort comparison pools used an oligo-labelled CD8a antibody and Drop-seq (JATS lines 68 and 156), and Table 2 gives the CD8a identifier `GCCAAT`. The Methods do not name the complete handle-bearing oligo separately for each pool. Assigning the CBMC v2_BC6 construct to those Drop-seq validation runs would therefore be an inference.

## 4. Construct-changing steps

2026-10-17 bounded pass: PMC JATS Online Methods lines 82–144 and author protocol (MOESM3, version 2017-07-11) lines 7–40, 196–335, 382–442 were re-read. The earlier model separated cDNA and ADTs by size *before* reverse transcription and built the ADT library as a dA-ligated TruSeq DNA library; both are corrected below.

1. 🟢 5′-amine/C12 oligos are biotinylated with EZ-link NHS S-S biotin (disulfide, cleavable in reducing conditions) and bound to LYNX streptavidin-conjugated antibodies (Methods line 84; protocol lines 7–40).
2. 🟢 About 500,000 cells are Fc-blocked (TruStain FcX, 10 min), stained 30 min at 4 °C and washed three times (Methods line 136; protocol lines 196–210).
3. 🟢 The 10x Single Cell 3′ v2 run follows the manufacturer's instructions; the bead oligo(dT) primer is extended by the reverse transcriptase on both mRNA and the DNA antibody tags (Methods lines 92, 144).
4. 🟢 cDNA is amplified for 10 cycles; 0.6X SPRI then keeps cDNA >300 bp on beads and ADTs (~180 bp) in the supernatant (Methods line 144; protocol lines 211–222). 🟡 The native antibody-tag strand/end state through the delegated 10x cDNA PCR remains unassigned by this model.
5. 🟢 The ADT supernatant is purified by two sequential 2X SPRI rounds to reduce primer carryover, which otherwise yields a ~140 bp TSO-RT-oligo product (protocol lines 239–268, Figure 3a).
6. 🟢 ADT PCR: KAPA HiFi, PE 1.0 P5 and a TruSeq Small RNA RPIx primer, 8–12 cycles, then 1.6X SPRI; the library is 170–180 bp (protocol lines 269–317). The Methods state 10 additional ADT cycles for the runs (lines 140 and 144), within the protocol range.
7. 🟢 The transcript library follows the 10x v2 manufacturer instructions (Methods line 144).

🟡 Derived ADT library (computed from the printed oligos in `cite_seq.ADT`): PE 1.0 P5 – TruSeq Read 1 – cell barcode – UMI – poly(dT) – complement of the ordered B (V) – antibody barcode reverse complement – complement of the shared A – reverse complement of `{{= cite_seq.ADT_HANDLE }}` – the remaining RPI arm – six-base RPI index – P7. The poly(dT) length is a vendor bead property, drawn as an inferred placeholder.

## 5. Evidence boundary

🟢 The defining paper and author protocol print the antibody oligos and ADT PCR primers. 🟡 This model has not assigned acquired read cycles (only “HiSeq 2500 Rapid Run with v2 chemistry”, ADT 10% / cDNA 90% of a lane, Methods line 140; protocol lines 330–335) or selected Index 1 / Read 2 sequencing primers for the small-RNA arm; the model therefore assigns no cycle counts and declares only TruSeq Read 1 on the ADT library. The precise antibody-barcode extraction/read-primer mapping remains an open model target. The complete 10x bead oligo and the cDNA amplification primers are separate vendor dependencies (CG00052 not re-read in this pass).

## Native graph migration, 2026-10-18

Read two published originals: article Online Methods (PMC5669064 JATS, antibody sequences/staining and Drop-seq/10x sections, 82–144) and the author protocol, printed version 2017-07-11 (MOESM3, staining and fractionation 196–268, ADT PCR/cleanup 269–317, sequencing 330–335, modified oligos 403–444). Ran the full bundle oligo scraper; verified its relevant pointers against those originals. The superseded preprint and additional original/vendor documents were not read this pass.

The graph now owns ADT amplification, 1.6X post-PCR cleanup and sequencing declarations locally. It retains the paper's 10 additional cycles within the author protocol's 8–12-cycle range; PE1.0 P5 retains both terminal phosphorothioate bonds, and RPI1 is a printed representative of the RPIx indexed-primer family. The first copied ADT is explicitly single stranded and retains the copied small-RNA handle. The preceding cDNA PCR receives the actual independent transcript/ADT copies, and the post-cDNA-PCR ADT strand/end state remains unassigned. Its final PCR therefore uses that incoming state with unavailable geometry, never the finished indexed ADT duplex.

Independent RNA, protein and ADT populations use separate Scene rows. Empty/unresolved droplets have unresolved recovery, rather than an invented physical-discard fate; staining washes and SPRI supernatants are physical removal steps. The ADT final cleanup is explicit, including its source 1.6X ratio, ethanol washes and eluate. Both transcript and ADT readout plans remain separate. The former shared transcript endpoint had silently installed a dA junction, an eight-base sample-index window and a full transcript/poly(dT) assembly. Those unsupported defaults are removed: the v2 manufacturer library-preparation route is an explicit bounded operation with unassigned native adapter/index/primer geometry. No generic PCR finisher stands in for that multistep route.

This capped migration does not establish information absence. The applicable vendor route, native pre-ADT-PCR ends and exact run/selected small-RNA primer assignments remain review/model targets.
