# methyl-ATAC-seq

Spektor et al. (2019), DOI https://doi.org/10.1101/gr.245399.118 . Private capped source-curation note; no finished schematic or accepted missing-information claim.

Published PMC6581052 main Methods102–116 and published Supplemental Tables S1–S3 / Figures S1–S6 inspected in the current pass. Table S3 rendered original page1 confirms individual internal methyl-dC tokens, including a line-wrapped token; ordinary source text extraction is not treated as a chemical discrepancy. Main112 specifies methylated Tn5 adapters, methyl-dCTP end repair, a saved ten-percent QC aliquot before bisulfite and a separate post-conversion PCR. Full source pointers, limits and local model will follow without claiming complete source review.

## Positive chemistry and local scope

🟢 The active library uses250000 source cells and methylated A/B adapters; lambda is an unmethylated conversion-control spike. Source lysis is explicitly the local modified-Omni procedure; inactivation is SDS/proteinase K followed by isopropanol/GlycoBlue recovery, then0.5x large-fragment bead removal/1.8x final recovery. This explicit local inheritance does not establish the original Omni paper's whole recipe. ATAC comparator Q5 PCR is separate from the methyl-ATAC library.

🟢 End repair uses Klenow exo-minus plus5-methyl-dCTP. After1.8x cleanup, ten percent is saved for QC PCR; the rest proceeds to historical EZ DNA Methylation-Lightning D5030T conversion, immediate PfuTurbo Cx13-cycle PCR and1.8x cleanup. Qubit is analytical; its assayed material is separated from production provenance. Source NextSeq500 acquisition is single-end75bp. Computational methylation extraction excludes terminal9nt at both ends; it does not chemically remove those nucleotides.

🟢 Table S3 source oligos (verbatim modifications; wrapped internal token rejoined only after rendered-original inspection):

```
Tn5MErev /5Phos/CTGTCTCTTATACACATCT
ME-A_mC T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG
ME-B_mC GT/iMe-dC/T/iMe-dC/GTGGG/iMe-dC/T/iMe-dC/GGAGATGTGTATAAGAGA/iMe-dC/AG
i7 CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG
i5 AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC
```

The exact study index members and lengths remain unassigned; the model attaches sample-index meaning to each placeholder. Known nucleotide bodies are imported from canonical Nextera/Illumina constants and source-token identity is checked. Synthetic protective adapter/repair methylation is distinct from biological methylation. The selected physical5mC-only relation produces uracil at unprotected C, and a complementary DNA descendant separately; it is not a universal statement about every bisulfite-resistant modification or conversion efficiency.

## Limits and reproducibility

The local graph covers source-known operations through acquisition with source-specific choices. It does not solve strand-specific tagmentation gaps, complete converted template/primer geometry, kit-internal conversion/recovery or all read cycles. `tools/chem.py` provides selected chemical classes and `tools/selftest.py` tests copying, mark constraints, saved-QC provenance and source transcription. No private preview or finished page is claimed. The source guard checks the drawn sequence states, not all source coverage.

Directly named Picelli2014 Tn5 production, Adey/Shendure2012 assembly, Corces2017 Omni, Lu2015 T-WGBS, Suzuki2018 X-WGBS and historical manufacturer originals remain deferred at the cap. The downloaded published table and figure PDFs were inspected as selected above; superseded preprint is not active authority. `ref/MANIFEST.md` and the source-review JSON contain document locators, exact scopes, actual failed inventory URLs and open targets. No two-reader missing-information audit or accepted absence conclusion.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 Published main104 specifies independent A/MErev and B/MErev anneals, followed by the actual25ulA/25ulB/50ulglycerol preparation and equal-volume Tn5 contribution. The graph now allocates reverse-strand contributions before separate annealing, retains both source-paired duplex references and the protected internal methyl-dC tokens, and selects the actual7ul contribution from the prepared loaded stock. Its copresent catalyst/adaptor/other classes feed tagmentation. The original free Tn5 input is not added again after loading. Exact occupied transposome distributions and free companions remain unassembled; paired reference duplexes do not establish a single universal occupied complex.

SDS lysis buffer and actual100ug proteinaseK are distinct contributions. End repair receives actual5UKlenowexo-minus,0.5mMeach dATP/dGTP/dTTP/5-methyl-dCTP and buffer. Later PfuTurboCx PCR supplies have unreviewed composition and do not inherit protective5mdCTP or the Omni Q5 recipe. Known final DNA-library kind survives cleanup/Qubit allocation; single-end NextSeq75 is structured while index members/cycles remain unassigned.

The saved ten-percent preconversion QC material and postconversion Qubit material remain separate from production. The original pending SourceItem stays required and uncovered. Historical Picelli/Omni/Zymo/Pfu references, whole converted substrates, source end/read geometry and complete diagnostic/QC details remain unresolved. No source coverage is generated and no ready page/builder is manufactured for this incomplete model. Source table wrapping was checked against the original and is not a publication error.
