# i-BLESS — research notes

Biernacka et al., *Communications Biology* (2018), [doi:10.1038/s42003-018-0165-9](https://doi.org/10.1038/s42003-018-0165-9), PMCID: PMC6208412.

## Sources read

- 🟢 Full text and Methods from the PMC article.
- 🟢 Supplementary Information from the PMC Cloud copy (`42003_2018_165_MOESM1_ESM.pdf`), especially Supplementary Table 4.
- 🟢 Rubicon Genomics, *ThruPLEX DNA-seq Kit Manual*, QAM-108-003 (revised September 2016), Appendix 1B and Fig. 9, [archived vendor PDF](https://www.takara.co.kr/file/manual/pdf/ThruPLEX-DNA-seq-Kit-Manual-QAM-108-003-updated-9_16_16.pdf).
- 🟡 [Takara's later ThruPLEX adapter FAQ](https://help.takarabio.com/36168/kb/article/110112/what-are-the-index-adapter-sequences-after-preparation-using-thruplex-dna-libraries) has a locally fetched 2023-updated landing page. Its substantive content/source-child identity requires reconciliation; the title and navigation alone are not evidence for historical 2018 reagent geometry.

## What makes this i-BLESS rather than BLESS

1. 🟢 Yeast cells are embedded in agarose beads before spheroplasting, lysis and protein removal. The beads permit reagent diffusion while protecting very small yeast nuclei and long DNA from handling-induced breaks.
2. 🟢 Break ends are blunted in the beads and ligated to the same biotinylated proximal hairpin used by BLESS.
3. 🟢 Encapsulated DNA is first sonicated, recovered from agarose and sonicated again to about 400 bp. Streptavidin capture therefore still separates the pre-labelled native end from a later fragmentation end.
4. 🟢 The captured free end is blunted, phosphorylated and ligated to the BLESS distal hairpin. I-SceI linearization, PCR and XhoI digestion retain the proximal and distal 11-nt tags.
5. 🟢 Final sequencing libraries were prepared with the ThruPLEX DNA-seq Kit and run 2×75 bp on HiSeq 2500 and HiSeq 4000.
6. 🟢 The 2016 kit guide says repaired DNA receives a blunt-ended, 5′-blocked stem-loop adapter; extension and cleavage of that adapter occur during PCR. Its single-index kit has an **8-nt Sanger index**. Appendix 1B gives all 12 index *bases* (for example, Tube 1 `ATCACGTT`), and Fig. 9 places the index at the library's 3′ end.
7. 🟡 Study-specific ThruPLEX product, revision, index member and exact final geometry remain open assignments. The 2016 guide also offers dual-index products (Appendix 1D), so its single-index example is not a universal study configuration. The schematic now uses an unassembled role-level endpoint and a conditional dated-guide reference, while retaining the exact public i-BLESS hairpins and insert-PCR primers. No missing-information claim is accepted.

Supplementary Table 4 prints the proximal and distal adapters, their modifications, and the two insert-amplification primers. It explicitly states that the adapters are identical to those of original BLESS.


## Source-owned reconstruction checkpoint

The model now declares the workflow locally. It no longer imports the generic
`dsb_capture` endpoint/PCR recipe or a generic TruSeq construction.

- 🟢 Methods labeling paragraph (article text 131): sonicate encapsulated DNA,
  recover with Zymoclean, sonicate again to about 400 bp, capture on streptavidin,
  blunt/phosphorylate and ligate distal hairpin, I-SceI linearize, PCR amplify,
  purify PCR products, then XhoI digest. The earlier graph omitted the first recovery
  and insert PCR; these were repository modeling errors, not publication errors.
- 🟢 Supplementary Table 4 (text 217–239) supplies both whole hairpins with 5′
  phosphates, proximal internal biotin-dT, and both whole insert-PCR primers.
  All four are now molecular graph inputs. Exact intermediate duplex assembly and
  primer-binding coordinates remain unassembled; the model does not substitute the
  final sequencing library as the initial insert-PCR template.
- 🟢 BLESS original Methods (PMC3651036 text 131–141) establishes related hairpin
  capture/release chemistry. Its actual Table 5 is cached as
  `nmeth.2408_supp_NIHMS449608-supplement-6.xls`, not `supplement-5.xls`.
  Its L1/L3/P1/P3 agree with the i-BLESS table. Its HaeIII fragmentation,
  nuclei preparation and TruSeq endpoint are not transferred to i-BLESS.
- 🟢 September 2016 ThruPLEX guide Principle/Workflow (294–327): repaired blunt DNA,
  blocked-5′ stem-loop ligation leaving a genomic-3′ nick, extension/cleavage and
  indexed amplification. This conditional reference replaces the earlier generic
  dA/Y-adapter diagram. Whole study kit geometry stays unassembled.
- 🟢 The guide distinguishes an eight-base reagent index from acquired cycles:
  the TruSeq LT IEM example reads only the first six bases (1458–1465).
  The page declares only the study's positive **2×75** genomic read lengths,
  leaving its index-cycle configuration unassigned.
- 🟢 Methods primary analysis (139) trims and selects proximal-tagged reads
  computationally. That is not physical purification of a molecular pool.

## Actual linked-source scope remains open

The original calls Mitra *et al.* 2015,
[10.1371/journal.pone.0120520](https://doi.org/10.1371/journal.pone.0120520),
for low-diversity acquisition. Its actual Methods (383–457) concerns image-template
cycles, cluster analysis and read preprocessing, with a separate HiSeq 2000/TruSeq
library example. That example does not establish i-BLESS ThruPLEX reagent identity.
Its Instant Sequencing page and cited historical Illumina technical documents remain
actual linked scope to resolve.

The identified Overhauser 1992 encapsulation predecessor
([10.1385/0-89603-229-9:129](https://doi.org/10.1385/0-89603-229-9:129))
has no retrieved full original in the current registered route. The original BLESS
supporting page `http://www.breakome.eu` and i-BLESS software links
`http://breakome.utmb.edu/software.html` and `http://breakome.eu/software.html`
remain actual targets, alongside deposited study SRP125409. Remaining original
supplement visual scope and vendor-linked documents have not passed a complete
independent two-reader gate. These are concrete obligations, not evidence that an
unassigned study value was unpublished.
