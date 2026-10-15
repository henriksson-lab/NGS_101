# scMethyl-HiC — single-nucleus source notes

Defining paper: Li et al., *Nature Methods* (2019), [doi:10.1038/s41592-019-0502-z](https://doi.org/10.1038/s41592-019-0502-z). This page covers the separately named **scMethyl-HiC** branch of the paper; [bulk Methyl-HiC](../methyl-hic__10.1038+s41592-019-0502-z/01_methyl-hic.html) has a different library-preparation path.

## Sources read and source scope

- 🟢 Original article, PMC6765429 JATS XML, Results paragraph on the single-cell adaptation (extracted line 67), sequencing analysis (line 113).
- 🟢 Original Springer Supplementary Protocol `41592_2019_502_MOESM1_ESM.pdf.txt`, oligo table lines 212–258; single-cell branch lines 640–795. The independent PMC copy `NIHMS1532780-supplement-5.pdf.txt` prints the same families near lines 59–105. The oligo scraper was run on the full source bundle; its family hits were verified in the original PDF text.
- 🟡 The article's single-cell Methods (extracted line 85) calls the conversion kit Zymo **D5020**, while the Supplementary Protocol reagent list (lines 202–204) calls it **D5023** and the procedure (lines 651–658) specifies the EZ-96 DNA Methylation-Direct deep-well format. The model records the shared kit family and does not choose a catalogue number; the authors' actual kit variant for each preparation remains unresolved.
- 🟢 [IDT xGen Adaptase Module successor guide](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/protocol/xgen-adaptase-module-protocol.pdf?sfvrsn=a5a7e007_10), 2022, local `xgen-adaptase-module-protocol.pdf.txt`: lines 150–164, 523–532 and 1235–1243 describe related tailing plus R2 stub attachment and an average eight-base synthetic tail. This later guide supports mechanism only; it does not establish the sequence or exact length of the 2019 Swift reagent product.
- 🟢 [Official 2021 Swift-to-IDT product notice](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/miscellaneous/swift-change-notification---090721.pdf?sfvrsn=a9f4e007_2), opening letter and Appendix A: maps Swift 33096 to IDT 10009826 and says formulation did not change during that branding transition. Both source reviewers followed and checked this product link; the notice does not print the proprietary adapter sequence or establish any earlier formulation history.
- 🟡 The linked 2019 Swift Accel-NGS Adaptase Module 33096 manual itself has not been obtained. Both source reviews therefore leave the historical junction unresolved.

### Source checklist and linked information

Two independent readers checked the original article, the Springer Supplementary Protocol,
its PMC copy, and each of the following local attachments: MOESM2 reporting summary
(brief single-cell workflow); PMC supplement 1 (figures), 2 (contact tables), 3
(per-cell metrics), 4 (compartment data), and 6 (background notes). The latter
five contain no additional oligo table or library construction details in the
sections inspected. Both readers also followed Swift product 33096 to the
official 2021 Swift-to-IDT mapping and checked the 2022 IDT successor guide.
The version-matched 2019 Swift manual remains unavailable, so the historical
Adaptase junction is **blocked**, not confirmed absent. Exact file locators and
the two located findings for every item are recorded in the repository's
`catalogue/source_reviews/scmethyl-hic__10.1038+s41592-019-0502-z.json`.

## Protocol path

1. 🟢 DpnII digestion, biotin fill-in and in situ proximity ligation generate contact junctions, as in the bulk branch (original article Results; supplementary protocol, in situ Hi-C procedure).
2. 🟢 **After proximity ligation**, single nuclei are sorted into individual 96-well bisulfite reactions (MOESM1 lines 640–660). This timing distinguishes the branch from bulk Methyl-HiC.
3. 🟢 Each well is bisulfite-converted before library construction. One of four P5L random primers, AD002/006/008/010, is used when four-plex pooling is planned (lines 665–668).
4. 🟢 Klenow exo-minus extends the random primers. Exonuclease I and shrimp alkaline phosphatase treat remaining primers, followed by SPRI cleanup. Four separately indexed plate reactions are combined if multiplexing (lines 671–730).
5. 🟢 Swift Adaptase Module 33096 treats the pooled single-stranded products; selected P5L and P7L primers perform indexing PCR (lines 735–790).

## Oligos and molecular boundaries

- 🟢 The four random primers are `/5SpC3/TTCCCTACACGACGCTCTTCCGATCT` followed by one of `CGATGT` (AD002), `GCCAAT` (AD006), `ACTTGA` (AD008), `TAGCTT` (AD010), then `(N1:25252525)(N1)(N1)(N1)(N1)(N1)(N1)(N1)(N1)` as written in the original table (lines 213–216). The first parenthesized synthesis symbol plus eight `(N1)` symbols is modeled as a random nonamer, not a UMI. The 5′ spacer is a modification, not a DNA base.
- 🟢 The source prints eight full P5L_D501–D508 and twelve full P7L_D701–D712 PCR primers (lines 218–258). For the representative D501 and D701, the model derives full oligos from canonical P5/P7 arms, the published eight-base index and their published annealing tails. D501 carries `TATAGCCT` after P5; D701 carries `CGAGTAAT` after P7. A source-text self-test checks **every** full family member against the original supplement.
- 🟡 The page shows exact source members AD002, D501 and D701 as examples, without asserting that this combination was used in a particular well. The P5 primer and copied random-primer handle overlap by 18 bases, giving a computable Read 1 site. The bisulfite-converted genomic contact insert is symbolic because its original DpnII `GATC` bases need not survive conversion unchanged.
- 🔴 The 2019 Adaptase-added P7 junction's exact sequence and length remain unresolved from the checked paper and supplement. The later IDT guide describes a variable synthetic tail of average length eight bases, but does not prove the historical Swift formulation. The model therefore uses an unknown-length inferred junction; complete Read 2/dual-index primer binding is not claimed.

## Readout

🟢 The paper uses paired-end contact reads; its single-cell analysis notes G-to-A treatment for Read 1 and C-to-T treatment for Read 2 (article extracted line 113). The source identifies P5L and P7L indexed PCR families and a four-way inline primer index, but does not assign one member combination to each nucleus in the checked material. The page computes only the representative Read 1 primer site; remaining read placement stays unresolved pending the historical adapter geometry and sample-index assignment.
