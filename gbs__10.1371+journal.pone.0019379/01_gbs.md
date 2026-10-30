# GBS — original genotyping-by-sequencing

Elshire et al.(2011), *A Robust, Simple Genotyping-by-Sequencing (GBS) Approach for High Diversity Species*, defining DOI https://doi.org/10.1371/journal.pone.0019379 .

## Inspected originals

🟢 Published PMC3087801 main Methods79–109, Fig1/2 captions and Discussion195–203; published TextS1 (`pone.0019379.s001.docx`, optimization/adapter-titration sections); published TableS1 (`pone.0019379.s004.docx`, complete96barcode list). Three originals substantively read within the capped pass. Mandatory scraper rerun after acquisition; full strings, lowercase placeholder notation and absence of adapter5′phosphate checked against those originals. FigureS1/S2 alias downloads failed; these uninspected electrophoresis figures are not claimed to omit chemistry. No recursive predecessor/vendor chase or two-reader absence audit.

## Source adapter chemistry

🟢 GBS uses two different, fully complementary adapter types, not a divergent Y adapter. All four ordering strands are **unmodified without5′phosphate**, as explicitly stated in Discussion197. ApeKI makes5′cohesive ends compatible with source `CWG`; W denotes A/T. The private component view selects CAG, one allowed member, without assigning its fraction or predicting which genomic sites were cut. Methylation sensitivity/sequence alone does not establish quantitative cleavage.

The barcode top has canonical `illumina.TRUSEQ_READ1` plus4–8nt sample barcode. Its bottom is `CWG` plus reverse-complement barcode and reverse-complement canonical body. CTCC, a4nt member from published TableS1, is selected for the component example; no experimental sample/well assignment is inferred. The common-adapter bottom is `C` plus canonical `illumina.PE_READ2_ARM`; top is `CWG` plus its reverse complement. PCRA is canonical `illumina.NEBNEXT_UNIVERSAL_PRIMER`; PCRB is `illumina.P7` plus `CGGT` plus the common-adapter bottom. Whole printed equations, common strands and PCR bodies are original-dependent checks that skip when originals are absent.

Three precisely scoped sourceguard exemptions are explicit: two CTCC barcode-adapter strands assembled from the main's equations and TableS1, and the common-top CAG member of the printed CWG class. They are source-defined component expansions rather than contiguous literal resolved strings. Strict global matching remains unchanged; no synthetic source twin or claimed study-specific molecule/ratio. Independent review verified the complete barcode equations,CTCC TableS1 membership,CAG as a valid CWG member,both common strands and both PCR primers,with no study-specific molecule/ratio inferred. This critical implementation review is not an absence audit. The common bottom and both PCR primers are checked directly.

🟢 Selected ligation seals the genomic5′phosphate–adapter3′OH junction; the opposing adapter5′OH–genomic3′OH junction remains unsealed. The local `NickedFragment` represents that selected class, not universal ligation success. Its selected PCR descendant is a different object: copying does not ligase-seal the original5′OH nick. Original templates and displaced/unligated adapter strands persist physically until actual recovery. Unknown whole genomic ends, displacement/copy products and exact run-primer geometry remain unassembled.

Barcode/common, barcode/barcode and common/common end-type classes remain possible; both orientations of the heterotypic class are retained. Input1:1adapter mixture is not a predicted50%yield or population distribution. Source homotypic flowcell/readout limitations are historical qualitative statements, not inferred physical purification of homotypic molecules or a universal modern-instrument rule.

## Operational route and fractions

🟢 Barcode and common stocks are annealed separately, quantified through separate PicoGreen analytical aliquots, mixed1:1 and dried into wells before genomic DNA is added and dried. Thus adapters are already present during ApeKI digestion75°C2h. T4 ligase/ATP is added afterward for22°C1h;65°C30min inactivation follows. The source explicitly distinguishes this sequential single-well route from simultaneous digestion/ligation because the enzyme temperatures differ. No fresh adapters are introduced after the digestion in the local graph.

Five microlitres from each50µl well enter the48/96plex pool, leaving independent remainders with unassigned later fate. Main105 specifies QIAquick PCR purification and50µl elution; Fig2 caption instead calls the step a size-exclusion column to remove adapters. Both source wordings are retained: the local model selects explicit main105, without asserting cartridge equivalence or treating the caption as a proven publication error. TextS1 titration cleanup also names Qiagen PCR cleanup; its separately optimized200ng/adapter-titration conditions are not copied into the100ng multiplex route.

Two microlitres of pooled DNA enter Taq PCR, with an independent eluate remainder. Initial72°C5min supports source fill/copy before98°C and18cycles; Discussion197 also allows displacement during assembly/early heating, so the model predicts no universal fill-time or efficiency. New copies, original nicked templates and all end-type backgrounds remain distinct. PCR purification/30µl elution precede a separate1µl Experion aliquot. The size profile/dimer assessment determines acceptance or reconstruction; it is **not a fragment-size excision**. Source170–350bp profile and∼128bp dimer observations do not become exact molecule lengths or a selection filter.

## Readout and provenance limits

🟢 Main109 reports GenomeAnalyzerII single-end86nt sequencing and inline4–8nt sample tags; no separate indexing read is added. Later64nt computational trimming/polyA padding, genotype filtering and imputation do not physically modify DNA. Discussion203's source explanation of homotypic limited-readout/cluster behavior remains separate from the local selected-copy relation and does not establish yields.

Main71/Discussion191 explicitly contrasts and develops the earlier Baird2008 RAD barcoding strategy, while removing its random shearing, size excision and divergent-Y path. This supports a scoped developmental relationship, not import of that molecular recipe. Historical Illumina/Bentley2008 and source product guides remain deferred. The private graph/source class model is incomplete; no finished schematic, accepted source absence or complete chemical-validation claim.
