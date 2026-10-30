# HAC-seq: m3C-selective cleavage with ordered demethylation controls

Cui et al. (2021), **Nucleotide resolution profiling of m3C RNA modification by HAC-seq**, [published defining article](https://doi.org/10.1093/nar/gkaa1186), [published PMC original](https://pmc.ncbi.nlm.nih.gov/articles/PMC7969016/).

Three originals inspected: published main Methods65/77/81/85 and Results95/103/Fig2caption107, including the rendered original Fig2 on PDFpage6; published supplementary workbook opening/header and selected cleavage-ratio rows; published supplementary DOCX caption text. The latter has a false-positive/validation figure caption, not a full supplementary protocol; its embedded plots were not adopted as inspected molecular evidence. The workbook is a results table, not an oligo table. All three count toward the cap. No source-era manufacturer original or cited predecessor was opened. Mandatory primer scraper hits and modifications were reviewed: tRNA isolation/Northern probes are separate assays, not the operational sequencing adapters. Their existence is not evidence that kit primers are missing or absent.

The preparation boundary differs from existing RNA chemical-cleavage methods in its high-salt m3C-selective treatment, pre/post-cleavage AlkB controls and adapter-compatible end selection. Source main95 describes hydrazine reactivity at m3C and U but m3C preference under3MNaCl. This is not a new protocol merely for a different tRNA target or computational cutoff. It is not modeled by importing HydraPsiSeq's end identities or AlkAniline-Seq's substrate chemistry.

🟢 The three groups in original Fig2 are important:

| Group | AlkB before fragmentation/HAC | Chemical HAC | AlkB after HAC |
| --- | --- | --- | --- |
| Ctrl | Yes | No | No second treatment shown |
| HAC | No | Yes | Yes |
| DM-HAC | Yes | Yes | No second treatment shown |

The named HAC group receives post-cleavage demethylation; Fig2 resolves the otherwise broad Methods phrase “HAC samples.” Each group is an independently prepared biological replicate, not three names for one physical RNA pool. AlkB-WT and AlkB-D135S are external reagents; their production is outside this library graph. Demethylation is variable, not assumed complete for every modification/site.

🟢 Upstream main65 describes TRIzol total RNA extraction and RiboMinus depletion using10µg totalRNA. The private graph starts from the independently prepared RiboMinus RNA stock, with its recovery and upstream kit chemistry unassembled. Main77 gives the total-RNA pre-demethylation scale as5µg in100µl, whereas each fragmentation allocation uses2µg. The graph preserves the remaining post-treatment stock without assuming perfect recovery or that one depletion produces all group aliquots.

🟢 Selected source preparation:

1. For Ctrl/DM-HAC, treat the RNA stock before cleavage with AlkB-WT160pmol, AlkB-D135S200pmol, RNasin and the source demethylation mix at room temperature in the dark for2hours. Quench with5mMEDTA and purify with acid phenol/chloroform. The2x buffer stock, not its final reaction concentration, is0.6MKCl/4mMMgCl2/0.1mM(NH4)2Fe(SO4)2·6H2O/0.6mM2-ketoglutarate/4mM ascorbate/0.1mgmlBSA/0.1MMESpH5.0. The HAC group omits this pre-treatment.
2. Allocate2µg RNA, fragment in10µlAM8740 reagents at70C10minutes and stop with20mMEDTApH8.0. Ethanol precipitate. Source nominal fragments are about200nt, not a perfect length distribution.
3. Apply Antarctic phosphatase, then T4PNK following manufacturer protocols, and purify with ZymoRNA-Clean&Concentrator5. Main103/Fig2 reports5′phosphate/3′OH ends before HAC. Exact starting fragment ends and enzyme-specific substrate transitions are unassembled. This is **before** chemical cleavage; no extra end repair is silently inserted afterward.
4. HAC/DM-HAC receive25µl ice-cold10percent hydrazine/3MNaCl on ice4hours, followed by ethanol precipitation. Resuspend the pellet in100µlH2O:glacialacid:aniline7:3:1 and incubate2hours at room temperature in the dark. Ethanol precipitate again and redissolve in water. Ctrl omits both chemical stages.
5. Only the named HAC group then receives post-cleavage AlkB treatment, using160pmolWT+200pmolD135S in100µl for2hours at room temperature in the dark, with source quench/purification. This preserves the control order in Fig2. No guarantee of perfect modification removal is asserted.
6. Allocate100ng of the recovered fragments from each group to historical NEBNextSmallRNA library preparation. The source specifies adapter ligation, RT and PCR and a local RT change to60C1hour. Exact kit revision, oligo modifications, junction geometry, cycle program and physical off-product recoveries are unassembled; the compound graph boundary does not assert simultaneous mixing of all kit stocks.
7. Sequence on IlluminaNextSeq500. Actual read lengths, sample indexes/index-cycle mask and complete primer-binding geometry are not established in this pass.

🟢 Source Fig2 states that chemical cleavage yields an adapter-compatible3′fragment with5′phosphate/3′OH, while the5′fragment has a damaged3′end. The exact damaged terminal group is not assigned from another protocol. Uncleaved RNA and blocked/compatible fragments coexist. Lack of adapter ligation is not physical destruction or an actual size selection of every blocked fragment: those populations remain until a source workup, and their detailed kit fates remain unassembled.

Source-demethylation resistance and false-positive validation limit interpretation; the graph predicts neither quantitative cleavage nor absolute m3C occupancy. Cleavage-ratio/site calling, reference alignment and significance thresholds are computational steps. The analysis follows m7G-TRAC-seq in the source; that does not authorize inheriting its wet chemistry or treating this as a generic kit-only variant.

This is a source-qualified private partial model with three explicit group constructors. Main and published supplements are available. Remaining source-era kit/reagent-production and exact molecular/readout questions are ordinary unfinished review, not actual unavailable documents or confirmed missing disclosure. No full source/absence audit or public schematic is claimed.
