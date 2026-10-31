# STAMP

Brannan et al. (2021), *Robust single-cell discovery of RNA targets of RNA-binding proteins and ribosomes*. DOI https://doi.org/10.1038/s41592-021-01128-0 .

## Source scope

🟢 Associated accepted-manuscript PMC8148648 overview/Methods124–202 inspected; final-VOR procedural applicability remains unverified. Published TableS11/MOESM11 is MD5-identical in the manifest to cached `Supplementary_Table_10.xlsx`; filename numbering is not its source identity. Its header and cloning oligos were inspected. Published reporting summary identity/date only inspected; unchecked OCR is not molecular evidence. Mandatory scraper ran on the actual bundle; five large data-workbook/text files were skipped by its default scan, not declared irrelevant or sequence-free. No absent oligo claim.

## Editing and branch boundaries

RBP/ribosomal-protein–APOBEC1 writes physical RNA C→U. Complementary firstDNA carries A at those selected sites; RNA/reference-oriented T calls are a separate representation. Constructors require explicit original C sites; they predict neither binding specificity nor editing efficiency. Ribo-STAMP is the ribosomal-fusion branch, not a separate method entry. Main56/124 explicitly describes removal/replacement of the DART YTH cassette to repurpose APOBEC1 for general RBPs. The local model owns this chemistry; it imports neither DART's methylation-targeting assumptions nor its library recipe.

Separate eCLIP validation140 is not an obligatory STAMP capture workflow. Source bulk totalRNA extraction is followed by either TruSeqStrandedmRNA or KAPARNAHyperPrep/RiboErase with500ngRNA. Kit names stay distinct despite source generic totalRNAseq wording. Historical internals stay opaque. Bulk acquisition144 is single-end100nt; eCLIP HiSeq platform values are not assigned to bulk STAMP.

Nanopore route: source20µgRNA polyA magnetic-resin selection (NEBE7490L),100ngpolyA RNA into SQK-DCS109 **direct cDNA** kit; PromethIONFLO-PRO002~48h/high-accuracyGuppy. It is not nativeRNA sequencing. PacBio route separately starts300ngtotalRNA into IsoSeqExpress cDNA; source13–15cyclePCR uses forward `GGCAATGAAGTCGCAGGGTTG` and reverse `AAGCAGTGGTATCAACGCAGAG`; SMRTbellExpress2.0 plus sample-specific BarcodedOverhangAdapterKit follows. Equimolar pooling, initial barcode-count rebalance and SequelII/SMRTcell8M/2.0chemistry/2.1polymerase/2hpre-extension/30hmovie are preserved. Do not force the Nanopore polyA frontend onto the separately described PacBio route.

Single-cell route202:72h1µg/ml doxycycline,TrypLE,count/resuspend1000cells/µl0.04%BSA/PBS, Chromium3′v3 GEM cell lysis/barcodedRT and cDNA amplification; mixedRBFOX2:TIA1 uses FeaturecDNAprimers1. Gene-expression fragmentation/end-repair/A-tail/adapter/sample-index continuation remains separate from the **saved Ampure cleanup supernatant** used for feature/template-switch-primer amplification and indexing. Other cleanup material remains unresolved, not ideally depleted/discarded. Both libraries use NovaSeq6000. Actual cellular/UMI/index sequences and run geometry are not invented from a modern kit. Feature-capture expression construct oligos in TableS11 are not independent RNAseq library adapters.

## Limits and ancestry

Five source preparation branches have explicit unresolved kit/end/strand geometry and one notes-only wrapper; no public finished schematic. Sample/RBP/long-read variants remain branches. Direct component-strategy ancestry is supported for DART→STAMP enzyme-construct repurposing (main56/124, YTH removal/replacement and shared source construct); it does not establish inherited library operations or equivalent targeting. TRIBE is compared as a conceptually related A→I method, not asserted a physical precursor to this library. Historical vendors and final-VOR reconciliation remain pending; no two-reader absence review or accepted missing-information claim.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 The associated accepted manuscript Methods136 explicitly places Agilent TapeStation quantification before 500ng RNA enters either named bulk library kit. The local bulk branches now separate analytical RNA material (amount unassembled), the actual 500ng library input, and untransferred RNA stock (amount and subsequent fate unassembled). This is an allocation trace, not a reconstructed TapeStation assay or an assumed quantitative recovery. The source bulk DNA library endpoint is typed as a library.

The required historical-library/final-VOR item remains uncovered. It cannot be satisfied by the endpoint label or a source-inventory list. Other branches retain their explicit limits: Methods184–192 requires actual independent technical preparations and subsequent physical pooling, not a one-product relabeling into replicates; Methods202 has separate expression and saved-supernatant feature preparations and library QA. Reconstructing these producer arrays and all historical kit interfaces is outside this bounded pass. Final-VOR applicability is unresolved, not a newly verified access failure. No public page or completion marker is introduced.

The closed source pass reread only associated manuscript identity1–27 and Methods124–144/172–202. Mandatory scraper pointers into other originals were not adopted as source facts; five default large-file skips do not establish absence. No original acquisition or missing-information audit was performed.
