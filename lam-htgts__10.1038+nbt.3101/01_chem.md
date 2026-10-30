# LAM-HTGTS: linear amplification and bridge capture

Frock et al. (2015 issue; online2014), **Genome-wide detection of DNA double-stranded breaks induced by engineered nucleases**. DOI[10.1038/nbt.3101](https://doi.org/10.1038/nbt.3101).

This private partial graph explicitly separates two source revisions/scopes: associated author-manuscript PMC4320661 supplies selected LAM operational Methods107, while publishedMOESM10 Table8 independently supplies modified bridge oligo inventory. Full published-main equivalence is unresolved. The older end-repair/A-tail emulsionHTGTS recipe in the adjacent paragraph is not merged into the LAM branch; the source says LAM/bridge bypasses those operations. The later2016 detailed LAM protocol is one follow-up, not a new entry and not an automatically inherited recipe.

## Inspected original sections and oligos

🟢 Associated manuscript selected LAM paragraph107 and adjacent older-route paragraph105 only for boundary comparison. 🟢 Published SI Table8 bridge rows, originalPDFp30 rendered `/tmp/lam19-bridge.png`, plus primer-name context. Mandatory scraper418 candidates; selected whole modified strings were verified against the original and have direct transcription checks.

- ssA-upper-TN: `GCGACTATAGGGCACGCGTGGTNNNNN-/3AmMO/`
- ssA-lower: `/5Phos/CCACGCGTGCCCTATAGTCGC-/3AmMO/`
- ssA-lower-A: `/5Phos/ACCACGCGTGCCCTATAGTCGC/3AmMO/`

The source common core is21nt; upper carriesTNNNNN after it and a3′amino blocker. Both lower alternatives retain5′phosphate and3′amino blocks. Published Table8 assigns lower to A-tailed DNA and lower-A to G/C/T-tailed or blunt DNA. These are end-conditioned inventory alternatives, not both strands ligated onto every substrate. Actual substrate-end distribution, adapter allocation and complete junctions remain explicitly unassembled. The partially random splint is not a molecular barcode/UMI.

## Source operational graph

1. Entry is isolated sonicated genomic DNA containing cellular bait/prey rearrangement junctions. Those junctions already arose in cells; the graph does not claim a new in-vitro chromosomal joining reaction. Genomic geometry and LAM-specific size allocation are unassembled.
2. A single5′biotin primer and Taq perform50linear-amplification cycles. Additional Taq is added for another50cycles. These are newly copied biotin-bearing DNA strands, distinct from original genomic DNA and unused primer; no one-cycle PEM-seq extension recipe is substituted.
3. Source MyOneC1 streptavidin binding precedes overnight on-bead bridge-adapter ligation. Exact wash/denaturation/strand recovery details remain unreviewed, so the retained population is coarse rather than assumed pure ssDNA. End-conditioned adapter allocation has an explicit unknown boundary and unassigned stock remainder.
4. Source selected ligation conditions are15%PEG8000/1mMhexamminecobaltchloride,25°C60min,22°C120min,16°C8–12h. Beads and other reaction components remain until a real source separation; no invented elution is inserted.
5. Nested-locus15cyclePCR adds sourceI5/I7-related sequence tails. EcoNI blocking digestion follows; only recognition-site-containing molecular classes are selected as cut, site-free classes remain intact and digestion fragments are **not physically discarded**. The source's intended suppression of germline amplification is not a perfect-removal claim.
6. P5/P7 step-outPCR10–15cycles reconstructs the MiSeq sequence tags. New daughter DNA does not inherit the manufacturing biotin/amino-blocker chemistry. Production cleanup/QC/loading remain an explicit unassembled handoff after the last known PCR.

Source acquisition is MiSeq paired-end; selected Methods do not yet establish modeled full primer landing/run lengths. Mapping and duplicate filtering are computations, not physical separations.

## Pending sources and ancestry

Defining published full Methods/version applicability, directly cited LAM-PCR and bridge-ligation predecessors, complete bait/template geometry and actual adapter allocation remain pending. Main Methods explicitly attributes these two molecular components, supporting component-level ancestry once exact parent identity/version is verified; it does not establish complete recipe identity with any existing PEM-seq, GUIDE-seq or nrLAM entry. No accepted absence or public completion claim follows this capped pass.
