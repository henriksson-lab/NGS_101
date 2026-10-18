# Long-read ChIA-PET

## 1. What it is

Long-read ChIA-PET enriches chromatin contacts associated with a chosen protein, joins paired genomic fragments through a bridge linker, and converts the products into longer paired-end tags.

## 2. Authoritative source

🟢 Li et al., *Long-read ChIA-PET for base-pair-resolution mapping of haplotype-specific chromatin interactions*, Nature Protocols (2017), DOI `10.1038/nprot.2017.012`.

The PMC full text, including the oligo list and stepwise procedure, was read.

## 3. Construct-changing operations

1. 🟢 Dual-cross-linked, sonicated chromatin is immunoprecipitated for the target protein.
2. 🟢 Fragment ends are repaired and dA-tailed.
3. 🟢 Bridge linker-F is `/5Phos/CGCGATATC/iBIOdT/TATCTGACT`; bridge linker-R is
   `/5Phos/GTCAGATAAGATATCGCGT`. Their 18-bp core leaves a 3′ T at each end.
4. 🟢 One linker bridges two dA-tailed chromatin fragments in a single ligation reaction.
5. 🟢 Nextera Tn5 tagmentation runs for 5 min at 55 °C.
6. 🟢 Streptavidin captures tagmented fragments retaining the internal biotin linker;
   on-bead indexed PCR uses no more than 13 cycles.

🟢 The amplification table names Nextera NPM, PPC and separate i7 and i5 index primers, but does not identify the index wells used for each library or print those ordered primer oligos. The canonical Nextera end structure in the schematic is supported by the named kit; a particular sample's index bases remain unspecified.

## 4. Final library and reads

The final library has canonical Nextera sequencing ends around `genomic tag — bridge —
genomic tag`. Paired 150- or 250-base reads approach the bridge from the outer Tn5 ends.
The paper accepts either bridge orientation in a read, trims it, and maps the two flanking
genomic tags as the contact pair.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC5537732 JATS88–92,120,222,714.

Added crosslink reversal/DNA release before tagmentation. Exact bridge bases and internal biotin position agree with original. Candidate Tn5 products now unresolved rather than presumed discarded. Qualified contact inference and recorded alternative paired 150/250-base reads.

Deferred: Detailed source steps outside inspected passages; exact selected sequencing/index mask; fuller library states and tagmentation products. This is a capped check, not a full source audit or confirmation of missing information.

## Batch02 native-graph migration (2026-10-18)

The generic Nextera constructor is removed. A locally declared conditional reference uses canonical sequence constants only for the source-named FC-121-1031 Nextera DNA kit; physical PPC/index primer members and index widths remain unassigned. These inferred segments, the final caption and the page caveat distinguish that reference from a reconstructed historical kit lot. Biological2×150/2×250 source alternatives remain separate from unassigned index acquisition.

Independent chromatin loci now occupy separate rows until true bridge ligation. Actual source dA-prepared fragment ends accompany the known double-ended bridge geometry. Crosslink reversal remains before tagmentation. Actual tagged/captured contact precursors, with unresolved gap/strand geometry, feed PCR; the expected final indexed duplex is not used as template. Source PPC cocktail, initial72°C3min fill/copy phase and the generally≤13cycle ceiling are explicit. Native source Zymo post-Tn5 recovery, Ampure post-PCR cleanup and300–600bp BluePippin selection are now explicit. Non-bridge fragments are removed by physical bead washes; unresolved Tn5 product classes are not declared discarded.

One original was read: published PMC5537732 JATS identity/overview92, source kit/linker reagents190/222, and procedural49–72 (55°C5min Tn5, Zymo recovery, M280 blocking/capture, PPC/i5/i7 PCR setup, initial extension, cleanups and sequencing). Mandatory bundle scraper hits were followed into the published source; no predecessor or newer kit revision was substituted. Full first-cycle copying/gap geometry, kit/index member/width and exact index acquisition remain deferred. This capped pass does not accept source absence or full molecular validation.
