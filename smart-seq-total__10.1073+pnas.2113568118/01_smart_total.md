# Smart-seq-total

Defining original: https://doi.org/10.1073/pnas.2113568118

Isakova et al.(2021), *Single-cell quantification of a broad RNA spectrum reveals unique noncoding patterns associated with cell types and states*.

## Sources read and versions

🟢 Published non-manuscript PMC8713755 article, Methods111–133, Results53–57 and limitations99. Original main PDF pages6–7 rendered/inspected for oligos and read conditions. Original SI PDF `pnas.2113568118_supp_pnas.2113568118.sapp.pdf`: Supplementary Protocol beginning (guide assembly/Cas9 digest/cleanup), FigS2/S3 captions439–469, and TableS2 original page37, image inspected and sequence/modification marks verified against extraction11491–11565. Full-bundle scraper run; library oligo hits verified in source originals. Result CSVs, remaining SI panels/tables and linked predecessor originals deferred. This is one capped pass, not a complete source/absence review.

## Two branches, native RNA and artificial tails

🟢 v1/v2 share single-cell well lysis, E.coli poly(A) polymerase addition and anchored oligo-dT reverse transcription with SSII/template switching. v2 also applies Vaccinia capping enzyme/GTP/SAM simultaneously during tailing, and uses a constrained UMI-containing primer. The graph explicitly declares these capping reagents but does not invent end eligibility, exact cap structure or conversion efficiency across cellular RNA classes. Circular RNAs are not silently opened or polyadenylated. Added tail length is unassigned. The same observed total A-tail length can arise from different endogenous/artificial contributions; the source lists loss of endogenous polyadenylation status as a limitation.

🟢 v1 oligo-dT: `5′-Biotin-CATAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGT30VN-3′`. v2 TableS2: `5′-biotin-TCGTCGGCAGCGTCAGTTGTATCAACTCAGACNNNHNNNHNNHNNHNNGCTGTAACT22VN-3′`, with the repeated Ts written out in the original. The16base constrained UMI uses H, not an unrestricted N16. It is a structured UMI feature. Biotin modification is retained as source annotation; molecular attachment geometry is not reconstructed.

## Modified TSO and glycosylase product

🟢 TableS2 v1 TSO: `5′-biotin-dUCGdUCGGCAGCGdUCAGdUdUGdUAdUCAACdUCAGACAdUrGrG+G-3′`; v2: `5′-biotin-AdUGGCdUCGGAGAdUGdUGdUAdUAAGAGACAGdUCdUrGrG+G-3′`. Local ordered residues retain DNA/deoxyU, terminal twoRNA-G and finalLNA-G identities. Copying makes complementary canonical DNA while the original modified TSO remains a separate molecule. RNA template/free-TSO populations are retained after RT; no undocumented RNA degradation step is inserted.

🟢 Directly cited NEB M0280 current product description was read as the third original solely for its catalytic relation: UDG is a monofunctional glycosylase releasing uracil by hydrolysis of the deoxyuridine N-glycosidic bond. Current lot formulation/unit-definition changes are not substituted for the2021 reagent. https://www.neb.com/en/products/m0280-uracil-dna-glycosylase-udg . The local selected path therefore creates AP sites and liberated uracil but preserves the TSO phosphate backbone. RNA-U is not silently treated as deoxyU. Subsequent backbone cleavage, fragment cleanup and AP-template traversal require explicit paths; source wording “TSO digestion” alone does not supply them. Reaction products stay molecular material before actual separation.

## Source notation and read-primer difference

The original main PDFp6/Methods117 prints v1 TSO with bare U at positions where TableS2p37 explicitly uses dU. Bare U in that DNA-oligo string can be shorthand for dU; the main text does not explicitly assign RNA sugar there. Thus this is notation ambiguity, not an established RNA-versus-DNA reagent conflict or publication error. The local v1 TSO object follows TableS2's explicit dU notation because glycosylase susceptibility depends on sugar identity; no assertion is made about the actual order sheet beyond those originals. A neutral graph-owned explanation preserves both written forms and this limit.

The main PDFp7/Methods133 specifies v2 custom Read1 `TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC`; TableS2p37 instead repeats the v1 custom primer `TCGGCAGCGTCAGTTGTATCAACTCAGACATGGG` for v2. Both original images were checked. The model declares no final sequencing geometry, so it does not silently select either value. Source main v1 readout is NovaSeq6000SE75/100 plus two12base indexes; v2 NextSeq50030+112 paired reads plus two12base indexes with a custom Index2 primer. Primer-role/acquisition reconciliation remains open.

## Downstream outline and lineage

🟢 Source cDNA PCR/cleanup is followed by direct indexing or a cited tagmentation branch, pooling, guide-directed Cas9 depletion of abundant library DNA (DASH), protease release and size cleanup. These are recorded but not assembled; depletion does not mean original cellular rRNA was absent. Main protease treatment50°C15min differs from the included DASH supplement56°C10min; no single reaction condition is silently selected for an unmodeled step. Exact index primers, complete PCR/adaptation geometry and selection/read output remain unresolved source items. Notes status; no finished schematic.

Main Introduction53 explicitly says inspired by Smart-seq2ref19; Results57 cites earlier total-RNA tailing/switching approaches21–22, and DASHref23. These support scoped strategy/component-lineage candidates, not wholesale inheritance. No predecessor constructors or family recipes are imported. Keep v1/v2 and direct-index/tagmentation routes as branches of this method, not artificial new catalogue entries.
