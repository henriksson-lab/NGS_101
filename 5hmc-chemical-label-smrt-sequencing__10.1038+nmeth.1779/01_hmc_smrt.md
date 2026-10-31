# 5hmC chemical-label SMRT sequencing

Defining DOI: https://doi.org/10.1038/nmeth.1779

## Sources read and scope

🟢 `nmeth.1779_PMC3646335.xml.txt`: Associated manuscript NIHMS465877; Results40–54; OnlineMethods68–84/94. Actual genomic labeling/capture/disulfide reduction/workup and historical A-tailing SMRTbell/PacBioRS route; synthetic/control branches kept separate.

🟢 `nmeth.1779_41592_2012_BFnmeth1779_MOESM18_ESM.pdf.txt`: Published supplement DOI/title identity; Fig1chemical classes and Table3 source8309–8374/footnotes. Synthetic specificity/SMRTbell/validation sequences inspected in their named scopes; no assignment of synthetic hairpins/overhang recipe to genomic libraries.

Capped pass: at most10minutes/3originals, published supplements preferred over associated-manuscript aliases where available. No linked-reference recursion, two-reader absence audit or finished-protocol promotion.

## Source chemistry and model

🟢 Selected genomic route: E14Tg2A mESC DNA sonicated100–500bp;200µgDNA,100µMUDP-6-N3-Glc,2µMβ-GT in50mMHEPESpH7.9/25mMMgCl2,37°C1h. Bio-RadMicroBioSpin6 cleanup precedes150µM dibenzocyclooctyne-disulfide-biotin37°C2h, then Qiagen nucleotide-removal cleanup.

🟢 MyOneStreptavidinC1 binding/washing uses100µlbeads/B&W plus0.01%Tween20. The source directly introduces the disulfide biotin linker to avoid interference with the SMRT sequencing protocol's separate polymerase-immobilization biotin. This is a new chemical reporter/recovery component, not ordinary original hmC-Seal recovery.50mMDTT100µl25°C2h gentle rotation cleaves the disulfide. Released DNA retains HS-N3-5gmC, a thiol-linker/azido-glucose adduct; it is not restored to native5hmC or converted into thymidine. Beads and cleavage remnants remain physically together until actual magnetic segregation. The graph explicitly collects the DNA supernatant separately from retained beads, then BioSpin6 removes DTT and MinElute/EB recovery follows. No perfect binding/release/recovery is claimed.

🟢 Genomic SMRTbells use the source-cited A-tailing method (Travers2010ref14, DOI10.1093/nar/gkq543), average insert~200bp. Actual PacBioRS small-insert sequencing is modeled as known acquisition with unresolved historical primer/hairpin/end/library geometry. The original describes circular consensus repeated interrogation of the same molecule and strand-specific IPD signatures. A predicted5hmC coordinate is a computational interpretation of a retained adduct's polymerase kinetics, not a substituted base or guaranteed pause at every chemical modification. Context and stochastic pauses remain relevant.

### Synthetic branch is separate

Published Table3 gives specificity oligos with `x` denoting5mC or5hmC, synthetic SMRTbell oligos with `x` denoting5hmC, and separate glycosylation/restriction-validation PCR primers. Main84 says the synthetic oligos are5′phosphorylated, annealed with four-base overhangs, T4-ligated25°C60min, treated with exonucleaseIII/VII37°C30min and QIAquick-cleaned. These are source facts for the synthetic branch, **not** a genomic hairpin sequence/ligation recipe. No exonuclease digestion or synthetic-end sequence is imported into the selected genomic graph. Thus no synthetic sequence is falsely made a genomic adapter candidate merely because it is available.

### Component ancestry and remaining scope

Main40/46/Methods68/94 explicitly uses the earlier Song hmC-labeling component (DOI10.1038/nbt.1732); its azido-glucosylation provenance is a component-level ancestry candidate to hmC-Seal. New disulfide-click/DTT chemistry and genomic SMRTbell acquisition are local to this original, not equal to the parent recipe. Source-cited Travers2010 A-tailing and historical PacBioRS standard protocols are ordinary deferred originals; their chemistry/geometry remains unreviewed. Published SI was available despite the challenged PMC alias; this is not a blanket claim that the chemistry supplement is blocked. No accepted absence or full corrected-VOR equivalence claim.

## Actual access barriers and deferred originals

- https://pmc.ncbi.nlm.nih.gov/articles/instance/3646335/bin/NIHMS465877-supplement-SI.pdf — PMC supplementary file; reCAPTCHA challenge page; the publisher's ESM files were fetched, so probably redundant

Manifest wording such as “probably redundant” is a fetcher hint, not an identity/equivalence conclusion. Exact inspected scopes above determine authority; ordinary deferred historical geometry is not an access failure.

## Validation

Protocol-owned `tools/chem.py` has a zero-argument graph and focused state/fate checks. It is an incomplete graph-only implementation, with no public page or completion marker. Source-known late operations/acquisition are retained; unassembled adapters, modifications and fates are not inferred from other families.


## Bounded publication review — 2026-10-30

The actual genomic labeling reaction now receives explicit β-GT catalyst and UDP6N3Glc donor. DBCO-disulfide-biotin and DTT are physical molecular contributions, distinct from their source incubation programmes. Reduction retains the DNA-bound thiol-label remnant, followed by actual magnetic separation and source cleanup; it does not restore native 5hmC or itself discard beads. Original `pending` remains unchanged and uncovered: genomic A-tailing SMRTbell construction per Travers2010 and historical PacBioRS standard-protocol applicability are unreviewed. Published SI Table3 synthetic hairpins, four-base overhangs, 5′ phosphates and exonuclease steps are not assigned to genomic libraries. The source 200µg labeling context and 160µg labeled input/1.3µg pulldown recovery context remain distinct; no allocation between those quantities or invented physical pool is supplied.

No additional originals were read or acquired during implementation. This capped pass accepts no missing-information claim and does not publish the graph.
