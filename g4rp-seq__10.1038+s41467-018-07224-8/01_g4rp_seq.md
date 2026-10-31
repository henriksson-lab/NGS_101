# G4RP-seq

Yang et al.(2018), *Transcriptome-wide identification of transient RNA G-quadruplexes in human cells*. DOI https://doi.org/10.1038/s41467-018-07224-8 .

## Original sources and scope

🟢 Published non-manuscript PMC6226477 provenance36; Results60–76, limitations116, Methods138–148 inspected. Published SI MOESM1 identity/selected BioTASQ characterization pp1/5–6, crosslinked-control Fig6p11, analysis Fig7p12 and in-vitro Table3p22 inspected. Those labeled FRET/validation oligos are not sequencing-library adapters. Mandatory full-bundle scraper run; source modifications inspected where relevant. SI source-data/other workbook originals remain unreviewed; one capped pass, no absence audit.

## Source-owned baseline route

🟢 Source baseline MCF7 cells are fixed1%formaldehyde5min25°C, quenched0.125Mglycine5min, scraped/resuspended in G4RP buffer and sonicated(Covarism220,10%duty2min). The buffer contains150mMKCl. Separate BRACO19/RHPS4-treated cell conditions and biotin negative controls stay separate experimental branches; the graph does not assume baseline cells received a G4-stabilizing drug.

BioTASQ100µM overnight4°C associates with RNA-containing complexes. Its biotin supplies streptavidin affinity; this is not enzymatic/covalent biotinylation of the RNA nucleotide. Source5%sonicate input-control aliquot is separate from the pull-down; its exact handling point relative probe exposure is not fully reconciled in this capped pass, so the local split is bookkeeping for an unresolved control route, not a biochemical exposure conclusion.

Streptavidin magnetic beads10µg2h4°C capture complexes, followed by four G4RP-buffer washes5min. The bead fraction is heated70°C1h for reverse crosslinking, then TRIZOL RNA extraction. Heating alone is not physical removal of beads/probe/protein, nor proof all crosslinks reverse. The graph retains these materials until extraction and gives actual separated fractions unquantified recovery.

Main Methods148 explicitly orders DNase treatment, brief thermal fragmentation, IlluminaRiboZero rRNA removal and IlluminaTruSeq RNA Library Prep. Exact historical kit editions/internal strand steps/oligos are unassembled. RiboZero is an opaque kit boundary: the graph retains an unresolved processed preparation and RNA/DNA/remnant fates, without assuming affinity capture, a specific internal reaction or immediate physical disposal. No modern RiboZero/TruSeq recipe is substituted. Sequencing is NextSeq500PE75; exact adapter/index/primer geometry remains unresolved. Input-control and pull-down libraries are compared computationally, not physically merged by normalization.

## Physical relation and interpretation

🟡 `SelectedPartition` represents explicit caller-supplied molecular IDs and a selected binding outcome, conserving identities between retained/flowthrough and eluted/residual fractions. It predicts no G4 affinity, capture specificity, equilibrium, crosslink success or recovery fraction. Native RNA coordinates are not deleted because read enrichment changes.

🟢 Main limitations acknowledge probe-induced stabilization/parallel-G4 preference and formaldehyde capture of indirect RNA–RNA interactions. Fixed RNA structures and published enrichment signals do not prove every recovered RNA directly contacts a G4 or that every G4 is recovered. The graph keeps these limits rather than adding a universal G-rich binding rule. The paper explicitly describes G4RP as a modified RIP strategy; rG4-seq is a cited alternative with a distinct RT-stall frontend, not an inherited complete library.

Open work: exact full capture/polymer/adduct topology, input-control timing, historical kit originals, terminal adapters/index/strand/read geometry. No full finished page or confirmed missing-information claim.

## Bounded unlinked-page adjudication (2026-10-30)

🟢 Main140 specifies formaldehyde crosslinking and glycine quenching as actual reagent contributions. The local graph now represents those molecules separately from the source time/temperature conditions, with a typed quenching operation. The actual working G4RP buffer exposes its co-present KCl, Tris, EDTA, DTT, NP40, RNase inhibitor and homebrew protease-inhibitor-cocktail classes from one physical buffer aliquot; no extra doses or undisclosed cocktail components are invented. Working-buffer volume and all component allocations remain unassembled. TRIZOL is typed as a molecular extraction reagent rather than a selection matrix.

The original unreviewed_kit SourceItem remains uncovered. Main148 names historical Ribo-Zero/TruSeq preparation and specifies two sequencing replicates for each condition, with an input for each condition. It does not justify manufacturing replicate libraries by relabeling one reacted product, nor does it establish two input replicates. The current five-percent input-control producer and full library arrays remain unassembled obligations. Exact input timing/exposure stays qualified; no alternate exposure is asserted from the current split location alone.

Main144's SSIII/random-hexamer/SYBR qPCR route is separate analytical validation and is not substituted for the NGS kit recipe. BioTASQ remains a noncovalent probe; source heating is neither physical bead removal nor proof of complete reversal. No runtime source inventory, opaque endpoint or full-method equivalence is used to satisfy the unresolved kit/read-geometry item. The closed pass inspected published main140/148/limits116 and selected published SI contexts only; no new originals, recipe defaults or finished public page are introduced.
