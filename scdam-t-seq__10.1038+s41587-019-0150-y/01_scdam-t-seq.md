# scDam&T-seq — research notes

## Sources read

- 🟢 Rooijers et al. 2019, DOI `10.1038/s41587-019-0150-y`, PMC full text and supplements.
- 🟢 Markodimitraki et al. 2020, DOI `10.1038/s41596-020-0314-8`, detailed protocol and supplementary oligo table.
- 🟢 Hashimshony et al. 2016 CEL-Seq2, DOI `10.1186/s13059-016-0938-8`,
  original article (predecessor of the RNA branch).

## Construct path

🟢 A Dam fusion deposits N6-methyladenine at nearby GATC sites in living cells. Cells are
sorted into wells containing a CEL-seq2 primer; its printed architecture is T7 promoter,
RA5, alternating 3-nt UMI / 4-nt cell-barcode parts, dT25 and anchoring V. Reverse
transcription is followed by second-strand synthesis in the same well.

🟢 DpnI cuts fully methylated GATC sites to blunt GA/TC half-sites. The detailed protocol
prints both strands of each DamID adapter. They contain a deliberately non-complementary
six-base fork, followed by a T7/RA5-bearing paired core, alternating UMI/barcode parts and
a terminal GA/TC blunt end. The bottom oligo is 5′ phosphorylated. Ligation of adapter
`GA` to genomic `TC` restores `GATC`; this exact junction and the unpaired fork are both
retained in the model. Supplementary Table 2 contains the complete 384-member adapter
family; the page uses the protocol's printed `TGCA`/`TATG` example.

🟢 The 2019 Rooijers Methods, “DamID adapters,” and its Supplementary Table 2
distinguish a four-base-fork adapter from the six-base-fork Dam-RING1B mESC variant.
The drawing selects the later detailed protocol's six-base-fork example; it should
not be read as the only 2019 adapter. The original CEL-Seq2 article describes a
six-base cell barcode in its own assay, whereas scDam&T-seq explicitly prints a
split eight-base barcode for its adapted primer.

🟢 CEL-seq2 and DamID use non-overlapping cell-barcode sets, but both precursors carry a
T7 promoter. Wells can therefore be pooled before joint linear amplification. The aRNA is
fragmented and reverse-transcribed using the exact
`GCCTTGGCACCCGAGAATTCCA-N6` randomhexRT primer. The printed RNA-PCR1 and indexed reverse
primer then install P5, the six-base library index and P7; their sequence overlaps define the representative endpoint. The native PCR input is single-stranded randomhexRT cDNA. Exact multiple-population first-cycle geometry is not supported by the current model; the page retains that actual input and marks this modelling limit instead of manufacturing a duplex or a complete RPi arm before PCR.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC6609448 JATS57,99–117; detailed protocol PMC7779467 JATS61,940–975,1102.

Public discrepancy notice compares 2020 overview transcript R1/R2 assignments with 2020 analysis instructions and 2019 Demultiplexing. Model follows RNA R1 identifiers/R2 transcript. Six-base-fork example remains explicitly the 2020 variant, not all 2019 preparations.

Deferred: Original variant oligo tables and supplementary manual; exact run and index acquisition for the selected 2020 example; overview discrepancy remains apparent rather than author-confirmed. This is a capped check, not a full source audit or confirmation of missing information.

2020 original Reagents (PMC7779467 JATS276–302) verifies the adapter strands and randomhexRT. An additional published discrepancy is now shown on the page: DamID oligos and the named barcode are TGCATATG, but the following split-layout example repeats RNA GATGTCAT. The model uses the printed DamID oligos, NNN-TGCA-NNN-TATG. The detailed DOI is registered for the source transcription guard.

The same 2020 original randomhexRT mix table (JATS878–884) lists a 20 µM primer stock and 13.33 mM final concentration, incompatible with dilution. The intended unit is unverified; no corrected concentration is assigned. This apparent published inconsistency is also shown on the public page.

## Bounded native graph migration (2026-10-18, batch 12)

This page selects the explicitly printed **2020 detailed-protocol example**, not a verified allocation to a particular 2019 experiment. Three originals were inspected: defining 2019 PMC6609448 Methods “Robotic preparation”, “DamID adapters”, “CEL-Seq primers” and “Raw data preprocessing” (lines 99–117); 2020 PMC7779467 Reagents (274–302), RT/second-strand procedure (568–630), fragmented-aRNA purification through PCR/library purification/sequencing (840–975); and 2020 MOESM1 supplementary-protocol DpnI/ligation section (text 275–330), which is bulk-method context and does not override the single-cell procedure. The mandatory bundle primer scrape was run; printed 2020 primer/adapter examples were checked directly. Superseded preprint sequences were excluded.

The source-owned graph replaces the generic `illumina_branch` TruSeq endpoints with endpoints assembled from the actual Small RNA RP1/RPi oligos. The printed RPi example index is `CGTGAT`, not an invented random six-base window. Both printed PCR primers have a terminal phosphorothioate bond (`*` before their final A); the ordering text and input polymers retain it. That physical primer bond stays on its incorporated original strand; its complementary copy does not inherit the modification. The DamID adapter bottom 5′ phosphate is likewise an input modification, not a copied base. Barcode and UMI parts have stable, distinct feature identities linked by branch-specific groups.

Native mRNA, genomic duplex DNA and the two derived molecular populations remain independent polymers. Transcript cDNA remains visible during DpnI digestion and DamID ligation; it no longer disappears from those states. Different indexed wells pool physically before joint T7 IVT; exact plate/member allocations remain unassigned. IVT creates RNA sequence copies and retains an explicit unresolved residual-DNA fate. Fragmentation includes identifier-bearing and other fragment classes. The source-stated 5-µL RT aliquot is allocated from the 12-µL eluate, and its remainder has an unresolved QC/remainder fate. Candidate other randomhexRT products are inferred and kept unresolved; they are not claimed physically removed. Computational demultiplexing/mapping filters are separate from bead-cleanup waste.

The actual incoming PCR state contains separate single-stranded RNA-derived and DamID-derived randomhexRT cDNAs. It contains the **21-nt** randomhexRT tail, not the complete 33-nt RPi arm. The rest of that arm is added by PCR. The representative final products retain the RNA poly(dT)/anchor or restored DamID GATC, the split identifiers, the RPi example index and source-derived P5/P7. Two final classes represent co-resident molecules, not separate physical tubes. Both final purifications are explicit.

The source describes NextSeq 500 paired-end 75-nt sequencing and branch-specific analysis. Exact sequencing-primer identity/run configuration is unassigned in this capped reading; no standard TruSeq R1/I1/I2/R2 shell is substituted. The existing published read-label, barcode-example and concentration discrepancy notices remain on the page. Pre-IVT and post-IVT bead reagent/dilution choices were not verified sufficiently in the inspected 2020 sections and stay unassigned; post-fragmentation and final 0.8-volume AMPure XP cleanups are source-stated.

Deferred: precise 2019 experiment-to-oligo mapping, exact pool composition, T7 initiation/opposite-end readthrough and fragment distributions, native multiple-population first-cycle PCR geometry, exact sequencing-primer/index acquisition setup, pre/post-IVT cleanup reagent/dilution, and detailed QC allocation of the fragmented-aRNA remainder. No source-absence claim is accepted, and existing blocked review status is unchanged.

### Code review corrections (2026-10-18; zero new originals)

`MolecularMixture` now preserves the actual source `Scene` components and exposes their namespaced semantic strand/partner inventory. Rendering uses those same objects; mixture construction no longer collapses them to rows. Fragmentation retains unplaced branch-specific barcode/UMI segment inventories from the actual incoming aRNA components. The inventory records identity somewhere in the fragmented population and does not assign all features to every fragment.

The main physical path retains possible original-RNA hydrolysates/residual RNA after RNase H and through DpnI/ligation, possible residual DNA templates and RNA hydrolysates in the IVT pool, the fragmented input RNA plus possible other cDNAs after randomhexRT, and those possible original ss RNA/DNA/nonselected classes in the PCR product and cleanup eluates. The former separate residual-DNA/other-cDNA terminals were removed: they did not establish a physical allocation away from the main tube. Only actual waste streams and the source-stated RT aliquot/remainder are separate allocations. These co-carried classes are inferred possibilities without measured abundance or claims of complete cleanup removal. Both final readout-class materials retain the mixed physical state; their separate endpoint illustrations identify co-resident sequence classes, not purified populations.
