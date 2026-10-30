# ini-seq

Defining DOI: https://doi.org/10.1093/nar/gkw760.

## Actually inspected sources

🟢 Published PMC5137433 Fig1/main54 and Methods62–80; published File010SI FigS1–S8 and TableS1(qPCR) scope. Whole-bundle mandatory primer scraper run; hits are reading aids and not independent evidence. Superseded preprints, when automatically fetched, are not read or used.

## Scope and boundaries

Selected DIG-labelled cell-free origin initiation, not EdU/BrdU live-cell timing. Standard assay conditions/nuclei preparation defer to cited sources; source ChIP-seq sample-prep kit version and molecular library/read geometry remain unreviewed. qPCR validation primers inTableS1 are not NGS adapters. No complete origin specificity or uniform DIG incorporation assumed.

## Readout

Source HiSeq2000; exact acquired read roles/cycles/index members unassigned. Enrichment near origins is not a direct physical read of replication initiation nucleotide.

The private graph records source-known operations, while unresolved exact molecular geometry is an explicit model boundary. This capped pass is not a complete source review or a missing-information audit.

## Source-known route and material retention

🟢 MainMethods64/68: source lateG1EJ30 nuclei are prepared under named predecessor methods and incubated15min in the cell-free replication assay with proliferatingHeLa extract andDIG-11-dUTP. Ice-coldPBS stops the reaction; centrifugation recovers nuclei. Sarcosinate buffer/proteinaseK2mg/ml55°C18h→phenol/chloroform→RNaseA200µg/ml37°C1h→secondphenol/chloroform preserves physical peptide/RNA products until their actual extractions. Sonication yields100–1000bp, then centrifugation removes debris; this is the source input pool.

🟢 AntiDIG/ProteinGDynabeads2hRT, sixLiCl/Igepal/deoxycholate washes and twoTrisNaCl washes precede bead transfer and twice15min65°CSDSeleution. Input andIPDNA receivephenol/chloroform andethanol recovery, then the source historicalIlluminaChIP-seq kit andHiSeq2000. The graph represents input/IP separately, and does not confuse this with qPCR validation. The kit boundary remains opaque: no end-repair/adapters/PCR/run defaults are borrowed merely because a kit is calledChIP-seq.

DIG is an affinity chemistry marker, not a molecular identifier. New nascent DNA remains separate from parental templates. Source records incomplete initiation and elongating-background nuclei, so no ideal specificity/completeness is assumed.

## Direct predecessor targets

MainMethods64/68 cite ref27Krude2000,10.1074/jbc.275.18.13699 (cell-free initiation/nuclei conditions) and ref28Krudeet al.1997,10.1016/s0092-8674(00)81863-2. These exact originals are deferred; matching replication assays do not establish a whole-library ancestry edge toOK-seq/Repli-seq. The named historicalIlluminaChIP-seq sample-preparation kit/version remains an unreviewed original target; no absent-sequence claim.
