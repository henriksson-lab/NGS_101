# ACE-seq — APOBEC-coupled epigenetic sequencing

## Source read

- 🟢 Schutsky *et al.*, “Nondestructive, base-resolution sequencing of 5-hydroxymethylcytosine using a DNA deaminase,” *Nature Biotechnology* (2018), DOI `10.1038/nbt.4204` (defining paper and open author manuscript, PMCID: PMC6453757).

## Molecular interpretation

- 🟢 DNA and methylated-lambda/5hmC-containing T4-phage controls are sheared to about 300 bp; Supplementary Protocol 2 recommends adding controls before shearing.
- 🟢 T4 β-glucosyltransferase converts 5hmC to glucosyl-5hmC. After heat denaturation and snap cooling, APOBEC3A deaminates C and 5mC while protected 5hmC resists conversion.
- 🟢 PCR and sequencing therefore read original C and 5mC as T, while original 5hmC remains C. No bisulfite treatment is used.
- 🟢 The defining whole-genome experiment prepared converted DNA with the Swift Accel Methyl-NGS kit and used single-end NextSeq sequencing.
- 🟢 Publisher Supplementary Protocol 2, “Quality control on phage spike-ins,” prints four PCR primers for the methylated-lambda and 5hmC-containing T4 controls. These assay conversion controls; they do not construct or sequence the library.
- 🟡 Exact kit adapters and sequencing primers remain unassigned pending applicable vendor-original review. Their role-level display does not establish that these sequences are unpublished.

## Final readout

ACE-seq has the same base-call polarity as TAB-seq but achieves it by selective enzymatic deamination after glucosylation, avoiding destructive bisulfite chemistry.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC6453757 JATS56–64,110,118,130.

Corrected deamination intermediate to U at original C/5mC positions and retained T as PCR readout. Recorded single-end NextSeq500 and 150-cycle kit without inventing acquisition mask. Withdrew unverified proprietary adapter-end claims.

Deferred: Original supplementary protocols/QC table and vendor Accel manual; exact adapter installation and acquisition mask. This is a capped check, not a full source audit or confirmation of missing information.

## Native-graph migration (2026-10-18)

Read PMC6453757 JATS Methods 110–118 and publisher MOESM6 Supplementary
Protocol 2 Materials, sample preparation, conversion, library construction and
spike-in QC (lines 1298–1534). The previously read manuscript supplement-1
library/QC passage (1290–1505) was checked against the publisher original;
the publisher supplement governs this migration. Three originals inspected.
The bundle scraper returned zero general hits; targeted `--find` located the
lambda primer in both supplements, and all four primers were read directly in
the publisher table. No absence follows from the zero-hit result.

The protocol now owns its base-state views, role-level library and page fields.
The graph includes controls before shearing, separate heat/snap-freezing and
A3A steps, and source-stated Zymo recovery before the delegated kit procedure.
Kit library construction, qPCR-selected amplification and final purification
remain one explicit unresolved boundary: no adapter-end geometry or actual PCR
precursor is assigned without the vendor original. The former shared helper's
18-base adapters and 26-base insert were schematic defaults, not source values;
role tokens now preserve unknown lengths. Deamination produces U; T is the
amplified readout.

The main Methods (JATS 118) specify **0.5% total by mass** for spike-ins, while
publisher Supplementary Protocol 2 Materials specify **0.5% of each** control.
The public page reports both values. The graph declares their presence without
choosing a quantity; the per-experiment amount remains unresolved.

Deferred: the applicable Swift Accel-NGS Methyl-Seq 30024 and indexing-kit 36024
vendor originals, exact intermediate/adapter/primer assignments, and actual
sequencing acquisition mask. Existing blocked missing-information status is
unchanged; this is a bounded migration, not a complete chemistry review.
