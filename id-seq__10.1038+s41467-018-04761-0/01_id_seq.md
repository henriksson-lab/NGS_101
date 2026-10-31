# ID-seq — two-stage antibody-reporter barcoding

van Buggenum et al. (2018), *Immuno-detection by sequencing enables large-scale high-dimensional phenotyping in cells*, [10.1038/s41467-018-04761-0](https://doi.org/10.1038/s41467-018-04761-0).

## Inspected originals and scope

Three published originals: PMC6008431 article, Methods **Conjugation of antibodies to dsDNA**, **Immunostainings and release of DNA tags**, **Sample barcoding and sequencing library preparation**, and **ID-seq data analysis**; published MOESM1 Supplementary Note 1 and Figure 1 caption; published MOESM6 Supplementary Data 3, README, production-template/FW/production-RV sheets and first well/plate-primer rows. The manifest identifies the article as published, not an author manuscript. The standard source fetcher and full-bundle oligo scraper were run; relevant sequences were verified in Data 3. Other retrieved originals have not received a chemistry review. This capped pass is not a complete source or missing-information audit.

## Source chemistry

🟢 The assay quantifies proteins in fixed **cell populations**, rather than capturing a transcriptome. Antibody dsDNA tags carry an antibody barcode of 10 nt and a 15 nt UMI. The source conjugation paragraph specifies antibody NHS-s-s-PEG4-tetrazine functionalization, N3-dsDNA with DBCO-PEG12-TCO, conjugation, quenching and purification. Its cited detailed predecessor is van Buggenum et al. (2016), [10.1038/srep22675](https://doi.org/10.1038/srep22675); that original remains deferred. The bounded graph now includes source-reported production/conjugation, staining and DTT release, while modified manufactured termini remain unassembled.

🟢 Immunostained cells are washed, then incubated with fresh 10 mM DTT in BBS pH 8.4 for 90 min at room temperature. Released tags are collected. Collected reporter DNA and residual cell/protein material are distinct physical outputs.

🟢 PCR1 directly introduces a well-specific barcode using PFU, a common forward primer, a well-specific reverse primer and spike-in DNA. Products are pooled per plate, column-purified, treated with Exonuclease I for 30 min at 37°C and heat-inactivated for 20 min at 80°C. PCR2 adds the Illumina index/second adapter using the same common forward and a plate-specific reverse primer. This differs from scID-seq's additional unbarcoded 15-cycle preamplification. The written PCR program contains an initial amplification cycle and a repeat instruction; the model does not silently reinterpret this as a separately sourced cycle count.

🟢 PCR2 products are column-purified and treated with Exonuclease I for 45 min. The **PKIS screen branch** additionally selects fragments greater than 150 bp with Zymo columns, followed by final Qiagen MinElute purification. Digestion products remain in the reaction until physical purification; heat inactivation is not a separation. Two symbolic well classes illustrate pooling without assigning experimental well identities or sample-specific oligos.

## Source oligos and identifier semantics

🟢 Data 3 production template FLD0001:

{{oligo: chem.production_template()}}

🟢 Common forward: `AATGATACGGCGACCACCG`. Production reverse: `TCGCTTATCTGTTGACTGAT`.

🟢 Well primer FLD0151:

{{oligo: chem.well_primer()}}

🟢 Plate reverse primer 25:

{{oligo: chem.plate_primer()}}

The model preserves the workbook's sequence orientation. Its plate-primer index is not presented as a sequencer-specific sample-sheet value. UMI, antibody identity, well sample identity and plate sample identity are structured independently. These are component examples, not an assertion that FLD0001, FLD0151 and plate25 were one assayed well.

🟡 The production **template oligo** ends with the 19-base anchor `ATCAGTCAACAGATAAGCG`; the production reverse primer contains a 20-base reverse segment. The template and its manufactured modified dsDNA therefore require explicit production/end analysis before they can be equated. The complete released reporter and exact PCR junction strands remain unassembled. This is a local scope boundary, not a publication discrepancy or confirmed missing disclosure.

## Output, lineage and deferred work

🟢 Published final libraries are approximately 185 bp and sequenced on NextSeq500; analysis demultiplexes with bcl2fastq and splits around a common anchor to recover UMI, antibody and well identity. Complete molecular product, sequencing-primer binding and run-cycle geometry are not yet modeled. The graph now traces actual loading allocations but still fails its required full-reporter/read-assembly source item; its preview remains private.

The existing scID-seq Transparent Methods explicitly cites this original for downstream barcode sequences. This supports **ID-seq → scID-seq**, component scope: inherited downstream barcode design, with a distinct single-cell preamplification/barcoding order. No predecessor constructor is imported. The 2016 conjugation original is a separate component ancestry candidate, not reviewed recipe inheritance.

Deferred: detailed 2016 conjugation/modified dsDNA production, complete source-specific reporter/end and final amplicon assembly, all remaining Data 3 families and SI graphics, exact sequencing/read layout, remaining linked material and independent source review. No missing-information claim is accepted.

## Bounded producer attempt and retained assembly hold

The 30 October 2026 implementation used only the same three closed originals:
the published main Methods 126–145, SI Note 1 and selected Data 3 components.
No detailed 2016 producer, vendor recipe or additional source was opened.
Source-positive physical stages now include the N3-dsDNA producer interface,
functionalization/conjugation/quench/purification, fixed-cell staining, DTT
release, distinct well contributions, template-stock allocations, two PCR
stages, per-plate pooling, final QA and NextSeq loading. Whole intermediate
structures, individual antibody/well identities and exact modified ends remain
explicitly unassembled. Source PFU-buffer BSA is a component of the same actual
buffer stock, not a separate new dose. Plate25 remains a reference inventory
primer and is not assigned to the experimental plate.

Runtime enumeration of all graph nodes was replaced by a curator-written
literal source target inventory. The actual required `unassembled_source` item
remains unchanged and unaccounted. It specifically includes full reporter/read
assembly; an opaque production interface or final-library label cannot satisfy
that molecular obligation. The normal gate now fails only that item. No
completion marker or public promotion was added.

A source-component assembly attempt shows the printed template's nineteen-base
anchor pairing with the last nineteen bases of the twenty-base production
reverse primer, leaving the primer's first T unpaired. This relation is computed
by a Scene and does not establish the chemically modified manufactured duplex
or the full final PCR/read structure. The N3 position and manufactured ends
remain unresolved without inheriting an unread predecessor recipe. The selected
reference templates/primers are not asserted to be one experimental sample.

The model hold is therefore concrete despite useful physical-producer progress.
It is not a newly observed access failure or a confirmed information absence.
The normal source item must remain until actual matching molecular/read scope
is reconstructed; it must not be renamed, excluded or mapped to an unrelated
opaque factory merely to pass the publication gate.
