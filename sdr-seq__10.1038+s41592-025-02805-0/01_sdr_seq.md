# SDR-seq

**Defining paper:** Lindenhofer et al., *Nature Methods* (2025), doi:10.1038/s41592-025-02805-0. The open article (PMC12510883 JATS) and detailed methods were read; the associated detailed protocol is doi:10.17504/protocols.io.6qpvr9q43vmk/v1.

## Molecular boundary

🟢 Fixed permeabilized cells undergo in-situ RT with custom poly(dT) primers that add a UMI, sample barcode and capture sequence to cDNA. gDNA remains in the same cell.

🟢 In Tapestri droplets, multiplex PCR copies a bead cell barcode onto targeted gDNA and cDNA amplicons through complementary capture-sequence overhangs.

🟢 Version 2 reverse primers distinguish the branches: gDNA primers carry Nextera R2N and RNA primers carry TruSeq R2, enabling separate sequencing-ready libraries.

## Evidence boundary

🟢 Target panels and sample-barcode oligos are published in supplementary tables. They vary with the chosen loci and samples; the page models the invariant UMI/barcode/overhang architecture.

## Oligo source audit

🟢 Supplementary Table 1 (proof-of-concept v1) and Supplementary Table 2 (v2) are in `s41592-025-02805-0_supp_41592_2025_2805_MOESM3_ESM.xlsx`, sheets `Supplementary_Table_1` and `Supplementary_Table_2`. These tables give actual primer sequences, not only an assay diagram. The invariant v2 droplet-amplification forward primer is `GTACTCGCAGTAGTCGACACGTC`; the v2 RNA reverse primer is `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT-cDNA_specfic` as printed, where the locus-specific suffix is a placeholder. The corresponding v2 gDNA primers are `GTACTCGCAGTAGTC-gDNA_specific` and `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-gDNA_specific`.

🟢 The same table prints complete v2 RT-primer examples, including `RT_primer_1_8` as `GTACTCGCAGTAGTCGACACGTCTCGCCTTANNNNNNNNTTTTTTTTTTTTTTTTVN`. Its sample-barcode choice changes across the table. Separate i5 primers and branch-specific i7 primers are listed there as full sequences; these are choices from the published index sets, not a single immutable library end.

🟡 The diagram's role-level representation remains appropriate for the variable target panel and barcode/index choices. The source does resolve the common primer handles above; any exact read-placement model should select a stated v1/v2 and sample-index instance rather than treating the variable rows as one molecule.

## Bounded pass 2026-10-17 (round 17)

🟢 Methods "SDR-seq in human iPS cells" (PMC12510883 JATS lines 207–225): glyoxal or PFA fixation/permeabilization; in-situ RT with Maxima H Minus RT and the Supplementary Table 1/2 RT oligos (48 wells × 10,000 cells); Tapestri v2 device; RNA primers from TAP-seq design (150–300-bp products), gDNA primers from Tapestri Designer; V1 gDNA and RNA reverse primers both carry R2N (POP only), V2 gDNA R2N and RNA R2; V1 libraries per the Mission Bio user guide, V2 libraries built separately with Supplementary Table 2 primers. Results "Droplet-based SDR-seq" (line 109): first droplet = lysis, proteinase K and reverse primers; second droplet = CS-overhang forward primers, PCR reagents and barcoding bead.

🟢 Supplementary Table 2 (MOESM3 sheet `Supplementary_Table_2`): the V2 RNA forward primer `Droplet_amplification_fwd` (`GTACTCGCAGTAGTCGACACGTC`) equals the 5′ handle of every V2 RT primer, so the RNA branch uses one universal forward primer on the RT handle while gDNA targets use `GTACTCGCAGTAGTC`-tailed specific forward primers. i5 primers `Nextera_R1N_S50x` are shared by both libraries.

🟢 Methods "SDR-seq data analysis" (line 237): read barcode structure is cell BC1, variable 14–17-bp linker, cell BC2, 15-bp constant linker, then sample BC and UMI (RNA only). 🟡 That the 15-bp constant linker is the 15-nt capture sequence is inferred.

Model corrections: RT-handle order fixed to CS handle → sample BC → UMI → T16VN → cDNA (previously UMI/sample BC were placed after the cDNA target or before the CS); cell barcode split into BC1/linker/BC2/constant linker with part features; proteinase K lysis step in first droplets; universal RNA forward primer as an exact input; emulsion-breaking cleanup; named i5/i7 primer families; Maxima H Minus RT.

🔴/discrepancy: Supplementary Table 3 lists the V1 POP RT primers as `RT_primer_1_8`–`RT_primer_4_8` (V2 names, with the `GACACGTC` spacer), whereas Supplementary Table 1 lists V1 `RT_primer_1`–`RT_primer_4` without it. Shown on the page; unresolved.

Not found in the inspected sources: sequencing instrument, read lengths and index cycles (main text and Peer Review File grep); Mission Bio bead oligo bases. Deferred: protocols.io 10.17504/protocols.io.6qpvr9q43vmk/v1 (not on disk), Mission Bio Tapestri v2/v3 user guides, SDRranger repository (github.com/hawkjo/SDRranger) for barcode backbone sequences.
