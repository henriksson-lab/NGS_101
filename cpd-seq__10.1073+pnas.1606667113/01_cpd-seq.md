# CPD-seq

Defining DOI: https://doi.org/10.1073/pnas.1606667113

Mao et al. (2016), Chromosomal landscape of UV damage formation and repair at single-nucleotide resolution.

## Source scope

🟢 Published PMC4987812 article page, Results/Figure1 and inline FigureS1 (text1120–1143), Materials and Methods1266–1268, inline SI Materials and Methods1277–1297: genomic preparation, Oligonucleotides, Validation, Library Preparation and Sequencing, Data Analysis. The inline SI is part of the inspected PMC page export; no equivalence to an uninspected separate PDF is asserted.

## Preparation and physical readout

1. 🟢 Yeast genomic DNA after source UV treatment/repair incubation. Sonicate near400bp; ethanol precipitation and1.2×AMPure selection; end repair and dA tailing with NEBNext modules. Ligate source trP1 T-overhang duplex16°C16h and AMPure-purify.
2. 🟢 TdT plus ddATP37°C2h blocks pre-existing free3′ ends; AMPure cleanup follows. This is not a prediction of complete blocking: the source discusses residual free3′OH background.
3. 🟢 T4endoV37°C2h creates single-strand breaks at CPD sites, physical purification, then APE1 digestion37°C1.5h and purification. The source identifies the resulting free3′OH immediately upstream of the lesion (FigureS1). The exact end between the two enzymes is not assigned from an enzyme name alone. Complementary strands, remaining lesion-bearing DNA and other cleavage products persist as molecular species until physical separation.
4. 🟢 Shrimp alkaline phosphatase37°C1h removes5′ phosphate; heat65°C15min and purify. Denature95°C5min, cool on ice, ligate source A adapter16°C16h at the newly available3′OH. The adapter's5′ phosphate is distinct from the dephosphorylated genomic5′ end.
5. 🟢 Biotin on the A bottom oligo captures adapter-containing molecules on streptavidin beads.0.15MNaOH releases the NONbiotinylated DNA strand for library copying. The biotin strand remains bead-associated; it must not be drawn as covalently retained on the released genomic template. Capture and alkali-release fates are explicit and do not assume quantitative purity.
6. 🟢 PrimerA copying gives dsDNA; AMPure cleanup, diagnosticPCR, then six-cycle libraryPCR with primerA/trP1. AMPure-purify, pool equal volumes of three differently barcoded libraries and sequence on IonProton. The selected local preview uses the A1 reagent branch; A2/A3 are alternative sample codes, not extra protocols.
7. 🟢 Demultiplex/trim barcode then align. Infer a candidate CPD dinucleotide adjacent to the read5′ end on the OPPOSITE strand. Excluding nondipyrimidine or pUC19 reads is computational filtering after sequencing, not molecular disposal. APE1-associated one-base offsets/background are source limitations, not silently repaired genomic positions.

## Oligos and model boundary

Source trP1-top: 5′-CCTCTCTATGGGCAGTCGGTGAT-phosphorothioate-T-3′; trP1-bottom: 5′-phosphate-ATCACCGACTGCCCATAGAGAGGC-dideoxy-3′. Source A1-top: 5′-phosphate-ATCCTCTTCTGAGTCGGAGACACGCAGGGATGAGATGGC-dideoxy-3′; A1-bottom: 5′-biotin-CCATCTCATCCCTGCGTGTCTCCGACTCAGAAGAGGATNNNNNN-C3 phosphoramidite-3′. The N6 end is a degenerate ligation/splint window; no UMI meaning is asserted. PrimerA: CCATCTCATCCCTGCGTGTCTCCGAC; primertrP1: CCTCTCTATGGGCAGTCGGTGATT.

The private model owns these exact source oligos and source-described end histories, but does not claim complete final strand assembly, run-primer identity/read length or barcode cycle boundaries. Historical kit/end chemistry remains unassembled. Any text-parser limitation on the written “phosphorothioate-T” is documented with a narrow source test, rather than changing the oligo. No finished publication marker or complete chemical-validation claim is emitted.

## Lineage and deferred work

Main Methods1268 explicitly adopts emRiboSeq with altered cleavage chemistry; ref19 resolves at1387 to Ding et al., NatProtoc2015 DOI10.1038/nprot.2015.099. That is a companion protocol, not the registered emRiboSeq Nature14183 defining original. Component ancestry must preserve this cited version; source-local oligos were read independently and not imported. Direct enzyme reference35, historical Ion sequencing documentation, complete duplex/junction validation and independent review remain deferred.

## Bounded publication follow-up (2026-10-30)

Main Methods1268 explicitly adapts CPD-seq from the emRiboSeq/EndoSeq protocol cited as reference19 (text1387), DOI10.1038/nprot.2015.099. That exact relevant original remains unreviewed. The required `companion_protocol` item is preserved and uncovered: the positive inline-SI summary1295 does not waive it or establish complete historical recipe equivalence. Source-known oligos/end histories and coarse DNA-copying/PCR steps remain private; full strand, primer and run assembly is unimplemented.

The limited local fix exposes dideoxyATP as an actual molecular substrate of the TdT end-blocking step, separately typed from the catalytic enzyme. Dose, packaging and quantitative blocking remain unassigned; no new fresh-stock recipe is inferred. All other source-scope obligations, background/offset qualifications and nonbiotin-release fates remain intact. This pass reread only the existing published main/inline-SI original; no acquisition or additional original. Notes status remains a source-scope hold.
