# CFMS mRNA display: single-round selection and NGS

Olson et al. (2012), **Single-Round, Multiplexed Antibody Mimetic Design through mRNA Display**. DOI: [10.1002/anie.201207005](https://doi.org/10.1002/anie.201207005).

This is a private incomplete source graph. The author-institution PDF `source_pdf_3f348b50ea6f42f28215.pdf` actually combines the printed final main article (Angew.Chem.Int.Ed.51,12449–12453) and published Supporting Information. Its source-index label and earlier discovery preview exposed only the SI; actual reading establishes both components. The associated author-manuscript title uses “via”; it is not used as a coequal recipe or claimed completely equivalent.

## Selected originals and oligos

🟢 Combined published main: title/DOI, scaffold and Figure1. 🟢 PublishedSI DetailedMethods950–1054 and TableS1, PDFp21/SIp15, rendered original independently inspected. Mandatory scraper10 candidates; its merged table candidate is not one660nt manufactured oligo. Table row layout establishes the separate primers.

Five source primer bodies are represented locally; all are plain DNA as printed in this selected table:

| Name | Source5′→3′ sequence |
|---|---|
| FnRTDS | AATGATACGGCGACCACCGAGATCTACACTGTGCGGTAGTTGATGGAGATCG |
| OligoC | AATGATACGGCGACCACCGAGATC |
| OligoD | GCAAGCAGAAGACGGCATACGAGA |
| FNBCSEQ | CCAGTTTCACCGTAGGTGATGCGGTAGTAGCGA |
| FNFGSEQ | TGTCGACTATACCATCACGGTGTACGCCGTCAC |

Canonical matching Illumina portions are imported; complete primer bodies are checked against the actual TableS1 text and original layout. The source pF30P puromycin linker and splint remain unassembled: a name is not a complete modified-oligo sequence. BC/FG randomized protein-coding loops identify variants, not molecules/UMIs.

## Source operational route and physical fractions

1. The BC cassette undergoes source staged KOD extension/amplification, while the FG cassette is independently extended/amplified. Both are spin-column purified, BsaI digested and purified again; T4DNA ligase combines them. Exact cohesive ends/junctions are not inferred merely from BsaI. Source2%agarose extraction/spin purification and a concentration assay precede dilution to the selected2.5e9-molecule aliquot. Remaining library allocations and an input-DNA reference are separate.
2. Source19cyclePCR scales the template, then phenol/chloroform extraction/ethanol precipitation precede T7 runoff transcription. RNA receives its own extraction/precipitation. RNA, pF30P and splint undergo T4DNA ligase treatment, denaturing4.5%ureaPAGE, electroelution and ethanol precipitation. Only then is800pmol allocated for two target selections; remaining preparation is not forced into this experiment.
3. Rabbit reticulocyte lysate translation30°C1h is followed by KCl500mM/MgCl2 60mM and30min room temperature to promote covalent puromycin-mediated fusions. Unreacted classes are not assigned ideal removal.
4. OligodT-cellulose capture and washes are followed by three500µL low-salt eluates, collected separately from matrix. SuperScriptII/FnRTDS creates encoding cDNA while original RNA and the peptide linkage remain. EDTA follows RT. Anti-FLAG/M2 capture and washes select the intended fusion population without asserting perfect removal of nonfunctional species; FLAG-peptide elution and physical collection follow.
5. BSA is added and1%purified fusion pool reserved for composition sequencing. The remaining population undergoes four neutravidin-bead preclears: retained **unbound** fractions continue and beads are physically removed. Equal fractions enter the two source target-bead selections; the model follows one representative target while keeping the other separate.
6. After target binding, CFMS withdraws30mL/h through magnet-bearing tubing, washes3.5min and evacuates buffer. Trapped beads remain until magnets are removed and100µLPCR mixture resuspends the beads.95°C2min dissociates cDNA; no invented bead-removal step occurs before source17cycle KOD PCR with OligosC/D. New DNA copies do not become covalent peptide–RNA fusions.

🟢 Amplified inputDNA, inputfusioncDNA and selected material were sequenced on GAIIx,6pM, paired-end42 cycles with FNBCSEQ then FNFGSEQ. Published expected read patterns extend beyond42 bases for the first pattern; the model records PE42 without silently changing the source's printed pattern or inferring a new run. Full postPCR loading workup and all reference-library amplification assignments remain unassembled.

## Ancestry and pending scope

The publication explicitly incorporates prior e10Fn3 scaffold/fusion and CFMS components (Olson2011; DetailedMethods/refs). These are limited scaffold and affinity-platform inheritances; no complete source recipe is imported from a generic mRNA-display family. The directly cited original, pF30P/splint and complete assembly geometry remain pending. No accepted absence claim or public promotion follows this capped pass.
