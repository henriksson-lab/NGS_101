# Quantum Barcoding (QBC) — source notes

**Defining source:** Roy et al., *Communications Biology* (2020),
[doi:10.1038/s42003-020-0896-2](https://doi.org/10.1038/s42003-020-0896-2)
(PMC7205613, CC BY). An Author Correction exists
(doi:10.1038/s42003-020-1019-9, PMC7260357); it was not inspected.

**2026-10-17 bounded pass** (round breadth-first-15). Read: PMC7205613 JATS
(`s42003-020-0896-2_PMC7205613.xml.txt`) Methods lines 216–258 (antibody
conjugation, cell preparation, in-cell RT, aminoallyl crosslinking, staining,
split-pool, post-split-pool cDNA preparation, PCR, primer design), 286–292
(SNAIL post-split-pool, sequencing); Supplementary Data workbook (MOESM3)
sheet 1 (QBC alignments and generic ordering forms) and sheet 2 (60 AHCA
sequences); MOESM2 (supplementary item list). Sheets 3–13, the SNAIL and
TotalSeq-B branches, the Supplementary Information PDF and the NEBNext vendor
guides were not read.

## What the sources establish

- 🟢 Fixed (1.6% formaldehyde), methanol-permeabilized cells are proteinase-K
  treated and reverse-transcribed in cell with ProtoScript II, target-specific RT
  primers (Supplementary Data 5–7) and aminoallyl-dUTP; DTSSP crosslinks the cDNA
  (Methods lines 224–230).
- 🟢 RT primers carry the splint-annealing anchor `5′GCTCCCTGTCTGACG` on their
  5′ end, terminating in a 5′ phosphate (Methods, “Primer design”).
- 🟢 Antibodies are conjugated to La4FB
  `5′-(4FB)-AAGAGCTAGTTATTGCTCAGCGGAATAAAGCTGATGGAGTTCGTGACTGG(Propyl)-3′`;
  each conjugate is pre-hybridized to its own 5′-phosphorylated AHCA oligo,
  stained, then refixed in formaldehyde (Methods lines 218, 234).
- 🟢 Splint6 is hybridized once (37°C, 1 h); four rounds then split cells into 60
  wells, anneal and ligate a subcode (T4 DNA ligase M0202, 0.3 mM ATP, PEG buffer,
  15 min RT), wash and pool (Methods lines 238–240). Subcode lists are
  Supplementary Data 3 (SC1, SCN, SC2) and 4 (SC3).
- 🟢 Supplementary Data 1 ordering forms: SC1 `p-GCGACTTCC xxxxxxx CCAACCACC-OH`,
  SCN `p-GTCTTTGCC xxxxxxx CCGTTGTGG-OH`, SC2 `p-GCGTTACCC xxxxxxx CCCATTCTC-OH`,
  SC3 `CCCTACACGACGCTCTTCCGATCT nnnnnnnnnn xxxxxxx GCTCACCAAGTCCTAGT-OH`
  (no 5′ phosphate). Each incoming subcode's 3′-OH is joined to the 5′ phosphate
  of the growing reporter, so SC3 ends the chain with the partial Read 1 handle
  (🟡 derived from the printed ends and the alignment).
- 🟢 Splint6 bases and spacers (Methods line 238):
  `CGTCAGACAGGGAGCGGTGGTTGG/iSpC3/GGAAGTCGCCCACAACGG/iSpC3/GGCAAAGACGAGAATGGG/iSpC3/GGGTAACGCACTAGGACTTGGTGAGC/3SpC3/`.
  Its first 24 nt dock the anchor and SC1 right arm; the last 26 nt dock SC2 and SC3.
- 🟢 The antibody amplicon alignment (Supplementary Data 1) is 204 nt; the
  sequenced block (UMI through AHCA code) is 133 nt.
- 🟢 Targeted-RNA samples: RNase H, proteinase K lysis, then second-strand
  synthesis with RNA-specific primers (Supplementary Data 8–10) carrying
  `5′GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, Herculase II, 1–20 cycles; some
  mixtures also got thermolabile Exo I to remove unligated subcodes; then
  NEBNext E7600S indexed PCR (12–30 cycles). Protein-only samples go straight to
  PCR with NEBNext E7335 primers and PfuUltra II (24–28 cycles).
- 🟢 Runs: MiSeq, HiSeq and NextSeq, 151–258-base single- or paired-end; paired
  reads were overlapped with FLASH.

## Model corrections in this pass

- The RNA amplicon now ends in the reverse complement of the second-strand
  primer's Illumina arm; it previously inherited the AHCA-only `GAACTCCATCAGC` tail.
- A second-strand synthesis step joins the released cDNA to the RNA PCR template.
- Subcode inputs: SC3 is no longer called phosphorylated. The ligation geometry
  now uses each subcode's 3′ arm as the splinted ligating end, not the 7-nt code.
- Splint6 scene labels follow the docking arms, and it now shows the 3′ C3 spacer.

## Published discrepancies (shown on the page)

- AHCA anchor: Methods line 234 prints `5′-(phos)-CTCCCTGTCTGACG(x9)…`, while
  Supplementary Data 1/2 and “Primer design” give `GCTCCCTGTCTGACG`. The model
  uses the 15-nt form, which the Splint6 5′ arm exactly complements (🟡 derived).
- Splint6 5′ end: Methods prints no 5′ modification; Supplementary Data 1 prints
  `p-`. Bases and spacers agree. Not resolved.

## Open questions

- Not inspected: vendor index layout for NEBNext E7335 versus E7600S. Both libraries
  are drawn with placeholder dual 8-nt indices; acquired index cycles are not
  stated in the inspected passages.
- Supplementary Data 3–11 (subcode lists, RT/second-strand primers, SNAIL probes)
  were not transcribed.
- The TotalSeq-B linker/Phi29 extension and SNAIL in-situ amplification variants
  are not drawn.
- The Author Correction has not been checked for chemistry content.
