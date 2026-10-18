# Paired-seq — research notes

## Sources read

- 🟢 Zhu et al., *Paired-seq: simultaneous profiling of transcriptome and accessible chromatin in millions of single cells*, *Nature Structural & Molecular Biology* (2019), DOI [10.1038/s41594-019-0323-x](https://doi.org/10.1038/s41594-019-0323-x); open manuscript [PMC7231560](https://pmc.ncbi.nlm.nih.gov/articles/PMC7231560/).
- 🟢 Supplementary Table 1 gives pMENTs (`5Phos/CTGTCTCTTATACACATCTddC`), AdaptorA, Linker-R02/R03/R04, Blocker-R02/R03, Terminator, Anchor, PA-F/PA-R and the indexed PCR primers.
- 🟢 Supplementary Table 2 gives all eight modality-specific BC1 oligos and the 96-member BC2, BC3 and BC4 plates, including round-4 phase blocks.
- Copies read (2026-10-17 bounded pass): author-manuscript JATS `s41594-019-0323-x_PMC7231560.xml.txt` Results line 65 and Methods lines 117–149 and 159 (transposome and barcode-plate assembly, tagmentation/RT, three ligation rounds with blockers/terminator, lysis, TdT tailing/Anchor extension/preamplification, restriction split and Tn5-P5, indexing PCR and run, Read-2 barcode positions); Supplementary Tables 1 and 2 as the publisher MOESM2/MOESM3 workbooks (the bundle holds them under `scg_41594_2019_323_MOESM2/3_ESM.xlsx`; the manifest records the Springer copies as identical), Table 1 lines 1–15, Table 2 lines 1–20, 117–119 and 214–252.

## Construct-changing steps

1. 🟢 Barcoded Tn5 transfers `common arm—BC1—NotI—ME` to accessible DNA. In parallel, the RT primer installs `common arm—BC1—SbfI—T16VN` on RNA-derived cDNA.
2. 🟢 Linker-R02 and Linker-R03 splint the published phosphorylated BC2 and BC3 oligos to the preceding 5′ arm. Blocker-R02 and Blocker-R03 are the exact reverse complements of their linkers.
3. 🟢 Linker-R04 splints the terminal oligo, which carries the TruSeq Read-2 arm, 10-nt UMI, BC4, a 0–3-nt phase block and its ligating arm. The Terminator is the reverse complement of Linker-R04.
   After 30 min of each of rounds 2 and 3, Blocker-R02/Blocker-R03 are added and ligation continues for 30 min; round 4 is quenched with the R04 Terminator plus EDTA (Methods line 133). The page shows them as round inputs.
4. 🟢 After lysis (NaCl, SDS, proteinase K, 55 °C 2 h; 1× SPRI), TdT adds poly(C); the nine-G tract in Anchor primes linear extension. PA-F and PA-R then preamplify the indexed molecules.
5. 🟢 The material is divided. NotI cleaves the DNA architecture in the RNA aliquot, while SbfI cleaves the RNA architecture in the DNA aliquot; cleavage does not prove physical removal. Tn5-P5 retagmentation and N5XX/P7XX PCR complete both sequencing libraries.

## Read structure

🟢 The defining run was paired-end `53 + 7 + 130` (Read 1 + i7 + Read 2). Read 1 is payload. Read 2 begins with the UMI and traverses BC4, BC3, BC2 and BC1; the phase blocks shift their first cycles across the published ranges 10–13, 47–50, 84–87 and 121–124. The page draws a published zero-phase member, for which those starts are cycles 11, 48, 85 and 124.

🟡 The run acquired 7 Index-1 cycles, but P7XX prints a 6-nt `NNNNNN` i7 window; the seventh cycle would read into the adjacent P7XX adapter base. This is a derivation, not a published conflict.

Model corrections (2026-10-17): the RNA library now keeps the poly(A) copied from the oligo-dT/VN primer between cDNA and SbfI; the base before the TruSeq Read-2 arm was labelled a "dA junction" but is the last base of the R04 Read-2 arm (`...CGATCT`), so the model now draws the R04 arm and the P7XX extension separately (sequence unchanged).

The N5XX oligo contains an 8-nt i5 sequence, but the stated defining run does not assign index2acquisition. The completed construct retains that sequence while the read table retains the three positive acquired reads and explicitly leaves index2acquisition/cycles unassigned.

## Open questions

- Tn5-P5 loading: Methods does not state which adapter is loaded; AdaptorA in Table 1 is the likely candidate but is not assigned.
- TELP predecessor (ref. 43) for the TdT/Anchor chemistry was not followed.
- Cleanupsteps now have separategraphproducts; quantitative recovery ofrestriction-cleaved opposite-modality fragments remains unassigned.
- Extended Data Fig. 1a (read layout drawing) was not inspected; the JATS text was used.

## Native graph migration — 2026-10-18 batch03

Three published originals were inspected: PMC7231560Methods and SupplementaryTables1/2. The sourceworkbooks happen to have scg_ cache names; MANIFEST confirms the corresponding publisheroriginals are md5-identical. Their identity is the published tables, not the upstream drawnpage. Local source-owned splints preserve BC1/2/3/4 features and Round4N10UMI. SourcepMENTs writes a physical terminalddC; the nontransferredstrand retains thatC but the transferredtaggedDNA does not copy it. A scoped sourceguard parserexception for thisinput is covered by an exact wholemodifiedpMENTs regression plus a sourcewhole transferredstrand check; scraper._letters otherwise drops the nucleotide as a modification. The Anchor9G tract does not establish a nine-C TdTtail. Native tail/register/extension remainsunassigned and actualAnchor-extension products feed PA-PCR. Restriction cuts damage the opposite-modalityarchitecture inside each mixed aliquot; they do not physically discard thosemolecules. Sourcecleanup,QCaliquot andsublibrary handling are separategraphstates. The exactprintedP7XX ends...CTTCCGATC, onebase shorter than canonicalRead2; priorcode silently appendedT. The correctedprimer matchesTable1, while the R04terminalT/fullRead2 endpoint junction remains a conditional polymeraseextension reference. That is a repositorytranscription correction, not a publishedconflict. SourceTn5-P5 is named; assigning Table1AdaptorA as its exactload remainsunresolved. Positive53+7+130acquisition and6nti7/8nti5 physicalprimerwindows remain, whileIndex2acquisition/cycles are explicitunassigned.

Curator-written local review targets identify the source/implementation scope; targets are not derived from graph contents. This capped pass does not establish publication-wide absence or complete missing-information review.

### Implementation review corrections

Implementation corrections align displayed P7XX with its printed32nt arm, remove the obsolete workflow, preserve source barcode/UMI chain features in Anchor/preamp and intact restriction-target states, and mark a new native P5-associated end/gap role after retagmentation while retaining carried cleavage products. Native tail length and full first-cycle/end assembly stay unassigned. These corrections reference the original source pass above; no additional originals were read.
