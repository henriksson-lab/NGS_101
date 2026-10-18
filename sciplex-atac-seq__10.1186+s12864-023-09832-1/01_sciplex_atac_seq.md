# sciPlex-ATAC-seq

**Published version (authority):** Booth et al. (2023) "High-capacity sample multiplexing for single cell chromatin accessibility profiling", *BMC Genomics* 24:737, doi:10.1186/s12864-023-09832-1 (PMC10696879). It supersedes the preprint, bioRxiv doi:10.1101/2023.03.05.531201 (which named the catalogue directory until 2026-10-17; the directory is now named for the published DOI). Read: published JATS (local `s12864-023-09832-1_PMC10696879.xml.txt`) Methods "Hash labeling nuclei" (line 175), "Co-capture … two-level sci-ATAC" (lines 179–197), "Sequencing" (lines 249–263) and hash preprocessing (line 273); published Additional file 2 (`…MOESM2_ESM.xlsx.txt`) Supp_Table4 lines 779–803; and the author [two-level protocol](https://www.protocols.io/view/sciplex-atac-2-level-cz56x89e.pdf) linked by the article (local PDF text lines 94–274). Earlier reviews read the preprint; it is now history, not current evidence.

## Molecular boundary (two-level sciPlex-ATAC)

🟢 Nuclei permeabilized in cold lysis buffer (or OMNI/Pitstop variants) absorb an unmodified sample-specific ssDNA hash, ~0.5 pmol per 1000 cells, 5 min on ice; they are then fixed in 1% formaldehyde for 15 min on ice, pooled across samples, washed and flash-frozen (published line 175).

🟢 Thawed nuclei (5000 per well) receive a well-barcoded capture oligo (25 µM), anneal at 55 °C 5 min, chill, and are extended with NEBNext High-Fidelity 2X PCR Master Mix at 55 °C 10 min. Indexed Tn5 (96 combinations of 8 N5- and 12 N7-loaded complexes, unloaded Tn5 from Diagenode) then tagments at 55 °C 15 min; EDTA/spermidine stops it (line 179).

🟢 Wells are pooled, DAPI-stained and FACS-sorted in limited numbers into 12 µL reverse-cross-linking buffer (EB, SDS, proteinase K); cross-links are reversed at 65 °C for 13.5 h. One PCR with well-specific indexed Nextera P5 and P7 primers and NEBNext HiFi master mix (72 °C 5 min; 98 °C 30 s; 22 cycles) amplifies both hash and ATAC products (lines 179–191; protocols.io step 7). Pooled products are concentrated (Zymo) and separated on a 1% agarose gel: hash library 199 bp, ATAC ~200–3000 bp (line 197).

## Oligos (published)

🟢 Hash: `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` + 10-nt hash ID + `B` + A30 (line 175).

🟢 Two-level capture: `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` + N8 UMI + 10-nt well barcode + `TTT+TTT+…+TTTVN`, where `+` marks the following base as LNA (line 179).

🟢 Indexed Tn5 transferred strands (Supp_Table4 lines 781–800): eight `TCGTCGGCAGCGTCTCCACGC[8]GCGATCGAGGACGGCAGATGTGTATAAGAGACAG` and twelve `GTCTCGTGGGCTCGGCTGTCCCTGTCC[8]CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`; loading strand `/5Phos/CTGTCTCTTATACACATCT` (line 803).

Open: the indexed Nextera P5/P7 PCR primers of the two-level route are named but their sequences were not located in the inspected published Methods or protocols.io PCR step (not a confirmed-absence claim).

## Read layout

🟢 ATAC library: custom primers Read 1 `GCGATCGAGGACGGCAGATGTGTATAAGAGACAG`, Read 2 `CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`, Index 1 `CTGTCTCTTATACACATCTGAGGCGGAGACGGTG`, Index 2 `CTGTCTCTTATACACATCTGCCGTCCTCGATCGC`; Read 1 51, Read 2 51, Index 1 8 + 27 dark + 10, Index 2 10 + 21 dark + 8 cycles, giving T7:P7:P5:T5 (lines 249–259).

🟢 Hash library: standard primers; the text says Read 1 (18 cycles) reads the 8-nt UMI plus the PCR (P5) barcode, Read 2 (≥16) the hash ID, Index 1 (10) the P7 PCR barcode and Index 2 (10) the capture-oligo barcode (line 261).

🟡 **Apparent inconsistency in hash-read assignment:** by the printed capture-oligo sequence, the Nextera Read 1 site (`…AGATGTGTATAAGAGACAG`) is followed by the N8 UMI and the 10-nt capture well barcode, so 18 Read 1 cycles would cover UMI + capture barcode, while Index 2 would read the P5 PCR primer index. The page shows this and assigns no hash-read cycles.

## Published discrepancies shown on the page

🟢 Capture oligo: published Methods line 179 prints LNA `+` marks; the article-linked protocols.io two-level procedure (lines 147–149) prints plain T30VN. The page displays the published LNA form; which formulation each published sample used is unresolved. (The preprint's plain-T form is an ordinary revision difference and is not used.)

🟢 Gel windows: published 199 bp / ~200–3000 bp (line 197) vs protocols.io ~200 bp / 300–1200 bp (lines 229–231). The graph labels the published windows.

## Open questions

🔴 The three-level sciPlex-ATAC3 branch (published lines 199–209: TDE1, PNK, splinted N7/N5 T7-ligase ligations, indexed P5/P7 PCR) is a distinct molecule and not drawn. The sciPlex-RNA hash predecessor and sci-ATAC-seq/indexed-Tn5 predecessors (refs 17, 23, 53) remain unread within this pass.


## Native graph migration (2026-10-18, batch10)

The two-level route now has a native source-owned graph: hash labeling/fixation,
pooling and recovery, freezing, first-well allocation, capture annealing/extension,
indexed transposition, EDTA quench, FACS allocation, reverse crosslinking, initial
extension, 22-cycle PCR, pooling, column recovery and separate gel fractions.
Independent hashes and genomic DNA are never drawn as one covalent polymer.
Original capture molecules carry nine LNA-T residues (T4,7,10,13,16,19,22,25,28);
ordinary PCR daughters do not inherit LNA chemistry. Native complete primer ends,
gap completion and the VN-anchor junction remain explicitly unresolved.

Physical hash/UMI/first-well/PCR-well identifiers remain distinct. ATAC index reads
have eighteen acquired bases each, with intervening dark cycles retained in the
source recipe. Hash Read1/Index2 barcode assignments remain disputed as described
above. The graph preserves source acquisition without inventing exact primer landing
sites on an unassigned final molecule.

This pass read three originals: published Methods169–197 and Sequencing247–271,
published Additionalfile2 Table4 (779–804), and the article-linked two-level protocol
with printed ID87966 (created18Sep2023, modified4Oct2023; capture on pp3–4 and gel
selection step8). The original has no printed DOI. It was not replaced by the
superseded preprint. Linked predecessors, three-level chemistry and sample-specific
PCR primer sequences remain follow-up targets; no publication-absence claim is accepted.
