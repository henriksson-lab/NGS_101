# scGESTALT — research notes

Primary protocol: Raj and Gagnon, *Nature Protocols* (2018), [10.1038/s41596-018-0058-x](https://doi.org/10.1038/s41596-018-0058-x) ([open manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279253/)). Defining study: Raj *et al.*, *Nature Biotechnology* (2018), [10.1038/nbt.4103](https://doi.org/10.1038/nbt.4103).

## 1. What it is

Single-cell GESTALT (scGESTALT) writes cumulative CRISPR–Cas9 indels into an integrated target array at multiple developmental times, expresses the edited array as RNA, and reads it with the same cell identifier as an inDrops transcriptome.

## Sources read

- 🟢 The protocol's experimental design, target-array capture, inDrops amplification, nested lineage PCR and run-configuration sections were read.
- 🟢 Round 2026-10-17-breadth-first-16 re-read PMC6279253 JATS (`.txt` lines 117–121 experimental design; 165–217 reagents; 259–273 primer table; 391–393 primer-mix setup; 741–941 Steps 83–104; 1043 anticipated library size; 1788–1790 Fig. 3b/c legends) and MOESM1 Supplementary Tables 1–2 (lines 60–112).
- The source fetch command for the defining study is `python3 tools/get_sources.py --doi 10.1038/nbt.4103`; the protocol above is the primary chemistry authority for this page.

## Construct-changing steps

1. 🟢 A Tol2-integrated transgene expresses a DsRed-associated CRISPR target array. Early guide injection and later inducible Cas9/guide expression create cumulative indels.
2. 🟢 The edited cassette is transcribed and polyadenylated, allowing it to be captured alongside cellular mRNA in inDrops.
3. 🟢 Reverse-transcribed molecules receive the inDrops cell identifier from the bead primer (T7 promoter–PE1–cell ID–UMI–oligo(dT), Fig. 3b legend), undergo second-strand synthesis and are linearly amplified by T7 IVT. The IVT product is purified with 1.3× AMPure XP (Step 83).
4. 🟢 9 µl of the purified RNA is chemically fragmented for the transcriptome library (Step 84 defers to inDrops ref. 73). Per the Fig. 3b legend, that branch is reverse-transcribed with a random-hexamer primer carrying PE2 and then PCR-amplified with P5/P7 primers; the primer table lists R1-PCRix/R2-PCR for "transcriptome or scGESTALT final sequencing library". There is no RNA-adapter ligation (a repository error corrected in round 16). The remaining RNA stays unfragmented for the lineage library.
5. 🟢 Lineage branch: random-hexamer RT with PrimeScript (70 °C 3 min; 30 °C 10 min, 42 °C 1 h, 70 °C 15 min; 1.2× AMPure). PCR round 1 uses GP6 + PE1Sa with Q5 and Q5 enhancer for 15 cycles at 61 °C (0.6× AMPure). PCR round 2 uses GP12 + PE1Sb with Phusion at 60 °C (0.6× AMPure). A final KAPA HiFi HotStart PCR on a 1:10 dilution uses R1-PCRix/R2-PCR: 2 cycles at 55 °C, then typically 6–8 at 65 °C (0.6× AMPure). The same R1-PCRix/R2-PCR pair as the matching transcriptome library must be used. Libraries are typically 500–800 bp.
6. 🟡 Published inconsistencies, shown on the page: (a) the primer table cites Supplementary Table 2 for the R1-PCRix indexes, but the primer-mix setup cites Supplementary Table 1; MOESM1 Table 1 lists sgRNA sites and Table 2 lists R1-PCR-ix1–24. (b) The round-2 thermocycler table numbers its rows 1, 2–9, 16, 17, so the cycle count is either 8 or 14. The model states no count.

## Identifier structure and read layout

- 🟢 Transcriptome libraries: NextSeq 75-cycle high output, standard primers; Read 1 61, Index 1 8 (cell-barcode part 1), Index 2 8 (library index), Read 2 14 cycles.
- 🟢 The lineage library run (MiSeq 300-cycle kit, 20% PhiX, standard Illumina primers) is: Read 1, 250 cycles (edited barcode array); Index 1, 8 cycles (cell-barcode part 1); Index 2, 8 cycles (library index); Read 2, 14 cycles (8-base cell-barcode part 2 followed by 6-base UMI).
- 🟢 The two cell-barcode parts form one combinatorial cell identifier. The UMI identifies lineage-transcript molecules; the library index is a sample/library index.
- 🟡 The upstream inDrops capture oligos and a fully assembled final-library strand remain role-labelled. The scGESTALT-specific amplification primers are transcribed below.

## Boundary

The public page documents scGESTALT's expressed lineage-cassette branch. GESTALT genomic-DNA recovery and later optimized scGESTALT cassettes are related protocols, not silently merged chemistry.

## Primer inventory from the primary protocol

The primary [Raj and Gagnon protocol](https://pmc.ncbi.nlm.nih.gov/articles/PMC6279253/) prints a reagent table of primers for the lineage branch. These sequences are 🟢 as printed; spaces within the long GP12 entry are merely line grouping in the source, and `N` represents the nine random bases.

| Oligo | Sequence, 5′→3′ | Role |
|---|---|---|
| inDrops_GP6 | `GAGGACTACACCATCGTGGAG` | First, long lineage PCR |
| inDrops_PE1Sa | `CTCTTTCCCTACACGACGCTGGGTGTCGGGTGCAG` | First, long lineage PCR |
| inDrops_GP12 | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNTCGAGCTCAAGCTTCGGAC` | Second, short lineage PCR; random bases add sequence diversity |
| inDrops_PE1Sb | `CTCTTTCCCTACACGACGCT` | Second, short lineage PCR |
| R1-PCRix, first indexed member | `AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC` | Final library PCR; Supplementary Table 2 prints this complete indexed oligo |
| R2-PCR | `CAAGCAGAAGACGGCATACGAGATGGGTGTCGGGTGCAG` | Final library PCR |

The same table also prints editing-assessment primers scGESTALT F `TCGAGCTCAAGCTTCGG` and scGESTALT R `CTGCCATTTGTCTCGAGGTC`. Those amplify the genomic cassette for QC and do not define the sequenced lineage-library ends. Supplementary Table 2 gives additional indexed members; the first member above is a concrete example, not a universal fixed index. The inDrops V3 capture primer/barcode oligos belong to the upstream inDrops protocol and must be reconciled before calling the full library inventory complete.


## Native graph migration — 2026-10-18 batch11

Source scope (2 originals, capped): PublishedPMC6279253 ReagentsprimerTable263–273, Procedure82–105 approximately731–941 andFig3b/c; publishedMOESM1 originalSupplementaryTable2ix1–24, includingix1–6 approximately91–116.

Source-owned graph replaces rna_cell_assays defaultendpoints/states. ActualfirstSS, sourcepositive secondstrand/T7promoterduplex, IVTRNA and original/copyfates are distinguished. Both9µl expression and5µl unfragmentedlineage allocations retain the full transcriptmixture ratherthan modalitypurification; unallocatedRNA isunresolved. SourcePrimeScript randomhexamerSS feeds actualnestedPCRinputs, followedsourcecleanups/finalindexPCR. GP12N9 remains diversity, notUMI; multipartcellID andN6UMI persist withsource8ntlibraryindex identity. Bothknown61/8/8/14 expression and250/8/8/14 lineage runs remain. Publicprimer-tablecrossreference andround2cycle-row discrepancies are preserved.

Deferred scope: inDrops citedsource73 nativecapture/bead/RT/secondstrand/expressionlibrary components not read thispass; reactionorder known fromdefiningFig3/Procedure isretained without exactnativegeometry. Lineageeditedarray sequence/extent is source/sampledependent; no generic44ntinsert default. Round2amplificationcount remains sourceambiguous; final55°Ctwo+65°Csix-to-eightcycles sourcequalified.

Verification: 5 originalprimer tests;19steps/31unresolved0contradictory afteroppositestrand/actual5µlsampling correction. Allfourpages build; sourceguard17checked/0failed/9unassigned; pairing0; review schema checked. No sourceabsence/status upgrade follows.
