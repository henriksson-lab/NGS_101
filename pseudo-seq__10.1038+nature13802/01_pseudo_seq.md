# Pseudo-seq

**Defining source:** Carlile *et al.*, *Nature* (2014), doi:10.1038/nature13802.
The open author manuscript (PMC4224642, NIHMS624208; JATS cached as
`nature13802_PMC4224642.xml`) was read; the publisher version of record is
paywalled (manual row in the source manifest). Nature supplementary information
(`MOESM55_ESM.pdf`, table legends) and source-data workbooks were inventoried.

**2026-10-17 bounded pass** (round breadth-first-15), JATS text lines:
Methods “Pseudo-seq Library Preparation” 56 (RNA isolation, polyA selection,
Zn fragmentation), 58 (CMC/reversal/mock), 60 (CIP/PNK, PAGE fractions,
adapter ligation), 62 (AMV-RT, RT primer, cDNA gel), 64 (CircLigase, Phusion
PCR primers, HiSeq 2000); Results 24/26 and Extended Data Fig. 2 legend
(alternative SuperScript III, fragment sizes, CMC dose); data deposition line 234
(GEO GSE58200). Read length and acquired cycles were not found in these passages.

## Evidence used

- 🟢 Yeast RNA was fragmented in 10 mM ZnCl2 at 94°C (5 min total RNA, 55 s
  polyA+ RNA); HeLa RNA in 10 mM Zn acetate at 60°C for 10 min. Fragmentation
  precedes CMC treatment, which is applied to RNA fragments.
- 🟢 CMC derivatization followed by alkaline reversal leaves Ψ–CMC and a parallel mock
  treatment supplies the background control.
- 🟢 Methods publish CIP/PNK dephosphorylation, 80–100/100–120/120–140-nt
  urea-PAGE fractions, T4 RNA ligase (M0204) ligation in ATP-free buffer
  (22°C, 2.5 h), AMV-RT at 42°C, gel size selection of truncated cDNAs, the modified 3′ adapter,
  RT primer, CircLigase step, both PCR primers and HiSeq 2000 endpoint.
- 🟢 Selected cDNAs terminate at the RT block one nucleotide 3′ of pseudouridine.

## Published library oligos

The open author manuscript's “Pseudo-seq Library Preparation” Methods
(PMCID PMC4224642) prints these oligos. Modifications and the unknown six-base
sample index are retained as written.

| Role | 5′→3′ ordered oligo | Evidence |
|---|---|---|
| 3′ RNA adapter | `/5Phos/TGGAATTCTCGGGTGCCAAGG/3ddC/` | 🟢 Methods; described there as pre-adenylated, although the printed IDT string says `/5Phos/` |
| RT primer | `/5Phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCGGTGGTCGCCGTATCATT/iSp18/CACTCA/iSp18/GCCTTGGCACCCGAGAATTCCA` | 🟢 Methods |
| PCR forward | `AATGATACGGCGACCACCGA` | 🟢 Methods |
| PCR reverse | `CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` | 🟢 Methods; `X` is a variable sample index |

**Published adapter discrepancy.** The same Methods paragraph calls the 3′
adapter “pre-adenylated” but writes its ordered form with `/5Phos/`. A 5′
phosphate and a 5′ adenylate are distinct chemical states. The model uses the
printed bases and modification string and follows the described ligation step;
the source does not establish whether the ordered or reacted adapter carried an
adenylate. 🟡 The ATP-free ligation buffer is consistent with an adenylated
donor; it is not a statement of the ordered modification.

🟡 The paper reports an Illumina HiSeq 2000 run and no separate custom
sequencing-primer oligo. Its PCR primers and adapter establish the library
chemistry, but the current role-level final construct does not yet compute an
exact sequencing-primer binding coordinate.

## Construct logic

Pseudo-seq reads a reaction endpoint rather than the modified base itself. The informative
coordinate is the end of a truncated cDNA, and enrichment over the −CMC control distinguishes
the chemical signal from ordinary RT stops.

## Open questions

- Sequencing read length, read type and sequencing primer are not stated in the
  inspected Methods; GEO GSE58200 run metadata was not inspected (deferred).
- The paywalled Nature version of record was not compared with the author
  manuscript.
- RT primer, adapter and circle junction geometry is still drawn at role level;
  the two iSp18 spacers printed in the RT primer are kept in the
  displayed oligo but not modelled as a molecular feature.

## Local validation migration (2026-10-18)

Re-read published author manuscript PMC4224642 Methods lines 54–64; ran the source
scraper and checked all four modified primer/adapter strings directly. The graph
separates Zn fragmentation, CMC derivatization, alkaline reversal, CIP/PNK repair,
RNA PAGE selection, reported ligation, AMV-RT, truncated-cDNA selection, CircLigase,
Phusion PCR and final gel purification/precipitation. The mock-treated preparation
remains the source comparison; it is not converted into CMC-modified RNA.

The adapter includes the terminal dideoxyC nucleotide. The RT primer preserves both
non-nucleotide `/iSp18/` spacers; its 22-base 3′ region is paired by computation with
the full printed adapter, including that ddC nucleotide. Exact PCR primers and the
six-base sample-index window are explicit inputs. The published phosphate versus
adenylate discrepancy remains visible and its reacted donor state unresolved.

The spacer-containing circle and exact PCR path remain unassembled. PCR uses that
actual circle boundary rather than the expected library as its geometry template.
Read type/cycle counts/sequence-primer landing remain unassigned. These source and
geometry boundaries remain open; no publication-wide absence conclusion is made.
