# sci-Plex — source notes

Defining paper: Srivatsan SR, McFaline-Figueroa JL, Ramani V, et al. "Massively
multiplex chemical transcriptomics at single-cell resolution." *Science* 367:45–51
(2020), [doi:10.1126/science.aax6234](https://doi.org/10.1126/science.aax6234);
author manuscript PMC7289078.

> **Evidence marking.** 🟢 read in an original source · 🟡 derived or inferred ·
> 🔴 not read (source unavailable to this review; not a claim of publication absence).

## Sources read

| Original | Where | Used for |
|---|---|---|
| Author manuscript JATS (`science.aax6234_PMC7289078.xml.txt`) | Introduction line 77, Results lines 83–89, data statement line 159, Fig. 1 legend line 432, Fig. 3 legend line 444 | workflow, sci-RNA-seq versions, transposase provenance |
| Author IDT order forms (`gh_plate1…12_for_ordering.xlsx.txt`, github.com/cole-trapnell-lab/sci-plex `reagents/IDT_hash_oligo_order_forms`) | all 1,152 rows (plate 1 well A1 transcribed) | hash oligo sequence |

🔴 The Science Supplementary Materials (Materials and Methods, Tables S1–S8;
PMC NIHMS1595611 supplements and `science.sciencemag.org/content/367/6473/45/suppl/DC1`)
remain gated (PMC download gate; see `MANIFEST.tsv`). They hold the RT, ligation and
PCR oligos, enzymes, conditions and sequencing run. GEO GSE139944 (named in the data
statement) has not been fetched.

## What the inspected sources establish

- 🟢 Cells are lysed in the presence of well-specific polyadenylated ssDNA "hash
  oligos" (unmodified ssDNA) and the nuclear suspensions are fixed with
  paraformaldehyde; the hashed nuclei are then pooled (Results line 85; Fig. 1B).
- 🟢 The barnyard pilot used a **two-level sci-RNA-seq** experiment; the 188-compound
  screen used **sci-RNA-seq3** (Results line 85; Fig. 3A legend). Because the hashes
  are polyadenylated they are combinatorially indexed like endogenous mRNA.
- 🟢 Transposase complexes are as described in Cao et al. (ref. 25), supplied by
  Illumina on request (data statement).
- 🟢 Every one of the 1,152 order-form rows has the form
  `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` + 10-base hash + literal `B` + 32 A
  (checked with a pattern match across all twelve plates). Plate 1, well A1:

```
sciPlex_0  5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTTCTCGCATGBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA-3′   (IDT, no modifications printed)
```

  The column header calls the arm "Nextera Read2". `B` is kept as printed; its
  IUPAC reading (C/G/T) is not stated by the form.

## Derived relationships

- 🟡 *(computed)* The hash arm equals `nextera.READ2_PRIMER` (s7 + mosaic end); the
  module raises if the transcription drifts. The oligo-dT-primed copy of a hash
  therefore ends in the reverse complement of that arm, so the hash-derived molecule
  carries the Nextera Read 2 site without tagmentation, and Read 2 reads the hash
  barcode as printed.
- 🟡 Transcript-derived molecules are sequenced in the same pool; their
  version-specific ends are not drawn here.

## Model scope and correction (2026-10-17)

The earlier model drew a canonical TruSeq dual-index library with a dA junction and
TruSeq Read 2 arm, plus an unconditional ligation-index step and a physical discard
of free hash. None is supported by the inspected sources: the hash carries the
Nextera Read 2 handle, ligation applies to sci-RNA-seq3 only, and the main text does
not describe removing free hash. The page now draws the hash oligo and the
hash-derived molecule with the version-specific sci-RNA-seq parts as bracketed role
tokens, one warned boundary for the sci-RNA-seq indexing steps, and unresolved fates
for free hash and transcript-derived molecules. Only the Read 2 site is computed.

## Open questions

- RT primer, sci-RNA-seq3 ligation hairpin, P5/P7 PCR primers, enzymes, cycle
  numbers and sequencing run for each experiment (Science supplement, unread).
- Whether free hash oligo is washed away before pooling (supplement, unread).
- GEO GSE139944 sample protocols (not fetched in this bounded pass).
- Whether the large screen used one or several hashes per treatment population
  (an earlier note said pairs; not found in the inspected main text).
- Meaning of the printed `B` (IDT mixed-base convention presumed, not stated).
