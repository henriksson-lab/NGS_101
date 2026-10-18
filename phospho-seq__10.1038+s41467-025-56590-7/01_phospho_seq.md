# Phospho-seq

**Defining work:** Blair J. D. et al., "Phospho-seq: integrated, multi-modal profiling of
intracellular protein dynamics in single cells", *Nature Communications* 16, 1346 (2025),
doi:10.1038/s41467-025-56590-7, open article PMC11794950. The directory is named for this
published DOI; the bioRxiv preprint doi:10.1101/2023.03.27.534442 (which named it until
2026-10-17) is superseded and no longer supplies
active claims; earlier source-review passes that read it are history.

## Sources read (bounded pass 2026-10-17)

- PMC11794950 JATS Methods, text lines 213–217 (antibody conjugation and pooling, tag
  oligos), 241–245 (Phospho-seq staining, 10x scATAC v1.1 modifications, tag PCR),
  251 (ASAP-seq comparison), 255 (sequencing), 281–287 (Phospho-seq-multi and its run).
  The hash-oligo poly(A) count was cross-checked in the article PDF text line 1948.
- Cited ASAP-seq article (Mimitou et al. 2021, PMC8763625) JATS line 234 and MOESM1
  lines 253–256 for the BOA/BOB sequences that Phospho-seq names but does not reprint.
- `scrape_primers.py` run on the published JATS text; `--find` for BOA/BOB in the ASAP-seq
  bundle (three exact hits each).

## Molecular boundary

🟢 Antibodies are click-conjugated (mTz-PEG4 antibody, TCO-PEG4 oligo; 2–4 oligos per
antibody), pooled and freed of unconjugated oligo by 40% ammonium sulfate and 50 kDa
MWCO filtration. Fixed, permeabilized cells are blocked for 30 min with 1 µg TotalSeq-A
CD298/B2M hashing antibodies and 100 µg of a 3′-blocked oligo (organoid runs also add
free TotalSeq-A spike-in oligos), pooled, stained 1 h RT with the SSB-preincubated
antibody pool, washed and resuspended in 10x nuclei buffer.

🟢 30,000 cells are tagmented (7 µl ATAC buffer B + 3 µl ATAC enzyme) and the 10x scATAC
v1.1 protocol continues with ASAP-seq changes: 0.5 µl of 1 µM each of BOA and BOB is added
to the barcoding reaction; a 5 min 40 °C step starts the GEM incubation. After silane
elution in 43.5 µl, 40 µl goes to SPRI and 3 µl is set aside. The first-SPRI supernatant
is saved: it holds the hash and antibody tags while ATAC DNA stays on the beads.

🟢 Supernatant + 32 µl SPRI (5 min), 80% EtOH, elute 42 µl EB; combine with the 3 µl
aside and split into TotalSeq-A and TotalSeq-B PCRs: 50 µl 2× KAPA HotStart ReadyMix,
2.5 µl 10 µM P5, 10 µM P7, 22.5 µl H2O, 22.5 µl input; 95 °C 3 min, 14–18 × (95 °C 20 s,
60 °C 30 s, 72 °C 20 s), 72 °C 5 min; 1.6× SPRI, elute 20 µl EB.

## Oligos (published, verbatim)

- 🟢 Hashing tag: `/5AmMC12/CCTTGGCACCCGAGAATTCCA[15 X]B` + 28 A + `*A*A` (30 A in total).
- 🟢 Antibody tag: `/5AmMC12/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[10 N UMI][15 X barcode][9 N UMI]GCTTTAAGGCCGGTCCTAGC*A*A`.
- 🟢 Blocking oligo: 29 N + `/3ddC/` (described as a 30-nt 3′-blocked oligo).
- 🟢 P5 `AATGATACGGCGACCACCGA`; P7 TotalSeq-A `CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`; P7 TotalSeq-B `CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGC` (X = i7).
- 🟢 (ASAP-seq) BOA `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT/3InvdT/`; BOB `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC/3InvdT/`.
- 🟢 Phospho-seq-multi only: HTO additive `CCTTGGCACCCGAGAATT*C*C`, ADT additive `GTGACTGGAGTTCAGACGTGTGC*T*C`, SI-PCR P5 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGC*T*C`.

Model corrections in this pass: the hashing-tag poly(A) was drawn as 30 A + `*A*A`
(published: 28 A + `*A*A`) and the blocking oligo as 30 N + `/3ddC/` (published: 29 N +
`/3ddC/`). These were repository transcription errors, not publication discrepancies.

## Sequencing

🟢 Pooled 25% ADT, 10% HTO, 65% ATAC; NextSeq 550, 75-cycle High Output: 34 bp Read 1,
8 bp i7, 16 bp i5, 34 bp Read 2. ATAC by CellRanger-ATAC v2.0.0; tags by alevin.

🟢 Phospho-seq-multi (10x Multiome; not drawn): only BOB is spiked in (TSA tags are captured
by direct poly(A) capture); additives in pre-amplification; ATAC/ADT libraries 50 bp Read 1,
8 bp i7, 86 bp Read 2, 16 bp i5 (8 dark cycles); cDNA 40/8/44.

## Remaining questions

- 🟡 The tag-library cell-barcode end and the bridge-capture junctions depend on the 10x
  scATAC v1.1 gel-bead oligo, which is not printed; the role order drawn is inferred.
- The 10x scATAC v1.1 user guide (steps 3.1o, 3.2d and the ATAC library) and the author
  step-by-step PDFs (phospho-seq.com protocol_2 and protocol_3; manifest rows record curl
  error 6) were not read.
- The CITE-seq hyper-conjugation PDF (conjugation chemistry only) was not read.
- Model correction: earlier versions of this page and note cited "Chung et al."; both the
  preprint and the published article list John D. Blair as first author (Satija lab).

## Local validation migration (2026-10-18)

Re-read published PMC11794950 Methods lines 211–245 (printed tags, blocking oligo
and separate TotalSeq-A/B PCRs) and the cited ASAP-seq PMC8763625 Results lines
122, 132, Methods line 234 and Extended Data Fig. 1b legend line 332. Ran the bundle
oligo scraper and verified the printed nucleotide/modification passages.

🟢 The two tag endpoints now retain different identifier provenance. The printed
TotalSeq-B tag has a 15-base antibody barcode between a 10-base and a 9-base UMI
component, represented as one multipart molecular identifier. The TotalSeq-A
hashing tag has a 15-base hash barcode and no intrinsic UMI; its captured library
acquires the 10-base BOA-derived unique bridging identifier (UBI). This inherited
BOA mechanism is explicitly scoped to the bridge method Phospho-seq cites, rather
than supplying a generic antibody UMI to every endpoint. The blocking oligo's
terminal `/3ddC/` is included as its C nucleotide with a blocked end annotation.

🟡 Tag/bead junction and kit ATAC geometry remain unassembled; each PCR limitation
uses its actual incoming branch template. Source coverage is an explicit local
inventory of the inspected Methods' operations and named kit delegations, never
created by iterating graph IDs. The vendor guide itself remains uninspected.
Candidate non-library tagmentation and empty/ambiguous GEM outcomes have unresolved
physical fates. This migration is not a completed chemistry/source-absence audit.
