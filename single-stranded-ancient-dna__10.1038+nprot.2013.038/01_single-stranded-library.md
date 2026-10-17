# Single-stranded library preparation for ancient or damaged DNA

## 1. What it is

This protocol recovers individual highly degraded DNA strands instead of requiring an intact duplex with two usable ends.

## 2. Sources and evidence

- 🟢 Gansauge MT, Meyer M. *Single-stranded DNA library preparation for the sequencing of ancient or damaged DNA.* Nature Protocols 8:737–748 (2013). [doi:10.1038/nprot.2013.038](https://doi.org/10.1038/nprot.2013.038).
- 🟢 A phosphorylated, biotinylated adapter is ligated to the 3′ ends of heat-denatured DNA, captured on streptavidin beads, copied, and joined to a second adapter by blunt ligation.
- 🟢 Main article read (round deferred-followup-03): maintainer hand-downloaded PDF
  `nprot.2013.038.pdf` (13 pages; printed title, authors, DOI and "corrected after print
  15 May 2013" checked; ends with the Corrigendum page). Sections read in its `.txt`:
  Introduction/overview lines 85–105 and 270–290, Applications 450–463, Reagent setup and
  Table 1 541–748 (Table 1 also inspected on the rendered p. 741), Procedure Steps 1–34
  lines 749–1110, Anticipated results 1222–1232, Corrigendum 1352–1366. Earlier passes
  (below) used only the public article page, supplements and the 2017 oligo table.
- 🟢 Fig. 1 (publisher image `41596_2013_Article_BFnprot2013038_Fig1_HTML.jpg`, inspected;
  legend "A, deoxyadenine; Btn, biotin; Pho, phosphate; U, deoxyuracil") maps the procedure:
  Steps 1–4 leave strands without the drawn 5′ phosphate and broken at the deoxyuracil, then
  separated; Steps 5–7 join the 5′-Pho, spacer-and-biotin adapter to strand 3′ ends; Steps
  8–15 capture on streptavidin beads and extend a primer with a 5′ tail, the copy carrying a
  3′ A; after Steps 16–19 that end is blunt; Steps 20–24 ligate a duplex adapter whose short
  strand has a blocked 3′ end; Step 25 releases the molecules from the beads; Steps 26–29
  are qPCR and Steps 30–32 indexing PCR (drawn with a 7-N index on the P7-side primer only).
  The enzymes for Steps 1–4 and 16–19 are not named in the figure.
- 🟢 Supplementary Fig. 1 (MOESM473; rendered page inspected) draws the double-indexed
  single-stranded library as P5 · 7-nt index · `ACACTCTTTCCCTACACGACGCTCTTCC` · insert ·
  `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` · 7-nt index · P7, with arrows labelled "CL72 seq",
  "Read 2 seq", "Read 3 seq" and "Read 4 seq". The caption says the P5 adaptor has a 5-bp
  deletion relative to a regular TruSeq library.
- 🟢 Table 1 prints CL72 "Sequencing primer" `{{= ss_adna.CL72 }}`, and Step 34 says
  "replace the sequencing primer of the first read by the custom primer CL72" (10 µl of
  100 µM in 1,990 µl hybridization buffer). This replaces the earlier Fig. 1-based
  inference with the same bases.
- 🟡 The index-2 primer `{{= ss_adna.READ4_INDEX2 }}` remains inferred from the
  Supplementary Fig. 1 "Read 4 seq" arrow: Table 1 and Step 34 print no index-2 oligo and
  defer to Illumina's multiplex sequencing instructions. (An earlier note named it
  "Gesaffelstein"; that is not supported by any inspected original.)
- 🟢 Table 1 (corrected) prints, matching the 2017 table: CL9
  `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`; CL53 `CGACGCTCTTC-ddC`; CL73
  `[Phosphate]GGAAGAGCGTCGTGTAGGGAAAGAG*T*G*T*A` (`*` phosphothioate); CL78
  `[Phosphate]AGATCGGAAG[C3Spacer]10[TEG-biotin]`; IS4 "P5 primer without index"
  `{{= ss_adna.IS4 }}`; positive-control template CL104 and qPCR oligos CL105–CL108,
  IS7, IS8 (quantification only, not drawn). Indexing primers: "See Supplementary Table 2".
- 🟢 Gansauge et al. 2017, *Nucleic Acids Research* 45:e79, Supplementary Table 2
  (PMC5449542 `gkx033_supp.pdf`) lists, each marked "already published in Gansauge &
  Meyer 2013": CL9 `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` (extension primer), CL53
  `CGACGCTCTTC-ddC` and CL73 `Pho-GGAAGAGCGTCGTGTAGGGAAAGAG*T*G*T*A` (double-stranded
  ligation adapter; `*` phosphorothioate) and CL78 `Pho-AGATCGGAAG[C3Spacer]10-TEG-biotin`
  (single-stranded ligation adapter). These are not interchangeable with that paper's new
  T4-ligase splinter oligos.
- 🟢 Published correction (publisher article page, Change history, 15 May 2013): "the
  sequence of oligonucleotide CL78 reported in Table 1 as “(AGATCGGAAG[C3Spacer]10[TEG-biotin]
  …)” is incorrect. The correct sequence is “[Phosphate]AGATCGGAAG[C3Spacer]10[TEG-biotin]
  …”." The 2017 table agrees; the model uses the phosphorylated form and the page shows
  both values as a source discrepancy.
- 🟢 The 2013 Supplementary Table 2 (MOESM478) prints P5 indexing primers
  `AATGATACGGCGACCACCGAGATCTACAC[7 nt]ACACTCTTTCCCTACACGACGCTCTT` and P7 indexing primers
  `CAAGCAGAAGACGGCATACGAGAT[7 nt]GTGACTGGAGTTCAGACGTGT`, differing only at the index
  window; there is no single universal index-primer pair.

## 3. Interlocks (🟡 computed)

- CL9 ends in the reverse complement of the CL78 bases, so it primes copying from the
  captured adapter; CL78 plus the rest of the CL9 complement forms the P7-side adaptor
  `AGATCGGAAG…` of Fig. 1 (checked in `tools/ss_adna.py`).
- CL53 (12 nt including the 3′ ddC) is complementary to the first 12 nt of CL73, so the
  duplex adapter has a blunt end formed by the 5′-phosphorylated CL73 and the
  dideoxy-blocked CL53. Only CL73 can be joined to the copied strand's 3′ end
  (stated in the article Introduction); its reverse complement is the P5-side adaptor `TACACTCTTTCC…CTTCC`.

## 4. Procedure (🟢 Steps 1–34 of the corrected article)

- Steps 1–2: endonuclease VIII (0.12 U/µl) with optional Afu UDG (0.02 U/µl) in
  CircLigase buffer II with MnCl2, 37 °C 1 h (UDG optional; endonuclease VIII cleaves
  abasic sites). Steps 3–4: FastAP 1 U, 37 °C 10 min, then 95 °C 2 min and ice-water.
- Steps 5–7: 20% PEG-4000, 0.125 µM CL78, CircLigase II 5 U/µl, 60 °C 1 h; stop solution.
- Steps 8–11: 95 °C 1 min, bind to MyOne C1 streptavidin beads 20 min at RT; wash buffers A
  and B.
- Steps 12–15: CL9 2 µM, 65 °C 2 min then ice; Bst 2.0 polymerase (24 U) added at 15 °C,
  ramp 1 °C/min to 37 °C, 5 min at 37 °C; wash A, stringency wash at 45 °C, wash B. The
  Introduction says this wash removes excess primer, including primer bound to unligated
  adapter, which would otherwise cause adapter dimers.
- Steps 16–19: T4 DNA polymerase (5 U) in Buffer Tango, 25 °C 15 min, EDTA; washes.
  The Introduction says Bst 2.0 lacks 3′–5′ exonuclease activity and leaves 3′ overhangs.
- Steps 20–23: CL53/CL73 adapter (hybridized 95 °C 10 s, 0.1 °C/s to 14 °C) at 2.5 µM,
  5% PEG-4000, T4 DNA ligase 10 U, 1 h RT; washes. "Only one adapter strand is ligatable,
  whereas the other is blocked by a 3′-terminal dideoxy modification."
- Steps 24–25: 25 µl TET, 95 °C 1 min, supernatant contains the library molecules. The
  Introduction notes that a fraction of biotinylated molecules, including the original
  template strands, is also released but lacks priming sites and is not amplified; the
  graph records that as a computational/amplification fate, not a physical separation.
- Steps 30–32: AccuPrime Pfx, 400 nM P7 indexing primer and P5 indexing primer or IS4,
  95 °C 2 min, then cycles of 95 °C 15 s / 60 °C 30 s / 68 °C 1 min; MinElute or AMPure XP
  purification, 20 µl TE. Step 33 Bioanalyzer; Step 34 sequencing with CL72.

## 5. Molecular path

Every denatured template strand is independently captured through its 3′ end. Copying converts it to a duplex only after capture, avoiding the need for complementary ancient strands to survive together.

The dual-adapter molecule before indexing PCR is not the final library: it carries no
P5, P7 or index windows, and the CL53 strand is not ligated.

## 6. Open questions

- Exact index-2 sequencing primer oligo (inferred from Supplementary Fig. 1; the protocol
  defers to Illumina instructions) and the acquired read/index lengths (not specified;
  the worked example is a paired-end MiSeq run with overlap merging).
- The qPCR quantification branch (Steps 26–29, IS7/IS8 and CL107/CL108 assays) and the
  CL104 positive control are not drawn.
- PCR cycle number is sample-specific (chosen from the Step 29 qPCR plots); Supplementary
  Fig. 4 gives 11 and 16 cycles only for its worked example.


## Critical re-review, 2026-10-17

Released original template strands travel with the supernatant into PCR. The source-supported absence of amplification sites is distinguished from physical discard; later recovery remains unresolved.

Inspected scope: Existing source-qualified Introduction and Procedure Steps 24–25 declarations; no new reaction-source assertion. This bounded correction does not complete the linked-source or two-reader absence audit.
