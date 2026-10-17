# ChIL-seq — research notes

Authoritative protocol: Handa *et al.*, "Chromatin integration labeling for mapping
DNA-binding proteins and modifications with low input", *Nature Protocols* 15,
3334–3360 (2020), [10.1038/s41596-020-0375-8](https://doi.org/10.1038/s41596-020-0375-8).
Predecessor method: Harada *et al.*, *Nat Cell Biol* 21, 287–296 (2019),
10.1038/s41556-018-0248-3.

## Sources read

- 🟢 Defining Nature Protocols article, publisher PDF (printed doi:10.1038/s41596-020-0375-8,
  30 pp.; hand-downloaded by the maintainer, MANIFEST row `s41596-020-0375-8.pdf`; identity
  checked against the printed DOI, title and authors 2026-10-17). Read in
  `s41596-020-0375-8.pdf.txt`: Fig. 1 legend (lines 500–555), Fig. 2a legend (705–707),
  library preparation/sequencing overview (943–953), Reagent setup (1192–1310), Procedure
  Steps 7–84 (1500–1990; Box 1 lines 797–844 and Box 2 lines 880–917 only for QC context).
- 🟢 Harada 2019 predecessor main text, publisher PDF (printed doi:10.1038/s41556-018-0248-3;
  hand-downloaded, MANIFEST row `s41556-018-0248-3.pdf`): Methods "ChIL probe" / "ChIL
  protocol" / "ChIL reaction" (`.txt` lines 2505–2630) and sequencing (2673, 3080).
- 🟢 Harada 2019 predecessor Supplementary Information (`related_s41556-018-0248-3_..._MOESM1_ESM.pdf.txt`
  lines 25–56, 91–94, 142–163) and its Supplementary Table 4 (`..._MOESM6_ESM.xlsx.txt`
  lines 1–17), fetched from the Springer supplement URLs 2026-10-17. MOESM2–5 and 7–9
  (reporting summary; tag counts, correlations, P values, predictions, coverage, source
  data) were judged non-chemistry from their legends and not opened.
- 🟢 Author Correction doi:10.1038/s41596-023-00940-6 (Nat Protoc 19, 1288; 2024):
  replaces a mismatched merged image in Fig. 2c only; no chemistry change.
- 🟢 Supplementary Information (MOESM1 PDF): Supplementary Tables 2–4 (`.txt` lines
  193–445; Table 4 checked against the rendered PDF page 6) and captions of
  Supplementary Figs. 1–7 (lines 683–788).
- Not read: Box 3 probe-preparation steps beyond the annealing/conjugation lines
  checked, the Data analysis section, Supplementary Data MOESM2–5, the Protocol Exchange
  step-by-step protocol cited by Harada 2019 (ref. 36), and the SMART-Seq v4 kit manual.

## Procedure (defining article)

- 🟢 Probe binding (Step 7A): primary antibody 2 µg/ml, RT 2 h; three PBS washes;
  secondary-antibody ChIL probe 1 µg/ml in 0.1× Blocking One-P/PBS with 0.5 M NaCl, 4 °C
  overnight (option B: primary-antibody ChIL probe at 2 µg/ml in 1 M NaCl). Step 8: three
  20-min ice-cold PBS washes.
- 🟢 In situ transposition (Steps 11–17): 50 µl Tn5 (90 ng/well) in dialysis buffer, RT
  10 min; then 50 µl of 20 nM annealed MEDS-B, RT 1 h; PBS and dialysis-buffer washes;
  100 µl 1× TAPS-DMF (10 mM TAPS pH 8.5, 5 mM MgCl2, 10% DMF), 37 °C 1 h; 0.2% SDS, RT
  10 min; three PBS rinses. Tn5 is not pre-loaded: it binds the probe's ME first and
  MEDS-B second.
- 🟢 Fill-in (Steps 18–22): 1× T4 DNA ligase buffer with 0.1 mM each dNTP, T4 DNA ligase
  (2 U/µl) and T4 DNA polymerase (0.015 U/µl), RT 30 min; 0.2% SDS 10 min; PBS rinses.
  Harada 2019 describes the same "sealing reaction" (200 U ligase + 1.5 U polymerase in
  ligase buffer, no dNTPs listed; lines 2587–2590); that is a predecessor version, not
  treated as a discrepancy. Gap and overhang lengths at the junction are not printed, so
  the model keeps that step role-level.
- 🟢 In situ transcription (Steps 23–25): T7 RNA polymerase 2 U/µl, 2 mM each NTP, 0.4 U/µl
  RNase inhibitor in 1× T7 buffer, 37 °C overnight.
- 🟢 RNA recovery (Steps 26–41): DNase I in the well, 37 °C 30 min; supernatant + RLT +
  ethanol onto an RNeasy MinElute column; elution in 10 µl water, re-applied once.
- 🟢 First strand (Steps 43–48): Read2+N11 primer replaces the 3′ SMART-Seq CDS Primer IIA;
  72 °C 3 min annealing; SMARTScribe RT with 5× Ultra Low First-Strand Buffer and the
  "SMART-Seq v4 oligonucleotide" (2.4 µM, a template-switching oligo; bases proprietary,
  not printed); 42 °C 90 min, 70 °C 10 min. The overview (line 949) says that above ~10k
  cells the kit components can be used without template switching.
- 🟡 Ad1_noMX anneals within the probe-derived Read 1/ME sequence of the cDNA, so the
  template-switched cDNA end lies outside the PCR product; the TSO is declared as a
  reagent but its product is not drawn.
- 🟢 PCR (Step 49A): SeqAmp, Ad1_noMX and Ad2_index at 0.24 µM; 95 °C 1 min, 12–17 ×
  (98 °C 10 s, 65 °C 30 s, 68 °C 3 min), 72 °C 10 min; 17 cycles for single cells/100/1k,
  14 for 10k, 12 for 100k. Option B (mtChIL-seq) uses Ad1_mt with 50 °C annealing.
- 🟢 Cleanup and size selection: two 1:1 AMPure XP cleanups (first with 10× lysis buffer;
  Steps 50–61); optional MinElute pooling (Steps 63–71, not drawn); 2% E-Gel Size Select II
  recovery of 250–350 bp and MinElute (Steps 72–81); Bioanalyzer and Takara qPCR
  quantification (Steps 82–83).
- 🟢 Sequencing (Step 84): 51 cycles single-end on HiSeq1500/2500 or NovaSeq, 5.26% PhiX;
  dual-index sequencing only for mtChIL-seq. The index-read cycle count for the 8-nt Ad2
  index is not stated, so the page leaves Index 1 unfixed and Read 2 not acquired.
- 🟢 The overview (line 951) says "approximately 250–400 bp of DNA (Supplementary Fig. 5)"
  is analysed, while Step 72 recovers 250–350 bp; the first describes the Bioanalyzer
  profile after selection, so this is not treated as a conflicting instruction.

## Evidence used

- 🟢 Supplementary Table 4 prints the ChIL probe v1 strands (`ChILv1-Fw` 5′-azido,
  `ChILv1-Rv` 5′-phosphate/3′-TMR, 70 nt each), the multitarget v3 probes
  (`ChILv3.1`/`ChILv3.2`, 77 nt; v3.2-Rv carries 3′-Alexa Fluor 488), the MEDS-B
  transposome oligos (`MEDSB-Rv` 5′-phosphate), the `Read2+N11` first-strand primer,
  and the `Ad1_noMX`, `Ad1_mt` and indexed `Ad2_index` PCR primers. The table colours
  the probe into T7 promoter, partial Read 1, ME, v3 index, PCR-primer sequence, extra
  sequence and spacer components.
- 🟢 ChIL probes are made by conjugating azido ChIL-DNA to DBCO-modified IgG
  (Supplementary Fig. 1b); the standard transposition evaluation uses a primary
  antibody and a secondary-antibody ChIL probe (Supplementary Fig. 2, "as described
  in Box 1").
- 🟢 Integration mostly requires Tn5 plus MEDS-B; ChIL-DNA ends alone give weaker
  signal (Supplementary Fig. 2a). Probe-to-genome joining is assayed by genomic PCR with
  `ChILprobeT7-Fw` and `LINE-Rv` (Supplementary Table 2).
- 🟢 T7 RNA polymerase transcribes in fixed cells; only the 92-bp RT-qPCR product was
  obtained, "suggesting that T7 RNA polymerase was not highly processive" there
  (Supplementary Fig. 4d).
- 🟢 For mtChIL-seq, cDNA is made "as per the original ChIL-seq method" and the PCR uses
  `Ad1_mt`, which anneals to the 3′ region of the cDNA, with `Ad2_index`
  (Supplementary Fig. 6).
- 🟡 `Ad1_noMX` ends in the probe's partial Read 1 plus the first ME bases, and
  `Read2+N11` carries the MEDS-B forward strand, so the PCR product is a single-index
  Nextera shell with the probe-derived ME on the Read 1 side and the N11 random-priming
  site next to the Read 2 ME. Read 1, Index 1 and Read 2 sites are computed from that
  library. The N11 is a priming site, not a declared UMI.

- 🟢 Predecessor only (Harada 2019 Supplementary Information; the defining Procedure now
  supplies the drawn values): ChIL probes are IgG conjugated with DBCO-(PEG)5-NHS ester, then
  Cu-free click to the 5′-azide of pre-annealed double-stranded ChIL DNA, with one to
  four DNAs per IgG at a 1:2 IgG:DNA ratio (Supp Fig. 1a,d; MOESM1 lines 27–36);
  default conditions were 0.1 µl Tn5, integration at 37 °C, PEG5 linker and probe
  incubation in 0.5 M NaCl (Supp Table 4 lines 4–17); libraries were ~300 bp
  (~200–400 bp) on a Bioanalyzer (Supp Fig. 2a; lines 49–50). The defining Procedure
  uses the same 90 ng Tn5, 37 °C integration and 0.5 M NaCl probe incubation.

## Model corrections (2026-10-17 bounded pass)

- Removed an unsupported extension step and its invented "ChIL-probe extension
  architecture" oligo input. The interval between integration and in situ T7
  transcription was then an explicit unresolved step (supplement-only pass); the
  defining Procedure has since resolved it as a fill-in (below).
- Removed the stray "ChIL-probe end" placeholder from inside the final insert, and
  added the N11 site.
- Gave the binding, transposome and PCR inputs their Supplementary Table 4 oligos.
- Listed the v3 probes and `Ad1_mt` as oligos for the mtChIL-seq branch, which is not
  drawn.

## Source slips (supplement captions; library unaffected)

- Supplementary Fig. 4e says the Gapdh/Neurod6 promoter primers are in "Supplementary
  Table 4"; they are printed in Supplementary Table 3 (lines 248–262).
- The Supplementary Fig. 2 caption says "Tn5-mediated in situ transcription" for what is
  a genomic-PCR transposition assay.

## Model corrections (deferred follow-up 03, defining article read)

- Replaced the warned "magic happens" boundary between integration and transcription
  with the Procedure's in situ T4 DNA polymerase + T4 DNA ligase fill-in (Steps 18–22).
- Split transposome assembly into Tn5 first and MEDS-B second, with TAPS-DMF/SDS conditions.
- Added DNase I digestion and RNeasy MinElute RNA recovery, the SMART-Seq v4
  template-switching oligonucleotide and SMARTScribe conditions at first strand, SeqAmp
  PCR conditions, two AMPure cleanups and E-Gel 250–350-bp size selection with
  MinElute recovery.
- Declared the acquired reads: Read 1 = 51 cycles single-end; Index 1 acquired for
  demultiplexing with no printed cycle count; Read 2 not acquired.

## Open questions

- Gap/overhang geometry at the Tn5 integration junction before the fill-in (not printed).
- SMART-Seq v4 oligonucleotide sequence (proprietary kit component; kit manual not read).
- Index-read cycle count for the 8-nt Ad2 index.
- The mtChIL-seq final library and how the v3 probe index is read (it is not drawn).

## Protocol boundary

ChIL-seq is kept separate from CUT&Tag because its antibody oligo and T7 transcription
are part of the defining molecular path.
