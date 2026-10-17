# TAIL-seq — research notes

Defining source: Chang et al., 2014, *Molecular Cell* 53:1044–1052,
[doi:10.1016/j.molcel.2014.02.007](https://doi.org/10.1016/j.molcel.2014.02.007).

## Sources read

- 🟢 Published article PDF (cell.com, 9 pages; printed title, authors and DOI checked),
  hand-downloaded by the maintainer 2026-10-17: Results "TAIL-seq: The Challenges and
  Solutions" (features 1–8), Figure 1 legend and Experimental Procedures.
- 🟢 Supplemental Information mmc1 (Elsevier CDN; printed title checked): Extended
  Experimental Procedures "Construction of TAIL-seq library", "Sequence processing" and
  "Oligonucleotide sequences for TAIL-seq" (text lines 1105–1118, 1243–1258, 1349–1535).
- 🟢 Author-deposited GEO records [GSM1242324](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1242324)
  and [GSM1242325](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1242325)
  (GSE51299). Their protocol and processing text repeats the mmc1 wording.
- Not read: mmc2/mmc3 (data tables); Illumina HiSeq 2500 reagent documents.

## Oligonucleotides (mmc1, verbatim)

- 🟢 5′ adapter (RNA): `5′-GUUCAGAGUUCUACAGUCCGACGAUC-3′`.
- 🟢 3′ adapter: `5′-/5App/CTGACNNNNNNNNNNNNNNNTGGAATTCTCGGGTGCCAAGGC/iBiodT//iBiodT//3ddC/-3′`
  (legend: /5App/ 5′ adenylation, internal biotinylated dT, 3′ dideoxy-C).
- 🟢 RT primer: `5′-GCCTTGGCACCCGAGAATTCCA-3′`.
- 🟢 PCR primer (forward): `5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA-3′`.
- 🟢 PCR primers (reverse) 1–9: `5′-CAAGCAGAAGACGGCATACGAGAT[index]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′`
  with indexes CGTGAT, ACATCG, GCCTAA, TGGTCA, CACTGT, ATTGGC, GATCTG, TCAAGT, TACAAG.
  Primer 1's last three bases follow the page-14 number in the extracted text. The graph
  therefore uses primer 2 as the representative, so that the source guard can match it.
- 🟢 Spike-ins 0–118 (synthetic poly(A) controls) carry a fixed `CTGACGAGCTACTGTTGG…`
  3′ flank instead of the degenerate window; they are not drawn.

## How they interlock (computed in `tools/tail_seq.py`)

- 🟡 The RT primer is the exact reverse complement of the 3′-adapter arm
  `{{= tail_seq.RTP_SITE }}`. The forward primer ends on the first 21 nt of the 5′-adapter
  DNA copy. The reverse primer's 3′ 21 nt equal the RT primer without its 5′ G. At that
  position the reverse primer has T, so the amplified library carries A opposite where the
  adapter had C.
- 🟡 Final library (top strand): P5 · 5′-adapter copy · insert · poly(A) · terminal
  additions · `CTGAC` · N15 · `{{= tail_seq.RPI_ARM }}` · index (rc) · P7′.
- 🟡 No sequencing primer is printed. The page derives Read 1 from the 5′-adapter sequence
  (read 1 must start in the insert) and Read 2 from the RT-primer sequence. That placement
  makes read-2 cycles 1–15 the degenerate bases, as mmc1 states, and puts `GTCAG` at
  cycles 16–20 (mmc1 searches for it starting between cycles 14 and 16).

## Step by step

🟢 Total RNA (TRIzol, RNeasy MinElute; >200 nt fraction) → Ribo-Zero depletion (Epicentre
MRZH11124) → 3′-adapter ligation → partial RNase T1 digestion at low concentration (Ambion
AM2283; cleaves after G, preserving poly(A)) → streptavidin pull-down (Invitrogen 11206D) →
PNK phosphorylation (Takara 2021B) → gel purification at 500–1000 nt → 5′-adapter ligation →
RT (Invitrogen 18080085) → Phusion PCR (Thermo F-530L) → PCR products purified again.
Ligases, reaction conditions and PCR cycle numbers are not stated.

🟡 RNase T1 leaves 5′-OH on downstream fragments, which is why the PNK step precedes 5′
ligation. This is a mechanistic reading; the source states only the order.

## Read layout

🟢 HiSeq 2500, 51 × 251 paired-end, with 20% PhiX and 1% spike-in mixture (mmc1;
Experimental Procedures). Read 1 (51 nt) identifies the transcript. Read 2 begins with the
15 degenerate bases (deduplication with read-1 cycles 21–35), then `GTCAG`. The article's
"231 nt from the 3′ end" equals 251 cycles minus those 20 adapter bases (🟡 arithmetic).
Tail length is called from relative T fluorescence with a Gaussian-mixture HMM, because
base calls are unreliable through long T stretches.

## Open questions

- Ligases and conditions for both adapter ligations; RT and PCR conditions and cycles.
- Whether the 6-nt index was read, and which reverse primer each sample used.
- The sequencing primers actually loaded (derived here, not printed).

## Review trail

- 2026-10-16: first passes used a third-party article mirror (now HTTP 403).
- Round 2026-10-17-breadth-first-20: no original retrievable; no chemistry change.
- Round 2026-10-17-deferred-followup-02: GEO GSM1242324/5 used for reaction order,
  reagents and the read-2 adapter layout; inferred TruSeq amplicon shell.
- Round deferred-followup-03 (2026-10-17): article PDF and mmc1 read. Replaced the
  inferred TruSeq shell, inferred adapters and primers with the mmc1 oligos. Both
  ligations now have known single-strand end geometry (5′-adenylated, ddC-blocked
  biotin DNA 3′ adapter; RNA 5′ adapter). Reads 1/2 come from adapter-derived sites;
  Index is no longer listed as acquired.


## Critical re-review, 2026-10-17

Restored nucleotide-bearing /iBiodT//iBiodT//3ddC/ as TTC in the 3′ adapter state, including its blocked end. Variable native poly(A)/terminal additions no longer assert fixed lengths. RNA intermediates are single-stranded scenes. RT has a distinct source-limited cDNA state; PCR unavailable geometry refers to that actual input. The modified adapter suffix lies beyond the RT binding arm and is not appended to the final amplicon. Fixed the shared scraper terminal-modification parser and added a regression.

Inspected scope: Published mmc1 .txt oligonucleotides 1349–1365, modification legend 1534–1535, sequence processing 1243–1258; full-bundle primer scraper run. This bounded correction does not complete the linked-source or two-reader absence audit.
