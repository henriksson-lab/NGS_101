# s3-WGS — single-cell whole-genome sequencing by symmetrical-strand combinatorial indexing

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available. Relationships marked 🟡 *(computed)* were worked out with
> `lib/` while writing this note; they are not yet asserted in a self-test, because this
> protocol has no `tools/` module yet.

**s3-WGS** — Mulqueen RM, Pokholok D, O'Connell BL, Thornton CA, Zhang F, O'Roak BJ, Link J,
Yardımcı GG, Sears RC, Steemers FJ, Adey AC. "High-content single-cell combinatorial
indexing." *Nature Biotechnology* 39, 1574–1580 (2021).
doi:[10.1038/s41587-021-00962-z](https://doi.org/10.1038/s41587-021-00962-z)
(PMID 34226710, PMC8678206). Data: GEO GSE174226.

The same paper introduces three assays on one chemistry: **s3-ATAC** (catalogue slug
`s3-atac__10.1038+s41587-021-00962-z`), **s3-WGS** (this note) and **s3-GCC** (WGS plus
chromatin conformation; not a separate catalogue entry). Stepwise protocols are on
protocols.io: s3-WGS dx.doi.org/10.17504/protocols.io.beb3jaqn, s3-GCC
dx.doi.org/10.17504/protocols.io.beb4jaqw (not fetched, see below). s3-WGS is the "s3"
upgrade of the authors' earlier **sci-DNA-seq** (Vitak et al. 2017, *Nat Methods*, the
paper's ref. 10), whose xSDS nucleosome-depletion step it keeps.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/s3-wgs__10.1038+s41587-021-00962-z/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41587-021-00962-z_PMC8678206.html.txt` | author manuscript (PMC) incl. Online Methods | design rationale, s3-ATAC and s3-WGS methods, sequencing |
| `s41587-021-00962-z_41587_2021_962_MOESM6_ESM.xlsx` (+ `.txt`, made by hand with `tools/doctext.py`) | Supplementary Tables 1–8 | **every oligo** (Tables 2–5) and the stepwise molecular layout (Table 1) |
| `s41587-021-00962-z_41587_2021_962_MOESM1_ESM.pdf.txt` | Supplementary Information (figures) | only Supp. Fig. legends; no chemistry |
| `upstream_s3-WGS.html.txt` | scg_lib_structs s3-WGS page | secondary account; it is two sentences and points to the s3-ATAC page |
| `../s3-atac__10.1038+s41587-021-00962-z/upstream_s3-ATAC.html.txt` | scg_lib_structs s3-ATAC page (sibling source dir) | the upstream drawing the s3-WGS page defers to; checked below |
| `MOESM7–13` | source data for figures | not chemistry |

Could not be used:

| What | Why |
|---|---|
| protocols.io s3-WGS (dx.doi.org/10.17504/protocols.io.beb3jaqn) | not fetched by `get_sources.py`; the stepwise bench protocol — volumes per well, buffer recipes beyond the methods 🔴 |
| `MOESM2_ESM.pdf` | the reporting summary; its text twin is empty (print-to-PDF image) |
| `supp_NIHMS1707134-supplement-*` | each is an HTML error page (~21 kB), not the file |

---

## 1. What it is

A **two-level combinatorial indexing** single-cell whole-genome method (sci-), plate based,
no beads, no UMI. 🟢 (Online Methods)

- **Level 1**: nuclei are tagmented in 96 wells, each with a Tn5 loaded with a **different
  8-nt barcode** (`Tn5_96plx_idx`).
- **Level 2**: nuclei are pooled, flow-sorted in small numbers into a PCR plate, and each
  PCR well gets a unique **i7 × i5** primer pair.

Cell ID = Tn5 barcode + i7 + i5 🟢 (upstream s3-ATAC; the paper's "cellID", demultiplexed
by scitools).

What makes it **s3** ("symmetrical strand sci") 🟢 (main text, Fig. 1b; Supp. Table 1):
standard Tn5 libraries need *different* adaptors at the two ends, so a fragment tagged
with the same adaptor at both ends (about half, in the classic A/B Tn5 picture — see
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)) is lost to suppression PCR.
s3 instead loads **one** adaptor species, so every fragment has the same end at both
sides, then **switches** one end to a second adaptor after tagmentation:

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Single-adaptor Tn5** carrying a **deoxyuridine immediately 5' of the ME** | — (here the U is a stop for the gap-fill polymerase) |
| 2 | Gap fill with a **uracil-intolerant** polymerase (Nextera NPM) — it stops at the U, so the barcode and TruSeq handle are **not** copied onto the other strand | — |
| 3 | **Adapter switching** with a 3'-inverted-dT-blocked **LNA** ME′ oligo carrying an A14 (Nextera s5) 5' overhang, run for **10 thermocycles** | — |
| 4 | PCR with a **uracil-tolerant** polymerase (NEBNext Q5U) reads through the U | — |
| 5 | Read layout uses **stock sequencing primers** (Nextera read 1 / index 2, TruSeq read 2 / index 1) | — |

**s3-WGS vs s3-ATAC** (the only differences the paper gives) 🟢: nuclei are **lightly
formaldehyde-fixed** and **SDS nucleosome-depleted** before tagmentation, so Tn5 reaches
the whole genome instead of open chromatin; nuclei are loaded more dilute (500 vs 1,400
nuclei/µL); and the indexing PCR needs **fewer cycles (13–15 vs 16–18)**, which the
authors attribute to more capture events per nucleus. Everything from tagmentation to
sequencing is "as described for s3-ATAC". The upstream s3-WGS page says the same in one
sentence (crosslinked, nucleosome-depleted nuclei; library identical to s3-ATAC) — agrees.

## 2. Oligos

🟢 Verbatim from Supplementary Tables 2–5 (`MOESM6_ESM.xlsx`, sheets SuppTable2–5), as the
paper writes them. IDT notation: `/ideoxyU/` = internal deoxyuridine, `/5Phos/` = 5'
phosphate, `+N` = LNA base, `/3InvdT/` = 3' inverted dT. The tables split each oligo into
named columns; the name of each column is given as a comment.

Transposome (Supp. Table 2; 96 oligos `SBS12_18_UME_sci_1..96`, wells A1–H12, one barcode
each; first and last shown):

```
SBS12_18_UME_sci_1    CGTGTGCTCTTCCGATCT GAACCGCG /ideoxyU/AGATGTGTATAAGAGACAG
SBS12_18_UME_sci_96   CGTGTGCTCTTCCGATCT ATTGTGAA /ideoxyU/AGATGTGTATAAGAGACAG
                      [Truseq_R2_SBS12_partial] [Tn5_96plx_idx] [U-ME]

Nextera_Mosaic_End_REVCOMP   /5Phos/CTGTCTCTTATACACATCT
```

Adapter-switching oligo (Supp. Table 3; called "A14-LNA-ME" in the methods):

```
A14_ME LNA   TCGTCGGCAGCGTC AGATGTGTA+TA+AG+AG+AC+AG/3InvdT/
             [Nextera_R1_A14] [U-ME]
```

Indexing PCR, i7 side (Supp. Table 4; 32 primers: `PCR_i7_P7.S701..S708`,
`PCR_i7_T267..T278_IPE2F`, `PCR_i7_T179..T190_IPE2F`; first shown):

```
PCR_i7_P7.S701   CAAGCAGAAGACGGCATACGAGAT TCGCCTTA GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
                 [i7_Flowcell_Primer] [i7_pcr_idx] [Truseq_R2_SBS12_full]
```

Indexing PCR, i5 side (Supp. Table 5, titled "i7 Molecule side" in the source — a copy-paste
slip, the columns are i5; 64 primers `PCR_A..H_i5_A..H`; first shown):

```
PCR_A_i5_A   AATGATACGGCGACCACCGAGATCTACAC CCTTAAGA TCGTCGGCAGCGTC
             [i5_Flowcell_Primer] [i5_pcr_idx] [i5_Nextera_A14_partial]
```

### How they interlock — 🟡 (computed with `lib/`)

- `Truseq_R2_SBS12_partial` (18 nt) is the **3' 18 nt of `illumina.TRUSEQ_READ2`**
  (34 nt); the missing 5' part is `GTGACTGGAGTTCAGA`. The i7 PCR primer's
  `Truseq_R2_SBS12_full` is exactly `illumina.TRUSEQ_READ2`, so the PCR primer anneals
  over the 18-nt partial and **adds** the remaining 16 nt.
- `U-ME` minus the U is exactly **`nextera.ME`** (19 nt); `Nextera_Mosaic_End_REVCOMP` is
  exactly **`nextera.ME_RC`** — the usual 5'-phosphorylated bottom strand of a Tn5
  adaptor.
- Transposon oligo length: 18 + 8 + 1 (dU) + 19 = **46 nt**.
- `A14_ME LNA`, read as bases, is exactly **`nextera.ADAPTOR_S5`** = `nextera.READ1_PRIMER`
  (`S5` + `ME`, 33 nt). "A14" is the Nextera s5 / read-1 entry point. Its 3' half (the ME)
  is identical in sequence to the Tn5-loaded ME, so it pairs with the **copied ME′** at
  the 3' end of each gap-filled strand; the 3′ inverted dT keeps it from being extended
  itself, and the 6 LNA bases (all in the ME half) raise its Tm so it outcompetes the
  intramolecular ME/ME′ snap-back. 🟡 (the LNA rationale is the paper's 🟢; that the
  ME′ it pairs with is the gap-fill copy follows from Supp. Table 1)
- i7 primer = **`illumina.P7`** + 8-nt index + **`illumina.TRUSEQ_READ2`**, 66 nt.
  Index orientation: `S701` carries `TCGCCTTA`, whose reverse complement `TAAGGCGA` is the
  i7 base sequence the sequencer reports for Illumina index N701/S701 — i.e. the primer
  holds the **reverse complement** of the named index, as in Nextera kits. Checked for
  S701 only; the `T…_IPE2F` names do not say what is read. 🟡
- i5 primer = **`illumina.P5`** + 8-nt index + **`nextera.S5`**, 51 nt; identical to
  `nextera.n5xx_primer(index)` — a standard Nextera N5xx-type primer. 🟡
- Combinatorial space as listed: 96 Tn5 barcodes × 32 i7 × 64 i5. 🟡 (count of the
  tables; how many were used per experiment is in the paper's plate layouts, Fig. 2a/4b)

## 3. Step by step

Conditions 🟢 from the Online Methods ("s3-ATAC library generation" for the shared steps,
"s3-WGS library generation" for the nuclei). Molecular states 🟢 from Supp. Table 1
unless marked.

0. **Transposomes**: 96 uniquely indexed complexes assembled "using previously-described
   methods" (ref. 11 — the loading protocol itself is not repeated 🔴); diluted to 2.5 µM in
   50 % glycerol, 100 mM NaCl, 50 mM Tris pH 7.5, 0.1 mM EDTA, 1 mM DTT; −20 °C.
1. **Nuclei** (cell lines: PDAC patient-derived lines, GM12878): wash 2× in ice-cold PBS
   (adherent lines trypsinised with TrypLE, 15 min 37 °C), pellet, wash twice in NIB-HEPES
   (10 mM HEPES-KOH pH 7.2, 10 mM NaCl, 3 mM MgCl₂, 0.1 % IGEPAL CA-630, 0.1 % Tween,
   protease inhibitors), count; aliquots of 10⁶ nuclei.
2. **Light fixation**: 10⁶ nuclei in 5 mL NIB-HEPES + 246 µL 16 % formaldehyde
   (**0.75 % f.c.**), 10 min, 50 rpm orbital. Pellet; **quench** by resuspending in 1 mL
   NIB-Tris (same recipe with 10 mM Tris-HCl pH 7.4 instead of HEPES — Tris is the
   quencher).
3. **Nucleosome depletion (xSDS)**: wash in 1× NEBuffer 2.1, resuspend in 760 µL, add
   40 µL 1 % SDS (≈ 0.05 % f.c. 🟡 computed), 37 °C 20 min on a ThermoMixer (speed given
   as "300 rcf" in the methods, presumably rpm 🟡). Pellet, 50 µL
   NIB-Tris, count, dilute to **500 nuclei/µL**.
4. **Indexed tagmentation**: 420 µL nuclei + 540 µL 2× TD buffer (Nextera XT); 8 µL per
   well (≈ 1,750 nuclei/well for s3-WGS 🟡 computed; the "~5,000" in the methods is for
   s3-ATAC at 1,400 nuclei/µL) + 1 µL 2.5 µM indexed transposome; **55 °C 10 min**, then
   ice. Every fragment now carries the **same** adaptor at both ends:
   `SBS12-partial · Tn5 idx · U · ME` joined to the 5' end of each strand, with a gap
   (9 nt for Tn5 🟡, general Tn5 knowledge; Supp. Table 1 just draws `<GAP>`) opposite the
   non-transferred ME′.
5. **Pool, DAPI, sort**: pool wells, add 2 µL 5 mg/mL DAPI, sort single nuclei (or a set
   number per well) by DAPI on a Sony SH800 into a 96-well plate with 9 µL 1× TD per well;
   spin 500 rcf 5 min.
6. **Denature** Tn5 and nucleosomes: 1 µL 0.1 % SDS per well (~0.01 % f.c.).
7. **Gap fill, uracil-intolerant**: 4 µL NPM (Nextera XT), **72 °C 10 min**. The
   polymerase fills the 9-nt gap and copies the ME, but **stops at the dU**, so the
   3' end of each strand gains **ME′ only** — not the barcode or the SBS12 handle. 🟢
   (Supp. Table 1: "Use Uracil to block any pcr amplification, just want to close the gap")
8. **Adapter switching**: 1.5 µL 1 µM A14-LNA-ME; 98 °C 30 s, then **10 cycles** of
   98 °C 10 s / 59 °C 20 s / 72 °C 10 s (still NPM). The LNA oligo anneals with its ME to
   the new 3' ME′ and templates extension, adding **s5′ (A14′)** to that 3' end. Its own
   3' end is blocked, so it is not extended. Result per strand:
   `5'- SBS12-partial · Tn5 idx · U · ME · gDNA · ME′ · A14′ -3'` 🟢 (Supp. Table 1,
   "Pre-PCR Molecule"). Supp. Table 1 draws only the top strand; that the bottom strand is converted the same
   way follows from the symmetric adaptor 🟡. The main text states that all fragments end up
   with both adaptors and that the repeated thermocycling maximises reverse-adaptor
   incorporation 🟢.
9. **Quench SDS**: 1 % Triton X-100 (volume not stated in the methods 🔴). Plates can be
   stored at −20 °C for weeks.
10. **Indexing PCR, uracil-tolerant**: 16.5 µL sample + 2.5 µL 10 µM i7 primer + 2.5 µL
    10 µM i5 primer + 3 µL water + 25 µL NEBNext Q5U 2× + 0.5 µL 100× SYBR Green I.
    98 °C 30 s; **13–15 cycles for s3-WGS** (16–18 for s3-ATAC) of 98 °C 10 s, 55 °C 20 s,
    72 °C 30 s, read, 72 °C 10 s; stopped by qPCR inflection. The i5 primer (s5 tail)
    anneals to A14′ first and Q5U **reads through the dU**; the i7 primer (full TruSeq
    read 2 tail) anneals to the copy of SBS12-partial. 🟢 (Supp. Table 1)
11. **Clean-up**: pool 25 µL/well, QIAquick column, elute 50 µL; 1× SPRI (Mag-Bind
    TotalPure), elute 31 µL; Qubit, TapeStation D5000; dilute to 1 nM using the
    100–1000 bp range.

## 4. Final library — 🟡 (assembled with `lib/` from the oligos above; agrees with Supp. Table 1's last row and with the upstream s3-ATAC drawing)

Top strand, 5'→3' (P7 side first, as Supp. Table 1 writes it):

```
5'- P7 · i7′(8) · TruSeq Read 2 (34) · Tn5 idx (8) · T · ME (19) · <gDNA> · ME′ (19) · s5′/A14′ (14) · i5′(8) · P5′ -3'
```

Same molecule from the P5 end (the orientation upstream draws):

```
5'- P5 · i5 (8) · s5/A14 (14) · ME (19) · <gDNA> · ME′ (19) · A · Tn5 idx′ (8) · TruSeq Read 2′ (34) · i7 (8) · P7′ -3'
```

- 164 nt of adaptor in total around the insert (24 + 8 + 34 + 8 + 1 + 19 + 19 + 14 + 8 + 29).
- The position of the dU becomes an ordinary **T** (A on the other strand) after Q5U
  copying — Supp. Table 1 writes it `/T/` and `/A/`; upstream writes the single base `A`
  between ME′ and the Tn5 barcode. Agree. 🟢 / 🟡
- Unlike a standard Nextera library, the two ends are **asymmetric in kind**: the P5 end
  is Nextera (s5 + ME), the P7 end is TruSeq read 2 + Tn5 barcode + ME.

## 5. Sequencing

🟢 NextSeq 500 (High/Mid 150-cycle kits) and, for deeper runs, NovaSeq S2; paired end,
**read 1 = 85, read 2 = 85, index reads 10 + 10** (stated for s3-ATAC; s3-WGS is
"sequenced as described previously"). The paper states that stock primers suffice:
**Nextera read 1 and index 2, TruSeq read 2 and index 1** 🟢.

Read content 🟡 (computed by walking the primers along the construct above):

| Read | Primer | First bases read |
|---|---|---|
| Read 1 | `nextera.READ1_PRIMER` (s5 + ME) | **genomic DNA** directly |
| Index 1 (i7) | `illumina.INDEX1_PRIMER` | i7 as named (`S701` → `TAAGGCGA`), then 2 extra cycles |
| Index 2 (i5) | `nextera.INDEX2_PRIMER` (ME′ + s5′) | `CCTTAAGA` for `PCR_A_i5_A` on forward-strand instruments, its reverse complement `TCTTAAGG` on reverse-complement (NextSeq 500) chemistry |
| Read 2 | `illumina.TRUSEQ_READ2` | **Tn5 barcode (8) + T + ME (19)**, then genomic DNA from cycle 29 |

So read 2 carries the level-1 barcode in its first 8 cycles and loses 28 cycles to
barcode + ME. Upstream s3-ATAC lists the same four primer sequences and the same 8-nt
Tn5 barcode at the start of read 2 — agrees.

## 6. Checking upstream (scg_lib_structs) against the paper

- s3-WGS page: says the chemistry is s3-ATAC on crosslinked, nucleosome-depleted nuclei,
  with U-containing Tn5 adaptors, U-intolerant gap fill and U-tolerant PCR. **Agrees**
  with the Online Methods. It gives no conditions specific to WGS (fixation, SDS, PCR
  cycle count) — those are from the paper only.
- s3-ATAC page (the drawing s3-WGS defers to): every oligo sequence (Tn5 adaptor, ME
  bottom strand, A14_ME_LNA with its LNA and inverted-dT marks, i5/i7 PCR primer layouts)
  **matches** Supp. Tables 2–5. Its final library equals the one assembled in §4. Its four
  sequencing primers are not in the supplementary tables; they equal the stock
  `nextera.READ1_PRIMER`, `illumina.INDEX1_PRIMER`, `nextera.INDEX2_PRIMER` and
  `illumina.TRUSEQ_READ2` in `lib/` 🟡.
- Minor difference: upstream's step 2 reads as "sort nuclei into wells, then tagment";
  in the paper ~thousands of nuclei are **pipetted** into the tagmentation plate, and
  only after tagmentation are they pooled and **flow-sorted** into the PCR plate.
  Upstream's step 3 ("pool, redistribute, gap fill") matches the paper.
- Upstream adds the "only 50 % of a regular Tn5 library is amplifiable" framing; the
  paper makes the same argument. 🟢 both.

## 7. Open questions

- 🔴 The protocols.io stepwise protocol (beb3jaqn) was not read: per-well volumes for
  Triton X-100, number of nuclei sorted per PCR well for s3-WGS, and plate layouts beyond
  Fig. 4b.
- 🔴 How the 96 transposomes are loaded (annealing of `SBS12_18_UME_sci_N` with
  `Nextera_Mosaic_End_REVCOMP`, Tn5 source) — deferred to ref. 11.
- 🟡 How completely adapter switching converts both strands is argued (multiple cycles,
  LNA) but the efficiency is not reported as a number in the methods.
- 🟡 The ME/ME′ snap-back competition that the LNA is said to overcome is the paper's
  rationale; the 10-cycle program and 59 °C anneal are not otherwise justified.
- 🔴 s3-GCC (AluI digestion 2 h 37 °C, T4 ligation 16 °C overnight, then the same s3
  workflow) shares this chemistry; not covered by a separate note.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC manuscript and the Nature supplements, but ran
past a 10-minute timeout without writing `MANIFEST.tsv` (stuck converting the
18–39 MB source-data workbooks); the oligo workbook `MOESM6` had no `.txt` twin and was
converted by running `tools/doctext.py` on it alone. The NIHMS supplement links returned
HTML error pages saved under `.pdf`/`.xlsx`/`.zip` names. `tools/scrape_primers.py` on the
whole directory timed out (same large workbooks); on `MOESM6_ESM.xlsx.txt` alone it found
the Tn5 adaptor, the dU and inverted-dT modifications and recognised `TRUSEQ_READ2`,
`ME` and `ME_RC`. Reaction order and conditions were read from the methods by hand.
