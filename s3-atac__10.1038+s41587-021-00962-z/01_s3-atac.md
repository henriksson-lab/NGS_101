# s3-ATAC — symmetrical strand sci-ATAC by uracil-based adapter switching

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred (including
> the upstream scg_lib_structs page, a secondary source) · 🔴 not published / not
> available. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`).

**s3-ATAC** (and s3-WGS, s3-GCC from the same paper) — Mulqueen RM, Pokholok D,
O'Connell BL, Thornton CA, Zhang F, O'Roak BJ, Link J, Yardımcı GG, Sears RC, Steemers FJ,
Adey AC. "High-content single-cell combinatorial indexing." *Nature Biotechnology* 39,
1574–1580 (2021). doi:[10.1038/s41587-021-00962-z](https://doi.org/10.1038/s41587-021-00962-z)
· PMID 34226710 · PMC8678206 (author manuscript). Data: GEO GSE174226.
Step-by-step protocols on protocols.io: s3-ATAC
[10.17504/protocols.io.bd6wi9fe](https://dx.doi.org/10.17504/protocols.io.bd6wi9fe),
s3-WGS [10.17504/protocols.io.beb3jaqn](https://dx.doi.org/10.17504/protocols.io.beb3jaqn),
s3-GCC [10.17504/protocols.io.beb4jaqw](https://dx.doi.org/10.17504/protocols.io.beb4jaqw).
The catalogue (`catalogue/scg_lib_structs.tsv`) lists only this one paper for s3-ATAC.

Sources read (in `$CHEM_DATA/sources/s3-atac__10.1038+s41587-021-00962-z/`, never
committed):

| File | What | Used for |
|---|---|---|
| `s41587-021-00962-z_PMC8678206.html.txt` | PMC author manuscript, full text | rationale, **Online Methods** (all conditions) |
| `s41587-021-00962-z_41587_2021_962_MOESM6_ESM.xlsx.txt` | Supplementary Tables 1–8 (Springer copy) | **every oligo** (Tables 2–5), the stepwise molecular layout (Table 1), cellID definition (Table 6 info) |
| `upstream_s3-ATAC.html.txt` | Teichmann-lab scg_lib_structs page | second, secondary account — checked against Tables 1–5 below |
| `..._MOESM1_ESM.pdf.txt` | Supplementary Information (4 pages) | figures only (marker tracks, down-sampling); no chemistry |
| `..._MOESM2_ESM.pdf` | 4-page "Print To PDF" document, no text layer | not read (presumably the Reporting Summary) 🟡 |
| `..._MOESM3/5_ESM.zip`, `MOESM4,7–13_ESM.xlsx` | peak BED, s3-GCC compartment plots, DE tables, source data | not chemistry; not read |

Could not be fetched:

| What | URL | Why it matters |
|---|---|---|
| protocols.io s3-ATAC protocol | https://dx.doi.org/10.17504/protocols.io.bd6wi9fe | transposome assembly (annealing) conditions, nuclei per PCR well, exact volumes — 🔴 below |
| NIHMS supplement copies | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8678206/ (supplementary files) | reCAPTCHA page saved instead of the file; the Springer MOESM copies replace them |

---

## 1. What it is

A **sci-ATAC** variant (two-level combinatorial indexing: indexed Tn5 per tagmentation
well, then indexed PCR per sorted well) in which the transposome is loaded with **one
adapter species only**, and the second adapter is added afterwards by **adapter
switching**. Every tagmentation fragment then carries the same adapter at both ends and,
after switching, a different one at each end — so both strands of every fragment become
a library molecule. 🟢 (main text; Fig. 1)

The problem it solves 🟢: with an s5/s7 Tn5 mix, only fragments with two *different*
adapters amplify, which is about half of them (see
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)). The paper reports 6–13× more
usable unique reads per cell than snATAC, 10x scATAC and dscATAC on mouse whole brain.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Single-adapter, indexed Tn5** carrying TruSeq-Read-2 sequence + 8-nt index + a **dU** before the ME | sci-ATAC loads s5/s7 Nextera adapters with indices and needs custom sequencing primers |
| 2 | **Uracil as a polymerase stop**: gap fill with a uracil-*intolerant* polymerase (NPM) stops at the dU, so the fill-in copies only the ME, not the index or the adapter | dU-containing adapters (USER-cleavable) in NEBNext; here the dU is never excised, only read through later |
| 3 | **LNA, 3'-blocked switching oligo** (s5 + ME with LNA bases, 3' inverted dT) that anneals to the copied ME' and templates s5' onto the 3' end of every strand, over 10 thermocycles | splint/template-switch ideas; compare [template switching](../ref/concepts/template-switching.md) on RNA |
| 4 | PCR with a uracil-*tolerant* polymerase (Q5U) reads through the dU | — |
| 5 | Mixed TruSeq / Nextera adapters chosen so that **standard Illumina sequencing primers** work (Nextera Read 1 + Index 2, TruSeq Read 2 + Index 1) | most sci- methods need custom read/index primers |

## 2. Oligos

🟢 Verbatim from Supplementary Tables 2–5 (`MOESM6_ESM.xlsx`). The tables split each oligo
into named columns; here the columns are joined, separated by spaces. IDT notation as
given: `/ideoxyU/` = internal 2'-deoxyuridine, `/5Phos/` = 5' phosphate, `+N` = locked
nucleic acid, `/3InvdT/` = 3' inverted dT (Table 3 defines the last two).

```
Tn5 adapter, 96 wells (Table 2; Truseq_R2_SBS12_partial · Tn5_96plx_idx · U-ME)
SBS12_18_UME_sci_1    CGTGTGCTCTTCCGATCT GAACCGCG /ideoxyU/AGATGTGTATAAGAGACAG
 ...                  (idx for _1.._96, wells A1..H12; _96 = ATTGTGAA)

Annealed bottom strand (Table 2)
Nextera_Mosaic_End_REVCOMP   /5Phos/CTGTCTCTTATACACATCT

Adapter-switching oligo (Table 3; Nextera_R1_A14 · U-ME)
A14_ME LNA            TCGTCGGCAGCGTC AGATGTGTA+TA+AG+AG+AC+AG/3InvdT/

i7 PCR primers, 32 (Table 4; i7_Flowcell_Primer · i7_pcr_idx · Truseq_R2_SBS12_full)
PCR_i7_P7.S701        CAAGCAGAAGACGGCATACGAGAT TCGCCTTA GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
 ...                  S701..S708, T267..T278_IPE2F, T179..T190_IPE2F (last: GCTGGCAT)

i5 PCR primers, 64 (Table 5; i5_Flowcell_Primer · i5_pcr_idx · i5_Nextera_A14_partial)
PCR_A_i5_A            AATGATACGGCGACCACCGAGATCTACAC CCTTAAGA TCGTCGGCAGCGTC
 ...                  PCR_A..H_i5_A..H (last, PCR_H_i5_H: TGCGCTGA)
```

Table 5 is titled "Oligonucleotides used for PCR, i7 Molecule side" although it lists the
i5 primers — a copy-paste slip in the source. 🟢 (title) / 🟡 (reading)

### How they interlock — 🟡 (computed with `lib/`)

- **Tn5 adapter** = `illumina.TRUSEQ_READ2[-18:]` (the 3' 18 nt of the TruSeq Read 2
  primer, "SBS12 partial") + 8-nt Tn5 index + dU + **`nextera.ME`** (all 19 nt) — 46 nt.
  The dU sits immediately 5' of the ME: seen from a polymerase extending a genomic 3' end
  across the ME', it is the first base beyond the ME — which is what the paper means by
  "a uracil base immediately following the mosaic end".
- **ME bottom** = `nextera.ME_RC`, 5'-phosphorylated: the usual non-transferred strand.
- **A14_ME LNA** bases = **`nextera.ADAPTOR_S5`** exactly (s5 + ME, 33 nt), plus the
  3' inverted dT. The five LNA bases are at ME positions **10, 12, 14, 16, 18** (T, A, A,
  A, A): alternate bases of the 3' half of the ME. It is the *same sense* as the Tn5
  adapter's ME, so it pairs with the **ME'** (`ME_RC`) that gap fill writes on the 3' end of
  each strand, and its 5' s5 is a single-stranded template overhang.
- **i7 primers** = `illumina.P7` + 8-nt index + **`illumina.TRUSEQ_READ2`** (full 34 nt) —
  66 nt, all 32. The primer's 3' 18 nt *are* the Tn5 adapter's SBS12-partial, so the i7
  primer anneals to the complement of the Tn5 adapter and extends the read-2 handle to full
  length.
- **i5 primers** = `illumina.P5` + 8-nt index + **`nextera.S5`** (14 nt) — 51 nt, all 64.
  They anneal to the s5' written by adapter switching.
- Index orientation: the i7 is written as the **reverse complement** of the i7 read
  (`S701` contains `TCGCCTTA`, read as `TAAGGCGA` = Illumina N701). Index sets: 96 Tn5
  indices (all distinct, minimum pairwise Hamming distance 4), 32 i7 and 64 i5 (distinct,
  minimum Hamming distance 2 within each set).
- Cell barcode = **i5 + i7 + Tn5 index** — 🟢 (Table 6 info: "cellID ... i5_PCR_idx +
  i7_PCR_idx + Tn5_idx"). Maximum index space 96 × 32 × 64 🟡; the experiments used far fewer.

## 3. Step by step

Conditions 🟢 from Online Methods ("s3-ATAC Library Generation") unless marked.

1. **Assemble 96 indexed transposomes** "using previously-described methods" (Amini et al.
   2014): each Tn5 adapter annealed to the 5'-phosphorylated ME bottom strand and loaded
   onto Tn5. Diluted to **2.5 µM** in 50 % glycerol, 100 mM NaCl, 50 mM Tris pH 7.5,
   0.1 mM EDTA, 1 mM DTT; stored at −20 °C. Annealing and loading conditions 🔴 (not in
   the paper; probably on protocols.io).
2. **Nuclei**: frozen mouse whole brain or human cortex, minced on dry ice, dounced in
   NIB-HEPES (10 mM HEPES-KOH pH 7.2, 10 mM NaCl, 3 mM MgCl₂, 0.1 % IGEPAL CA-630,
   0.1 % Tween, protease inhibitor), 35 µm filter, two washes, diluted to 1,400 nuclei/µL.
3. **Tagment**: 420 µL nuclei + 540 µL 2× TD buffer (Nextera XT); **8 µL (~5,000 nuclei)
   per well + 1 µL of 2.5 µM indexed transposome**; 55 °C 10 min, then ice. Every fragment
   is now `adapter–gDNA–adapter`, the same Tn5 adapter at both ends, each strand carrying
   the transferred adapter at its 5' end, with the 9-nt gap and the annealed ME bottom
   strand opposite. 🟡 (mechanism: [Tn5](../ref/concepts/tn5-tagmentation.md))
4. **Pool, DAPI** (2 µL of 5 mg/mL), **flow sort** (Sony SH800) into 96-well plates with
   9 µL 1× TD buffer at 4 °C; gating on size, complexity, DAPI. Nuclei per PCR well: shown
   only in the Fig. 2a plate layout 🔴 (not stated in text).
5. **Denature** nucleosomes and remaining Tn5: 1 µL 0.1 % SDS (~0.01 % final).
6. **Gap fill with a uracil-intolerant polymerase**: 4 µL NPM (Nextera XT), 72 °C 10 min.
   Each genomic 3' end is extended across the ME (writing ME') and stops at the dU, so the
   Tn5 index and SBS12 are *not* copied. Each strand is now
   `5'-SBS12p · idx · dU · ME · gDNA · ME' -3'`. 🟢 (Table 1: "Use Uracil to block any pcr
   amplification, just want to close the gap"; segment list 🟡)
7. **Adapter switching**: add 1.5 µL of 1 µM A14-LNA-ME; 98 °C 30 s, then **10 cycles of
   98 °C 10 s, 59 °C 20 s, 72 °C 10 s**; hold 10 °C. In each cycle the LNA oligo's ME pairs
   with the 3' ME' of a strand (LNA raises its Tm so it outcompetes the strand's own
   partner), and NPM extends the strand over the oligo's 5' overhang, appending s5'
   (`GACGCTGCCGACGA`). The inverted dT stops the oligo itself from being extended. Result,
   for both strands: `5'-SBS12p · idx · dU · ME · gDNA · ME' · s5' -3'`. 🟢 conditions and
   rationale (main text, Table 1) / 🟡 segment list
8. **Quench SDS** with 1 % Triton X-100. Plates can be stored at −20 °C for weeks
   (Table 6 records fresh, 1-week and 10-week plates).
9. **Indexed PCR with a uracil-tolerant polymerase**: 16.5 µL sample + 2.5 µL i7 primer
   (10 µM) + 2.5 µL i5 primer (10 µM) + 3 µL water + 25 µL NEBNext **Q5U** 2× master mix +
   0.5 µL 100× SYBR Green I (50 µL). 98 °C 30 s; **16–18 cycles** of 98 °C 10 s, 55 °C 20 s,
   72 °C 30 s, read, 72 °C 10 s; stopped at inflection (qPCR on a CFX). The i5 primer
   primes on s5' and its extension **reads through the dU** (an A is put opposite it);
   the i7 primer then primes on the complement of SBS12-partial. 🟢 (Table 1: "Nextera i5
   Primes first and extends through Uracil") A "Half_Vol" PCR variant exists (Table 6);
   its volumes are not given in the text read 🔴.
10. **Clean-up**: pool 25 µL/well, QIAquick column, 1× SPRI (Omega Mag-Bind TotalPure),
    elute 31 µL; Qubit, TapeStation D5000, dilute to 1 nM on the 100–1,000 bp range.

### Variants in the same paper (s3-WGS, s3-GCC) — same oligos

- **s3-WGS**: nuclei fixed in 0.75 % formaldehyde 10 min, quenched in NIB-Tris, then
  **nucleosome depletion**: 40 µL 1 % SDS added to nuclei in 760 µL 1× NEBuffer 2.1,
  37 °C 20 min (0.05 % SDS final 🟡 *(computed: 40 µL × 1 % / 800 µL; the paper gives
  only the volumes)*); 500 nuclei/µL into tagmentation; everything after as s3-ATAC, PCR
  13–15 cycles. 🟢
- **s3-GCC**: the same fixed, nucleosome-depleted nuclei, **AluI** digest (100 U, 2 h
  37 °C) and **proximity ligation** with T4 ligase (16 °C, 14 h) *without* ligation
  junction enrichment, then as s3-WGS. 🟢 The library structure is that of s3-ATAC; only
  the insert differs (genome-wide, with ligation junctions). 🟡

## 4. Final library — 🟡 (computed with `lib/`)

Top strand, P5 side first, 5'→3':

```
P5 (29) · i5 (8) · s5 (14) · ME (19) · <genomic insert> · ME' (19) · A (1) · Tn5-idx' (8) · TruSeq-Read2' (34) · i7' (8) · P7' (24)
```

- `A` is the complement of the dU (it becomes T on the other strand after PCR —
  Table 1 writes the product with `/T/`). 🟢 (Table 1)
- Adapter overhead 164 nt (70 on the P5 side, 94 on the P7 side).
- Checks: the computed construct (with N for indices) is base-for-base identical to the
  upstream page's final library; and the final row of Supplementary Table 1 (written P7
  side first, with i7 = S701, Tn5 index = sci_1, i5 = `TAGATCGC`, insert omitted) is its
  exact reverse complement. 🟡 *(computed)*. Note: the example i5 in Table 1
  (`TAGATCGC`, the Illumina N501 index) is **not** among the 64 i5 primers of Table 5 —
  Table 1 is illustrative.

## 5. Sequencing

🟢 NextSeq 500 (High or Mid output, 150-cycle kits) and NovaSeq S2; paired end, **85 + 85
cycles, 10-cycle index reads**. Standard primers 🟢 (main text: TruSeq Read 2 and Index 1,
Nextera Read 1 and Index 2); which primer reads what 🟡 *(computed)*:

| Read | Primer (lib/) | Reads |
|---|---|---|
| Read 1 (85) | `nextera.READ1_PRIMER` (s5 + ME) | genomic insert from base 1 |
| Index 1 (10) | `illumina.INDEX1_PRIMER` | i7 (8; the reverse complement of the oligo's index) + 2 nt of P7' |
| Index 2 (10) | `nextera.INDEX2_PRIMER` (ME' + s5') | i5 (8) + 2 nt; orientation relative to the oligo depends on the instrument's i5 workflow (NextSeq reads it reverse-complemented) |
| Read 2 (85) | `illumina.TRUSEQ_READ2` | **Tn5 index (8)** · T (from the dU) · ME (19) · genomic insert — 28 adapter bases before the genome |

## 6. Check against the upstream scg_lib_structs page

| Point | Upstream | Paper / Tables | Verdict |
|---|---|---|---|
| Tn5 adapter, ME bottom, A14_ME_LNA, i5/i7 primers | sequences with `[index]` placeholders | Tables 2–5 | **agree** (bases and modifications) |
| Gap fill with uracil-intolerant NPM; PCR with Q5U | yes | Methods | agree |
| Read 1 / Index 1 / Index 2 / Read 2 primers | TruSeq R2 + Index 1, Nextera R1 + Index 2 | main text | agree |
| Cycles | 85 / 10 / 10 / 85 | Methods | agree |
| Final library | P5·i5·s5·ME·gDNA·ME'·A·Tn5idx·TruSeqR2'·i7·P7' | Table 1 (reverse complemented) | agree *(computed)* |
| Cell barcode = Tn5 + i7 + i5 | yes | Table 6 info | agree |
| Read 2 content | "first 8 cycles are the Tn5 barcodes" | — | agrees; upstream does not say that the next 20 cycles are T + ME |
| Adapter switching on both strands | draws the LNA oligo on both strands; "products ... are the same" | main text "all library fragments"; Table 1 draws one strand only | consistent |
| Conditions (volumes, temperatures, cycles, nuclei numbers) | not given | Methods | upstream silent |

No disagreement found.

## 7. Open questions

- 🔴 Transposome assembly (oligo annealing, Tn5 source and loading) — deferred to Amini
  et al. 2014 and protocols.io, neither read.
- 🔴 Nuclei sorted per PCR well (only in the Fig. 2a image) and the "Half_Vol" PCR recipe.
- 🔴 Whether leftover switching oligo matters in the PCR. It carries s5 (the i5 primer's
  3' end) and could otherwise act as a short forward primer; the inverted dT blocks its
  extension. 🟡 (reasoning; the paper only says the oligo is 3'-blocked)
- 🟡 Strand fate of the original ME bottom oligo after SDS + gap fill (displaced or
  ligated?) — not discussed; it does not enter the final product.
- 🟡 i5 read orientation for the NovaSeq runs (v1.0 vs v1.5 reagents) — not stated.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` timed out (killed after 10 min, no `MANIFEST.tsv` written) while
converting the very large source-data spreadsheets (MOESM4/10/12, 17–39 MB); the files
had been downloaded by an earlier run. The NIHMS supplement copies are reCAPTCHA pages
saved under `.pdf` / `.xlsx` names. The upstream page and Supplementary Tables 1–8
(`MOESM6`) were converted by hand with `tools/doctext.py`. `tools/scrape_primers.py` on the
whole directory also timed out (tens of MB of `.xlsx.txt`); run on the three relevant
`.txt` files it found every oligo segment. Because the tables split each oligo into
columns, the scraper reports segments, not whole oligos; `--find` nevertheless locates
the joined sequences across column breaks.
