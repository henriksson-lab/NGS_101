# BD Rhapsody WTA — microwell bead capture, random-primer second strand

> **Evidence marking.** 🟢 verbatim from a primary source (the paper, or a BD protocol
> document) · 🟡 derived, inferred, or taken only from the secondary scg_lib_structs page ·
> 🔴 not published / not available. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note; they are not yet asserted in a self-test, because
> this protocol has no `tools/` module yet (status `notes`).

**BD Rhapsody** (formerly BD Resolve) is a commercial BD Biosciences platform. The
catalogue keys it to its first use in a Nature Methods paper:

Yoon S-J, Elahi LS, Pașca AM, Marton RM, Gordon A, Revah O, Miura Y, Walczak EM,
Holdgate GM, Fan HC, Huguenard JR, Geschwind DH, Pașca SP. "Reliability of human cortical
organoid generation." *Nat Methods* 16, 75–78 (2019).
doi:[10.1038/s41592-018-0255-0](https://doi.org/10.1038/s41592-018-0255-0) · PMID 30573846 ·
PMC6677388.

(H. Christina Fan, a co-author, is first author of the CytoSeq paper below.) Papers
this one cites for the platform (not fetched, not in the catalogue row): Fan HC,
Fu GK, Fodor SP, "Combinatorial labeling of single cells for gene expression cytometry",
*Science* 347, 1258367 (2015), doi:10.1126/science.1258367 (the "stochastic barcoding" /
CytoSeq origin, their ref. 15); Birey F et al., *Nature* 545, 54–59 (2017) (their ref. 12,
earlier organoid data on the same platform). The catalogue lists no other papers for
this method.

**The paper contains no oligo sequences for the platform** — it is a user of a commercial
kit. The chemistry below comes from the BD protocol documents that scg_lib_structs
hosts, and from the scg_lib_structs page itself (whose author says the bead and primer
sequences are an "educational guess" from sequencing data).

### Sources read

All in `_data/sources/bd-rhapsody__10.1038+s41592-018-0255-0/` (never committed).
`upstream_data_*` were fetched by hand (curl) from the links on the upstream page, not by
`get_sources.py`.

| File | What | Used for |
|---|---|---|
| `s41592-018-0255-0_PMC6677388.html.txt` | the paper (PMC full text) | Methods "Single-cell gene expression (BD Rhapsody system)": capture, the (older) library description, sequencing |
| `s41592-018-0255-0_41592_2018_255_MOESM1_ESM.pdf.txt` | Supplementary Information (Springer) | qPCR primers only — not Rhapsody chemistry |
| `s41592-018-0255-0_41592_2018_255_MOESM2_ESM.pdf.txt` | Reporting Summary | software only |
| `s41592-018-0255-0_41592_2018_255_MOESM3_ESM.xlsx.txt` | Supplementary Table 3 (sample sheet) | nothing chemical |
| `upstream_BD_Rhapsody.html.txt` | scg_lib_structs "BD Rhapsody WTA" page | every bead/primer sequence, final library, read layout (secondary) |
| `upstream_data_GMX_BD-Rhapsody-WTA-alpha-Protocol_UG_EN.pdf.txt` | BD "Whole transcriptome analysis alpha protocol", 23-21179-00, 12/2018 | **Randomer sequence**, RPE, RPE PCR, index PCR, cleanups |
| `upstream_data_GMX_BD-Rhapsody-Single-Cell-Analysis-System-Instrument_UG_EN.pdf.txt` | BD Rhapsody instrument user guide, Doc ID 214062 Rev. 2.0, 02/2019 | lysis, bead retrieval, RT mix, Exonuclease I; bead-oligo cartoon |
| `upstream_data_BD_CLS1.txt`, `_CLS2.txt`, `_CLS3.txt` | cell-label lists hosted by scg_lib_structs | counting / overlap of the three CLS sets |

Not obtained (PMC download gate; the Springer copies above cover the same supplements
except possibly the figure DOCX):

| File | URL |
|---|---|
| Supplementary Figures (DOCX) | https://pmc.ncbi.nlm.nih.gov/articles/instance/6677388/bin/NIHMS1511424-supplement-Supplementary_Figures.docx |
| Supplementary Tables 1, 2, 4 (PDF) | https://pmc.ncbi.nlm.nih.gov/articles/instance/6677388/bin/NIHMS1511424-supplement-Supplementary_Tables_1__2__4.pdf |
| Supplementary Table 3 / Reporting Summary (PMC copies) | same `/bin/` directory |
| Fan et al. 2015 (Science; Cloudflare/paywall) | https://doi.org/10.1126/science.1258367 |
| BD's own bead-oligo / Universal Oligo / AbSeq primer sequences | 🔴 not published by BD |

---

## 1. What it is

Cells are loaded by limiting dilution into a cartridge of **>200,000 microwells**, then
barcoded **magnetic beads** are added to saturation so that each cell-containing well gets
a bead 🟢 (paper, Methods). Lysis in the well lets poly(A) RNA hybridise to the bead's
oligo-dT; beads are pulled out with a magnet and pooled, and **reverse transcription
happens on the pooled beads in one tube** 🟢. The cDNA stays covalently on the bead, so
beads can be archived and subsampled 🟢 (paper: ~67 % of beads used, rest archived; BD
alpha protocol: leftover beads stored at 4 °C up to 3 months).

The cell barcode is **split-pool combinatorial**: three 9-nt cell-label blocks (CLS1–3)
separated by fixed linkers, then an 8-nt UMI and oligo-dT 🟡 (upstream). Each of the
three lists hosted upstream has **97** distinct 9-mers, no sequence shared between lists,
so 97³ = **912,673** labels 🟡 (computed from the CLS files; matches the upstream figure).

| | Feature | Compare |
|---|---|---|
| 1 | Microwell + bead capture with a **magnet**, beads retrieved and processed in bulk | Microwell-seq, Seq-Well (also microwells; Drop-seq-type beads) |
| 2 | **Exonuclease I** after RT to clear unused bead oligo-dT | Drop-seq/Seq-Well do the same with Exo I |
| 3 | **No template switching**: the second strand is made **on the bead** by a randomer carrying a TruSeq Read 2 handle, extended by Klenow exo⁻ ("RPE") | 10x / Drop-seq use a TSO; here the cDNA 5' end is never captured |
| 4 | The library is 3'-anchored by the bead and cut at a **random** position by the randomer, so no fragmentation/tagmentation step | Tn5-based 3' kits |

Concept background: [reverse transcription](../ref/concepts/reverse-transcription.md)
(oligo-dT on a bead, random priming). Template switching
([concept](../ref/concepts/template-switching.md)) is **not** used — that is the main
difference from bead protocols that make full-length cDNA first.

### Two different library-prep descriptions — disagreement

- **The paper** (Dec 2018) says the whole-transcriptome libraries were made by
  "second strand cDNA synthesis, adaptor ligation, and universal amplification" with
  **22 PCR cycles**, followed by **random priming PCR** to enrich the 3' ends 🟢
  (paraphrased). That is a ligation-based WTA; no adaptor sequences are given 🔴.
- **The BD WTA alpha protocol** (12/2018) and the **upstream page** describe a different
  route: random-primer extension on the beads (RPE), RPE PCR (13 cycles), index PCR
  (8–9 cycles) 🟢 (BD protocol). No ligation, no 22-cycle step.

So upstream's drawn chemistry is BD's commercial WTA workflow, **not** necessarily what
Yoon et al. ran. The final library layout (bead side on Read 1, cDNA on Read 2) would
be the same in both, because the bead oligo is the same; the Read 2 adaptor of the
ligation route is unknown 🔴. The note below follows the documented RPE route.

## 2. Oligos

Sequences as written in the source. Upstream writes no modifications for any of them;
none are given anywhere else either.

**From the BD WTA alpha protocol** (third-party reagent the user orders) 🟢:

```
Randomer   5' TCA GAC GTG TGC TCT TCC GAT CTNNNNNNNNN 3'
```

(written with codon spaces in the source; 100 µM in 10 mM Tris, 0.1 mM EDTA pH 8.0; N
hand-mixed 25 % each 🟢)

**From the upstream page only** 🟡 (secondary; BD does not publish these — upstream
calls the bead sequences an educated guess from sequencing data):

```
V1 bead oligo        |--5'-CCCCCCTCTCTCTCTACACGACGCTCTTCCGATCT[CLS1]ACTGGCCTGCGA[CLS2]GGTAGCGGTGACA[CLS3][8-bp UMI](T)18 -3'
Enhanced bead (2022) |--5'-CCCCCCTCTCTCTCTACACGACGCTCTTCCGATCT[VB][CLS1]GTGA[CLS2]GACA[CLS3][8-bp UMI](T)18 -3'

Pre-amp forward      5'- GACGCTCTTCCGATCT -3'
Pre-amp reverse      5'- TCAGACGTGTGCTCTT -3'
Library Forward      5'- AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'
Library Reverse      5'- CAAGCAGAAGACGGCATACGAGAT[8-bp sample index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -3'

TruSeq Read 1        5'- ACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'
TruSeq Read 2        5'- GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -3'
Index read primer    5'- GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -3'
```

`|--` = bead surface. `[VB]` = "variable bases", none / A / GT / TCA, Enhanced beads only.
CLS1/2/3: 9 nt each, 97 per position (`upstream_data_BD_CLS*.txt`).

The BD documents name the PCR reagents only: **"Universal Oligo"** and **"BD AbSeq
Primer"** (RPE PCR), **"Library Forward Primer"** and **"Library Reverse Primer (1–4)"**
(index PCR) 🟢. Their sequences are not given 🔴. The instrument guide's cartoon of the
capture oligo reads `bead – Univ – CL – UMI – poly(T)` ("Univ: universal oligo; CL: cell
label") 🟢 — the same order as upstream's bead oligo.

## 3. How they interlock — 🟡 (computed with `lib/illumina`)

- **Bead oligo 5' part** (35 nt) = 13-nt `CCCCCCTCTCTCT` + the last **22 nt of
  `illumina.TRUSEQ_READ1`** (`CTACACGACGCTCTTCCGATCT`). The 13-nt C/CT head has no
  role in the final library (it is outside every primer) 🟡; its purpose is not stated 🔴.
- **Randomer** = the last **23 nt of `illumina.TRUSEQ_READ2`** + N₉ (32 nt). So every
  random-primed second strand starts with a partial Read 2 handle.
- **Pre-amp forward** (16 nt) = the last 16 nt of `TRUSEQ_READ1`, and the 3' end of the
  bead-oligo head. **Pre-amp reverse** (16 nt) = the first 16 nt of the randomer (a
  sub-sequence of `TRUSEQ_READ2`). A 16-nt pair on the two handles — the expected shape
  of BD's "Universal Oligo", which upstream *assumes* is this pair 🟡; BD's kit lists one
  "Universal Oligo" plus an AbSeq primer, sequences undisclosed 🔴.
- **Library Forward** (58 nt) = `illumina.P5` + `TRUSEQ_READ1` overlapping by 4 nt
  (`ACAC`) — identical to `illumina.NEBNEXT_UNIVERSAL_PRIMER` / `TRUSEQ_P5_FULL`. It
  restores the full Read 1 site that the 16-nt pre-amp primer had trimmed.
- **Library Reverse** (66 nt) = `illumina.P7` + 8-nt i7 + full `TRUSEQ_READ2` (34 nt).
  Its 3' 23 nt are the randomer's fixed part, so it primes on the RPE PCR product.
  Orientation of the 8-nt index as written relative to the i7 read: not checkable — no
  index sequences are given 🔴 (BD's "Library Reverse Primer 1–4" = four indices).
- **Index read primer** = reverse complement of `TRUSEQ_READ2` minus its 5'-most base,
  i.e. `illumina.INDEX1_PRIMER` — the standard TruSeq i7 read.
- **Linkers**: V1 linker 1 `ACTGGCCTGCGA` is 12 nt, linker 2 `GGTAGCGGTGACA` 13 nt.
  The Enhanced-bead linkers `GTGA` and `GACA` are both **sub-sequences of V1 linker 2**
  (positions 8–11 and the last 4 nt); `GTGA` does not occur in linker 1.
- **Read 1 length needed**: V1 CLS1 + L1 + CLS2 + L2 + CLS3 + UMI = 9+12+9+13+9+8 =
  **60 nt** — upstream's "at least 60 cycles". Enhanced: 43 nt plus 0–3 nt VB = 43–46.

## 4. Step by step

Steps 1–4 from the BD instrument user guide (Doc 214062 Rev. 2.0), steps 5–11 from the
BD WTA alpha protocol; conditions 🟢 unless marked. Strand bookkeeping 🟡.

1. **Load cells** (~10,000 in the paper 🟢; 5,000–10,000 per sample recommended for WTA
   🟢) into the cartridge microwells by limiting dilution, then load **Cell Capture
   Beads** to saturation, image with the scanner.
2. **Lyse**: 550 µL Lysis Buffer + DTT (75 µL 1 M DTT per 15 mL), room temperature
   **2 min**. Poly(A) RNA anneals to the bead `(T)18` 🟡 (the paper says only that poly(A)
   RNA hybridises to the beads).
3. **Retrieve beads** magnetically into a 5 mL tube with 4,950 µL lysis buffer; wash.
4. **RT on beads**, 200 µL cDNA mix (RT buffer, dNTP, 0.1 M DTT, Bead RT/PCR Enhancer,
   RNase inhibitor, reverse transcriptase), **37 °C 20 min, 1,200 rpm** ("shaking is
   critical"). First strand: `|--bead oligo · (T)18 · cDNA →`, the cDNA covalently on the
   bead with cell label and UMI at its 5' end 🟢 (paper: tagged on the cDNA 5' end =
   mRNA 3' end).
5. **Exonuclease I**, 200 µL, **37 °C 30 min**, then **80 °C 20 min** to inactivate.
   Removes bead oligos that captured nothing (single-stranded 3' ends) 🟡 (purpose not
   stated in the BD text).
6. **Remove the RNA**: beads 95 °C 2 min, magnet, discard supernatant. The bead now holds
   single-stranded cDNA 🟡.
7. **Random priming**: 174 µL Random Primer Mix (water, NEBuffer 2, 20 µL Randomer);
   95 °C 2 min (no shaking) → 37 °C 5 min → 25 °C 5 min (both 1,200 rpm). N₉ anneals at random positions
   along the bead-bound cDNA.
8. **Extension (RPE)**: + 26 µL (dNTP, Bead RT/PCR Enhancer, **Klenow exo⁻**), 37 °C
   **30 min**, 1,200 rpm. The randomer extends toward the bead and copies the bead oligo,
   so the new strand is `5'- R2 handle(23) · N₉-insert … · (A)18 · UMI' · CLS3' … CLS1' ·
   R1-handle' · AGAGAGAGGGGGG -3'` (complement written 5'→3') 🟡. Only extensions that
   reach the bead oligo carry the Read 1 handle.
9. **Release**: wash, resuspend in 200 µL Tris-Tween20, **95 °C 2 min**, magnet, keep the
   **supernatant = RPE product** (beads can be stored). Clean up with **1.6× SPRIselect**
   (320 µL to 200 µL), elute 30 µL.
10. **RPE PCR**: 50 µL PCR MasterMix + 10 µL Universal Oligo + 10 µL BD AbSeq Primer +
    30 µL RPE product, split in two; 95 °C 3 min; **13 ×** (95 °C 30 s, 60 °C 1 min,
    72 °C 1 min); 72 °C 2 min. 1× SPRI. Expected broad 200–2,000 bp; quantify the
    150–600 bp region. Upstream's drawn product: `GACGCTCTTCCGATCT · CLS1 · L1 · CLS2 ·
    L2 · CLS3 · UMI · (dT) · insert · AGATCGGAAGAGCACACGTCTGA` (top strand) 🟡 — the bead
    head's first 19 nt are lost because the forward primer is only 16 nt.
11. **WTA Index PCR**: dilute to 2 nM (150–600 bp peak); 10 µL + Library Forward +
    one Library Reverse (1–4) in PCR MasterMix (50 µL); 95 °C 3 min; **9 cycles (1–2 nM)
    or 8 cycles (>2 nM)** of 95 °C 30 s / 60 °C 30 s / 72 °C 30 s; 72 °C 1 min.
12. **Double-sided SPRI**: + 60 µL water (110 µL), take 100 µL, + 60 µL beads (0.6×, keep
    supernatant), + 10 µL beads to the supernatant (0.7× cumulative), keep beads, elute
    30 µL 🟢 volumes / 🟡 ratios. Final libraries **~300–700 bp**, >1 ng/µL 🟢.

## 5. Final library — 🟡 (upstream, checked by computation)

V1 beads, top strand 5'→3' (segment list):

```
5'- P5 (29) · TruSeq Read 1 (33, sharing ACAC with P5 → 58 nt together)
    · CLS1 (9) · linker 1 ACTGGCCTGCGA (12) · CLS2 (9) · linker 2 GGTAGCGGTGACA (13)
    · CLS3 (9) · UMI (8) · (dT)18 · <cDNA insert, antisense to mRNA>
    · revcomp(TruSeq Read 2) (34) · i7 (8) · P7' (24) -3'
```

Enhanced beads: the same, with `[VB] (0–3) · CLS1 · GTGA · CLS2 · GACA · CLS3 · UMI`.

Checks 🟡 (computed): upstream's written V1 library equals `illumina.TRUSEQ_P5_FULL` +
the barcode block above, and its 3' tail equals `revcomp(TRUSEQ_READ2)` + N₈ +
`illumina.P7_RC`. Fixed parts excluding the insert: **202 nt** (V1), which with the
upstream 300–700 bp final size leaves ~100–500 bp of cDNA.

The **Read 2 end of the insert is wherever the randomer landed**; the Read 1 end is the
poly(A) junction. Because the top strand is cDNA (antisense), **Read 2 reads the mRNA
sense strand** 🟡 (strand bookkeeping).

## 6. Sequencing

- 🟢 Paper: HiSeq 2500, **101 × 2**; ~38,000 reads/cell; molecular index coverage 3.2–3.5.
- 🟢 BD alpha protocol: NextSeq / MiniSeq High or Mid output, load 1–1.2 pM with **20 %
  PhiX**; 10,000 / 50,000 / 100,000 reads per cell for shallow / moderate / deep.
- 🟡 Upstream layout: **Read 1** from the TruSeq Read 1 primer, ≥60 cycles: CLS1 · L1 ·
  CLS2 · L2 · CLS3 · UMI (then poly(T)). **i7** 8 cycles from the TruSeq index primer
  (`illumina.INDEX1_PRIMER`). **Read 2** from the TruSeq Read 2 primer: cDNA. No i5 read
  (Library Forward has no index). Standard Illumina primers, no custom ones.
- With 101-nt Read 1 (paper), cycles 61–101 read into the `(T)18` and then poly(A)-
  adjacent cDNA 🟡.

## 7. Open questions

- 🔴 **What the paper actually ran**: its Methods describe ligation-based WTA
  (second strand, adaptor ligation, 22-cycle universal PCR, then random-priming PCR),
  which matches neither BD's 12/2018 RPE alpha protocol nor the upstream drawing.
  The adaptor and primer sequences of that route are unpublished.
- 🔴 The bead-oligo sequence itself (head, linkers, `(T)18` length, any anchor) is BD
  proprietary; upstream's version is inferred from reads. The 13-nt `CCCCCCTCTCTCT` head
  and the purpose of the "BD AbSeq Primer" in a WTA-only RPE PCR are unexplained.
- 🔴 "Universal Oligo" sequence: upstream assumes the 16-nt pre-amp pair; BD gives only
  the name (one tube, so either a single oligo or a premixed pair).
- 🟡 Strand displacement: Klenow exo⁻ has limited strand-displacement activity, so when
  several randomers prime the same cDNA, only bead-proximal products may reach the
  Read 1 handle — this would bias inserts toward the 3' end, consistent with the
  "3' WTA" description, but nobody states it.
- 🟡 Whether the 97-per-position CLS lists hosted upstream are complete and identical
  between V1 and Enhanced beads (upstream says the barcodes are the same).

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC text, JATS XML and the three Springer
supplements; the four PMC `/bin/` supplements stopped at the download gate (manual), and
the Europe PMC supplementary zip held only a placeholder GIF entry. None of these
contain Rhapsody chemistry. The useful primary documents (BD alpha WTA protocol,
instrument user guide, CLS lists) are linked from the upstream page but are **not
followed by `get_sources.py`**; they were fetched by hand with curl into the sources
directory and converted with `tools/doctext.py`. `tools/scrape_primers.py` found every
upstream sequence, but in listing mode it **missed the Randomer in the BD alpha
protocol** because it is written in codon-spaced triplets (`TCA GAC GTG …`); `--find`
does locate it there.
