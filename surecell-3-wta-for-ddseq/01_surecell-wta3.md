# SureCell WTA 3′ (Illumina / Bio-Rad ddSEQ): droplet 3′ RNA-seq with split-pool bead barcodes

> **Evidence marking.** 🟢 verbatim from a primary source (here: Illumina's own kit
> documents) · 🟡 derived or inferred · 🔴 not published. Relationships marked 🟡
> *(computed)* were worked out with `lib/`; the model is `tools/surecell.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence. **The checked kit guide does not print bead, TPP1 or Sequencing Primer SP bases** 🔴. Official MiSeq and later NextSeq custom-primer guides have been read for instrument loading, but do not give SP bases; NovaSeq, HiSeq and MiniSeq guides remain unchecked. `tools/surecell.py` therefore contains no unverified kit base: every kit-derived region is a bracketed, inferred placeholder. The one base-level
> fact that does have a vendor source is the N7xx index set, printed in Illumina's own
> *Adapter Sequences* document. Earlier revisions of this note printed the drawn page's
> reconstruction of the bead oligo, TPP1, the Read-1 primer and the spacers; those
> sequences are gone from it.

**Illumina Bio-Rad SureCell WTA 3′ Library Prep Kit for the ddSEQ System** — a vendor
kit (Illumina cat. 20014279 / 20014280), run on the Bio-Rad **ddSEQ Single-Cell
Isolator**. There is **no defining publication**; the reference is the kit
documentation: *SureCell WTA 3′ Library Prep Reference Guide*, Illumina document
#1000000021452 v01 (June 2017; v00 February 2017). Analysis was the BaseSpace
"SureCell RNA Single-Cell" app. The kit is **discontinued**; Bio-Rad's successor is the
*ddSEQ Single-Cell 3′ RNA-Seq Kit*, which is a different chemistry and is not covered
here.

Sources read (`/husky/henriksson/NGS101_papers/sources/surecell-3-wta-for-ddseq/`, never committed). Only the
first was fetched by `tools/get_sources.py`; the others were fetched by hand from the
Illumina support page and the `vals/umis` GitHub repository and added to `MANIFEST.tsv`:

| File | What | Used for |
|---|---|---|
| `upstream_SureCell.html.txt`, `scg_adapters_SureCell.fa` | scg_lib_structs method page, and that repository's own derived adapter FASTA (`chemistries/adapters/SureCell.fa` — not a mirrored vendor file) | **inspiration only, not sources** — nothing is cited from or verified against them |
| `scg_10000167449.pdf.txt`, `scg_Bio-Rad_ddSEQ_3p_scRNA_bead_structures.xlsx.txt` | original Bio-Rad ddSEQ 3′ scRNA kit guide and a separate scg-hosted bead-structure workbook of uncertain authorship | the **successor** kit has a different dual-index/deconvolution workflow. Its guide is source context; the workbook is not established as a vendor original and neither proves SureCell bases |
| `umis_SureCell_K562_R1.fastq` | 10 real SureCell Read 1s | **empirical** evidence for the read-1 anatomy and the two 15-nt spacers (§2) |
| `surecell-wta3-library-prep-reference-guide-1000000021452-01.pdf.txt` | Illumina Reference Guide v01 | **all steps and conditions**, kit contents, index set, sequencing primer |
| `surecell-wta3-library-prep-checklist-1000000021454-00.pdf.txt` | Illumina checklist | same steps, condensed (agrees with the guide) |
| `surecell-wta-3-pbmc-demonstrated-protocol_1000000044179-00.pdf.txt` | PBMC demonstrated protocol | variant: PBMC input, half tagmentation mix |
| `surecell-wta3-nuclei-demonstrated-protocol-1000000044178-00.pdf.txt` | nuclei demonstrated protocol | variant: nuclei instead of cells |
| `scg_illumina-adapter-sequences-1000000002694-14.pdf.txt` | Illumina *Adapter Sequences* v14, Nextera DNA Indexes table | exact N7xx index bases in adapter and sample-sheet orientation; does not disclose SureCell bead oligos |
| `umis_SureCell_transform.json` | `umis` read-1 parsing regex | independent check of the bead-oligo layout |
| `umis_SureCell_barcodes.txt` | 96 six-mers | the bead barcode whitelist |
| `umis_SureCell_K562_R1.fastq` | 10 real SureCell read 1s (K562) | read-1 length, the variable 5′ offset (§5) |

Official MiSeq Custom Primers Guide 15041638 v01 (pp. 2–5) and NextSeq 500/550 Guide
15057456 v03 (pp. 2–5) were read online: both describe loading custom primers, without
printing the SureCell SP bases. Their local PDF fetches failed DNS. Not yet checked 🔴:
NovaSeq, HiSeq and MiniSeq guides; the Bio-Rad *ddSEQ Single-Cell Isolator Instruction Manual*
(#10000069430); Illumina's ddSEQ technical / application notes. No vendor document
seen gives the proprietary SureCell bead, TPP1 or SP sequence.

---

## 1. What it is

A **droplet** single-cell 3′ RNA-seq kit. The ddSEQ cartridge co-encapsulates cells and
**barcode beads** in sub-nanolitre droplets (4 chambers per cartridge, ~300 cells per
chamber, 1200 per cartridge) 🟢 (Reference Guide, Introduction). Lysis and
**reverse transcription happen inside the droplet** from a bead-borne oligo-dT primer
carrying a cell barcode and a UMI; the emulsion is then broken and everything after
that is bulk: second strand by nick translation, **tagmentation** of the
double-stranded cDNA, and a PCR that keeps only the fragments carrying the bead end
🟢 ("How does the SureCell WTA 3′ Assay Work?", Appendix A).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | Cell barcode built **split-pool on beads**: three 6-nt blocks separated by two fixed 15-nt spacers, 96 × 96 × 96 = 884,736 combinations 🟡 | inDrop (two blocks, one spacer), SPLiT-seq / sci-RNA-seq (split-pool in cells, not on beads) |
| 2 | **No template switching**: the 5′ end of the cDNA is never needed, because the library is cut by Tn5 and only the bead (3′) end is amplified | inDrop and CEL-seq use IVT instead; 10x and Drop-seq use a TSO |
| 3 | **One-sided Tn5**: a P5 primer that reads into the bead handle plus an s7-only Nextera i7 primer — so no s5 is needed, and the read-1 primer is custom | Drop-seq / Nextera XT 3′ libraries (custom read 1 on the bead handle) |

Concept pages: [reverse transcription](../ref/concepts/reverse-transcription.md),
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). Template switching is
deliberately **not** used (see [template switching](../ref/concepts/template-switching.md)
for what it would have provided).

## 2. Oligos — 🔴 Illumina published none

The Illumina documents name the reagents — **3′ Barcode Mix**, **TPP1** (Tagment PCR
Primer), **DNA Adapters N7xx**, **Sequencing Primer SP** — and nothing else: a
`grep -E '[ACGT]{14,}'` over the extracted text of the Reference Guide, the checklist and
both demonstrated protocols returns **no hit at all**. The guide's only statement about
the custom primer is

> "A custom sequencing primer is provided for Read 1. Dilute the provided primer to the
> concentration …" — Reference Guide 1000000021452 v01 🟢

So: the bead-oligo 5′ handle, TPP1, the Read-1 primer SP, the transposome's adapter load,
the poly(dT) length and chemical modifications remain **unverified in the checked guides** 🔴, and
`tools/surecell.py` asserts none of them — each is a bracketed placeholder created by its
`inferred()` helper, which refuses anything that is not a bracketed description.

### The one vendor-published base-level fact 🟢

**DNA Adapters** supplied: N701–N707, N718 (8-sample kit); N701–N707, N710–N712,
N714–N716, N718–N724, N726–N729 (24-sample kit) 🟢 (Reference Guide, Kit Options).
Illumina *Adapter Sequences* v14 prints those N7xx bases in both orientations 🟢 — e.g.
N701 is `TCGCCTTA` in the adapter and `TAAGGCGA` in the sample sheet. That identifies the
index window only; the full SureCell index-PCR oligo around it stays 🔴.

### What the real reads establish — 🟡 empirical, not a vendor source

`umis_SureCell_K562_R1.fastq` holds ten real SureCell Read 1s, and `vals/umis` ships a
parsing regex for them. Both are third-party, not vendor documents, so nothing here is
🟢; but the evidence is sequencing data rather than someone's drawing:

```
read 1:  [0-? offset][6 nt]TAGCCATCGCATTGC[6 nt]TACCTCTGAGCTGAA[6 nt]ACG[8 nt]GAC(T)n…
```

- The two 15-nt constant spacers occur **18** and **14** times across the ten records
  (sequence plus header copy), i.e. in every parseable read 🟡.
- **Bead barcodes**: `umis_SureCell_barcodes.txt` lists 96 distinct 6-mers 🟡; in the nine
  parseable reads every barcode1/2/3 matches exactly one of them (N as wildcard), so one
  96-list serves all three positions 🟡 *(computed; 9 reads)*. Minimum pairwise Hamming
  distance within the list is **3** 🟡 *(computed)*.
- This fixes the *internal* anatomy of the bead oligo downstream of the handle, which is
  why the model draws the barcode region as one placeholder of unstated anatomy rather
  than claiming a resolved layout: the anatomy is supported, the bases of the flanking
  kit oligos are not. Promoting it would mean treating a third-party parser as an
  authority.
- **DNA Adapters** supplied: N701–N707, N718 (8-sample kit); N701–N707, N710–N712,
  N714–N716, N718–N724, N726–N729 (24-sample kit) 🟢 (Reference Guide, Kit Options).
  Illumina *Adapter Sequences* v14 prints the N7xx bases in both orientations 🟢:
  for example N701 is `TCGCCTTA` in the adapter and `TAAGGCGA` in the sample sheet.
  The named N7xx set is the Nextera XT/Nextera DNA i7 index family; this does not
  disclose the full SureCell index PCR primer or its chemical modifications. 🔴

### Primer roles the guide does fix 🟢

- The kit ships **one** custom primer, "Sequencing Primer (SP)", 50 µM, for **Read 1**
  only 🟢 — so Index 1 and Read 2 run on the instrument's standard primers 🟡 *(implied by
  the guide naming a custom primer for Read 1 alone)*. SP's bases are 🔴.
- The library is **single-indexed**: only an N7xx (i7) adapter is supplied, and no i5
  index is named anywhere in the kit documentation 🟢/🟡.
- The transposome is named "Nextera SureCell transposome" 🟢; which adapter it carries is
  not stated 🔴, so the model draws a `[transposome-derived end]` placeholder rather than
  asserting `nextera.ADAPTOR_S7`.

## 3. Step by step

Volumes and conditions 🟢 from the Reference Guide (agreeing with the checklist); the
molecular interpretation of each step 🟡 *(computed from the guide's reagent roles)*.

1. **Cells**: single-cell suspension in 1× PBS + 0.1 % BSA at **2500 cells/µL**
   (2250–2750), >95 % viable; 11,250 cells per sample, 4 samples per cartridge.
2. **Mixes** (on ice). *Cell Enzyme Mix* (red caps): Cell Suspend Buffer 60, DTT 8, RNA
   Stabilizer 6, RT Enzyme 13.2, Enhancer Enzyme 12 µL per cartridge. *Cell Suspension
   Mix*: 21.5 µL Cell Enzyme Mix + 4.5 µL cells per sample. *Barcode Suspension Mix*
   (blue caps): Barcode Buffer 60 + **3′ Barcode Mix** 60 µL per cartridge (the barcode
   beads). The RT enzyme is therefore already with the cells; lysis happens in the
   droplet. The identity of the "Enhancer Enzyme" is not given 🔴.
3. **Droplets** on the ddSEQ Single-Cell Isolator (Encapsulation Oil, Priming Solution);
   each sample's emulsion (35–40 µL) is split over two wells of a ddPCR plate.
4. **Lysis and RT in droplets**: 37 °C 30 min, 50 °C 60 min, 85 °C 5 min, 4 °C (lid
   105 °C, 50 µL). The bead oligo-dTV primes on the poly(A) tail and the first strand
   is `bead handle · bc1 · spacer1 · bc2 · spacer2 · bc3 · ACG · UMI · GAC · dT · cDNA`.
   🟡 What the 37 °C step does (lysis? releasing the oligos from the beads?) is not
   stated 🔴.
5. **Break emulsion**: 20 µL Droplet Disruptor, then 100 µL water, without mixing.
6. **First clean-up**: 90 µL Purification Beads (SPB) mixed into the aqueous layer only;
   two 80 % ethanol washes; elute in 35 µL RSB; the two wells of a sample are pooled
   (2 × 34 µL = 68 µL). The guide says this step removes **"unbound barcodes"** 🟢 —
   i.e. free barcode oligos, which implies they are in solution, not on beads, by this
   point 🟡.
7. **Second strand**: + 12 µL master mix (SSB 36 + SSE 18 µL per cartridge), **16 °C
   2 h**, lid off. The guide says it "removes the RNA template and synthesizes a
   replacement strand" 🟢 — which is the Gubler–Hoffman nick-translation operation, but
   Illumina names neither enzyme 🔴. Safe stopping point.
8. **cDNA clean-up**: 44 µL SPB (0.55× on 80 µL) 🟡 *(ratio computed)*, elute 11 µL,
   keep 10 µL. QC: Bioanalyzer HS, **>2 ng** cDNA needed (less over-tagments); typical
   cDNA 400–8000 bp.
9. **Tagment**: 30 µL Tagmentation Mix (TCB 88 + TCE 44 µL per cartridge) on the 10 µL
   cDNA, **55 °C 5 min** → 4 °C (remove within 6 min); stop with 10 µL Tagment Stop
   Buffer, 5 min room temperature. Transposome is the **"Nextera SureCell
   transposome"** 🟢; which adapter it carries is not stated 🔴, although a single-adapter
   (homodimer) load is what the single-indexed library implies 🟡. Each fragment that
   matters keeps the bead end on one side and gets ME + s7 on the other
   ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)).
10. **Library PCR**: + 30 µL Tagmentation PCR Mix (TPM), **10 µL TPP1**, **10 µL one
    DNA Adapter (N7xx) per sample** (100 µL). 95 °C 30 s; **15 cycles** of 95 °C 10 s,
    60 °C 45 s, 72 °C 60 s; 72 °C 5 min. TPP1 primes on the bead handle (copied into
    the second strand), N7xx on s7, so only **bead-end ↔ s7** fragments are amplified
    exponentially; internal s7–s7 fragments carry no P5 site 🟡. There is **no 72 °C
    gap-fill step**: the program starts at 95 °C 🟢. None is needed 🟡: of the two
    strands of a bead-end fragment, the one that Tn5 joined to s7-ME at its 5′ end
    ends at its 3′ end in the complement of the bead handle, so TPP1 primes on it in
    cycle 1 and copies s7-ME, creating the N7xx site; the gapped strand is simply lost.
11. **Library clean-up, double-sided**: 58 µL SPB (0.58×) → elute 51 µL, keep 50 µL; then
    add 30 µL SPB (0.6× on 50 µL) → elute 22 µL, keep 20 µL 🟢 volumes; ratios 🟡
    *(computed)*. Expected 2–10 nM, typical size ~300–1000 bp; normalise to 2 nM.

### Variants (demonstrated protocols) 🟢

- **PBMC** (#1000000044179): 3000 cells/µL, >80 % viable; all four chambers loaded,
  up to **two** samples per cartridge, ~500 cells per sample; cDNA yields as low as
  0.5 ng (mean 1.85 ng), and the tagmentation mix is scaled to the two samples (TCB 44 +
  TCE 22 µL per cartridge). Chemistry otherwise unchanged.
- **Nuclei** (#1000000044178): nuclei from cell lines, lysed on ice 10 min in 10 mM
  Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl₂, 0.1 % IGEPAL CA-630, SUPERase-In 0.2 U/µL,
  BSA 10 mg/mL; washed in PBS + 1 % DEPC and in storage buffer (same without IGEPAL);
  loaded at 2500 nuclei/µL. Chemistry otherwise unchanged.

## 4. Final library — 🟡 role-level layout, 🔴 most bases unresolved

Top strand, 5′→3′ (the strand whose 5′ end is P5):

```
5'- P5 (29) · GCCTGTCCGCGG (12) · bead handle (25) · [offset 0-3?] · bc1 (6) · spacer1 (15) · bc2 (6) · spacer2 (15) · bc3 (6) · ACG · UMI (8) · GAC · (dT)n V · cDNA (sense, 3' end of the transcript) · ME' (19) · s7' (15) · i7' (8) · P7' (24) -3'
```

Only the barcode/UMI interior of that line has evidence (the real reads, §2) and only
the i7 window has a vendor sequence; the 5′ region contributed by TPP1 and the 3′
transposome-derived region are 🔴, so `tools/surecell.py` draws them as
`[TPP1-selected barcode end]`, `[transposome-derived end]` and
`[N7xx adapter, including index]` and asserts no length. The lengths written above come
from the earlier drawn-page reconstruction and are kept only as an order-of-magnitude
sketch, not as a claim.

Strand bookkeeping 🟡 *(computed)*: the bead oligo *is* the first (antisense) cDNA
strand, and the top strand above is that strand with whatever TPP1 adds to its 5′ end.
So the cDNA written on the top strand is **antisense** to the mRNA, whatever a diagram
labelled "cDNA" suggests. Read 1 is primed on the bottom strand and copies the top
strand: barcodes, UMI, dT. Read 2 is primed on s7 + ME at the far end and copies the
bottom strand: it is **sense** to the mRNA, starting at the random Tn5 cut and reading
towards the polyadenylation site.

## 5. Read layout and sequencing

- **Read 1** (custom SP primer): reads barcode1, spacer1, barcode2, spacer2, barcode3,
  `ACG`, UMI, `GAC` and the start of the dT 🟡 *(from the real reads, §2)*. The ten umis
  example reads are all **68 nt** 🟡 — the only evidence for the read length, since
  Illumina publishes none 🔴.
- **Index 1** (standard Nextera i7 primer): **8 cycles**, the N7xx sample index. Nine of
  the ten example read headers carry `TAAGGCGA` (one has `TAAGGAGA`, one mismatch,
  presumably a read error) 🟡; `TAAGGCGA` matches the vendor N701 sample-sheet
  index exactly 🟢. No Index 2.
- **Read 2** (standard Nextera Read 2 primer): the cDNA, sense, starting at the Tn5 cut
  🟡; its length is not established by the checked guide 🔴.
- Illumina gives **no read lengths** in any document read 🔴; it requires Illumina
  Experiment Manager ≥ 1.13 / bcl2fastq ≥ 2.18 for "the appropriate UMI settings" and
  says SP is compatible "with this library and PhiX only" 🟢. Loading: NextSeq
  2.7–3 pM, HiSeq 2500 18–20 pM, HiSeq 3000/4000 350–400 pM, MiSeq 25–32 pM 🟢.

**A variable offset before barcode1.** 🟡 *(computed from
`umis_SureCell_K562_R1.fastq`)* Locating spacer1 in the ten example reads puts barcode1
at read position 1 in four reads, position 2 in two (preceding base `A`), and position 4
in three (preceding `GCA`); one read has no spacer1. The `umis` regex starts with a
greedy `(.*)` before barcode1 for the same reason. So Read 1 does **not** always start
with barcode1: the bead oligos carry a **0–3 nt (at least) phasing insert** before
barcode1, presumably to diversify the
first cycles of a low-complexity read. Its length set and sequences remain unverified from the checked guides
🔴. With a 3-nt offset only 3 nt of dT remain within 68 cycles. Other observations 🟡:
cycle 6 is `N` in every example read (a run-specific base-call dropout, not chemistry),
one read has `TAACTCTGAGCTGAA` for spacer2 (one mismatch), and one read lacks
spacer2–UMI entirely (a truncated bead oligo or synthesis error).

## 6. What each acceptable source supports

| Point | Illumina kit documents | real reads / `umis` |
|---|---|---|
| RT in droplets from bead oligo-dT | yes 🟢 | — |
| Second strand | "removes RNA, synthesizes replacement strand"; enzymes not named 🟢/🔴 | — |
| Nextera SureCell transposome | name only 🟢; adapter load 🔴 | — |
| PCR with TPP1 + N7xx, 15 cycles | yes 🟢 | — |
| Custom Read 1 primer only | yes 🟢; its sequence 🔴 | — |
| N7xx index bases | yes, *Adapter Sequences* v14 🟢 | N701 index seen in the read headers 🟡 |
| Bead layout bc1·sp1·bc2·sp2·bc3·ACG·UMI·GAC·dT | — | yes (regex + reads) 🟡 |
| Bead handle, TPP1, SP, spacer function, dT length, modifications | 🔴 | 🔴 |
| Read 1 starts at barcode1 | — | **no**: a 0–3 nt offset precedes it 🟡 |
| Read lengths | not given 🔴 | R1 = 68 nt in the ten example reads 🟡 |

## 7. Open questions

- 🔴 Every kit oligo sequence. Illumina printed none, and no acceptable source holds
  them; the model asserts none.
- 🔴 The phasing insert before barcode1: allowed lengths and bases.
- 🔴 Whether the barcode oligos are released from the beads (the "unbound barcodes" in
  the first clean-up suggests so), how, and what the 37 °C 30 min step and the
  "Enhancer Enzyme" do.
- 🔴 The transposome's adaptor load (s7 only, or s5 + s7 with the s5 ends simply lost).
- 🔴 The `(dT)` length, and the purpose of the `ACG` / `GAC` trinucleotides flanking the
  UMI (these two are visible in the real reads 🟡; everything 5′ of barcode1, including
  the bead handle and whatever TPP1 adds, is not).
- 🟡 That the library PCR relies on the un-gapped strand only (no gap fill) — inferred.
- 🔴 Full SureCell N7xx index PCR primer structure and modifications remain undisclosed;
  the index bases themselves are in Illumina *Adapter Sequences* v14.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` found the drawn page and the scg mirrors (the catalogue row has no
DOI: a vendor kit). The four Illumina PDFs were found on the kit's support "documentation"
page and the three `umis` files on GitHub, all fetched by hand into the source directory;
`tools/doctext.py` made the text twins. `tools/scrape_primers.py` found no proprietary
SureCell bead, TPP1 or SP sequence in the kit PDF; the separate Illumina Adapter
Sequences guide does print N7xx index bases. The scraper does not scan `.json` /
`.fastq` files, so the umis and real-read evidence was parsed by hand.
