# SureCell WTA 3′ (Illumina / Bio-Rad ddSEQ): droplet 3′ RNA-seq with split-pool bead barcodes

> **Evidence marking.** 🟢 verbatim from a primary source (here: Illumina's own kit
> documents) · 🟡 derived, inferred, or taken from a secondary source · 🔴 not published.
> **Illumina never published the oligo sequences.** Every sequence below comes from the
> upstream scg_lib_structs page (a secondary source), cross-checked where possible
> against the read structure that the `umis` package parses out of real SureCell
> reads; all sequences are therefore 🟡, not 🟢. Relationships marked 🟡 *(computed)*
> were worked out with `lib/` while writing this note; they are not yet asserted in a
> self-test, because this protocol has no `tools/` module yet.

**Illumina Bio-Rad SureCell WTA 3′ Library Prep Kit for the ddSEQ System** — a vendor
kit (Illumina cat. 20014279 / 20014280), run on the Bio-Rad **ddSEQ Single-Cell
Isolator**. There is **no defining publication**; the reference is the kit
documentation: *SureCell WTA 3′ Library Prep Reference Guide*, Illumina document
#1000000021452 v01 (June 2017; v00 February 2017). Analysis was the BaseSpace
"SureCell RNA Single-Cell" app. The kit is **discontinued**; Bio-Rad's successor is the
*ddSEQ Single-Cell 3′ RNA-Seq Kit*, which is a different chemistry and is not covered
here (upstream page, header note).

Sources read (`_data/sources/surecell-3-wta-for-ddseq/`, never committed). Only the
first was fetched by `tools/get_sources.py`; the others were fetched by hand from the
Illumina support page and the `vals/umis` GitHub repository and added to `MANIFEST.tsv`:

| File | What | Used for |
|---|---|---|
| `upstream_SureCell.html.txt` | scg_lib_structs method page (Teichlab) | **every oligo sequence**, library and read layout (secondary) |
| `surecell-wta3-library-prep-reference-guide-1000000021452-01.pdf.txt` | Illumina Reference Guide v01 | **all steps and conditions**, kit contents, index set, sequencing primer |
| `surecell-wta3-library-prep-checklist-1000000021454-00.pdf.txt` | Illumina checklist | same steps, condensed (agrees with the guide) |
| `surecell-wta-3-pbmc-demonstrated-protocol_1000000044179-00.pdf.txt` | PBMC demonstrated protocol | variant: PBMC input, half tagmentation mix |
| `surecell-wta3-nuclei-demonstrated-protocol-1000000044178-00.pdf.txt` | nuclei demonstrated protocol | variant: nuclei instead of cells |
| `umis_SureCell_transform.json` | `umis` read-1 parsing regex | independent check of the bead-oligo layout |
| `umis_SureCell_barcodes.txt` | 96 six-mers | the bead barcode whitelist |
| `umis_SureCell_K562_R1.fastq` | 10 real SureCell read 1s (K562) | read-1 length, the variable 5′ offset (§5) |

Not obtained 🔴: Illumina's *Adapter Sequences* document (#1000000002694; the old PDF
URLs 404, current version only as an HTML help site,
<https://support-docs.illumina.com/SHARE/adapter-sequences.htm>) — would give the N7xx
index sequences; the *Custom Primers Guides* for each instrument (read-primer
concentrations); the Bio-Rad *ddSEQ Single-Cell Isolator Instruction Manual*
(#10000069430); Illumina's ddSEQ technical / application notes. No vendor document
seen gives any oligo sequence.

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

## 2. Oligos

🟡 All sequences verbatim from the upstream scg_lib_structs page (secondary source); the
Illumina documents name the reagents (TPP1, DNA Adapters N7xx, Sequencing Primer SP) but
never print a sequence. The upstream page gives no chemical modifications; whether the
bead oligo is attached by its 5′ end through a linker, and whether it is released from
the bead, is not stated 🔴 (upstream marks it `|--5'`, i.e. bead-attached at 5′).

```
Beads-read1-oligo-dTV       |--5'- AAGCAGTGGTATCAACGCAGAGTAC[6-bp barcode1]TAGCCATCGCATTGC[6-bp barcode2]TACCTCTGAGCTGAA[6-bp barcode3]ACG[8-bp UMI]GAC(dT)V -3'
Spacer 1                    TAGCCATCGCATTGC
Spacer 2                    TACCTCTGAGCTGAA
Tagment PCR Adapter (TPP1)  AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
DNA Adapters (N7xx)         CAAGCAGAAGACGGCATACGAGAT[8-bp sample index]GTCTCGTGGGCTCGG
Read 1 sequencing primer    GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
Index sequencing primer     CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Read 2 sequencing primer    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
Transposon end (ME)         AGATGTGTATAAGAGACAG   (Nextera SureCell transposome; s7-ME assumed)
```

- Length of the `(dT)` stretch: not given 🔴.
- **DNA Adapters** supplied: N701–N707, N718 (8-sample kit); N701–N707, N710–N712,
  N714–N716, N718–N724, N726–N729 (24-sample kit) 🟢 (Reference Guide, Kit Options).
  The 8-nt index sequences are not in any source read 🔴 (they are the standard Nextera
  i7 set by name; Illumina Adapter Sequences document).
- The kit ships **one** custom primer, "Sequencing Primer (SP)", 50 µM, for **Read 1**
  only 🟢; Index 1 and Read 2 use the instrument's standard primers 🟡 (upstream, and
  implied: the guide says the custom primer is for Read 1).
- **Bead barcodes**: `umis_SureCell_barcodes.txt` lists 96 distinct 6-mers 🟡; in the
  nine parseable example reads, every barcode1, barcode2 and barcode3 matches exactly
  one of them (N as wildcard), so one 96-list serves all three positions 🟡 *(computed;
  9 reads only)*. Minimum pairwise Hamming distance within the list is **3** — one
  substitution can be corrected 🟡 *(computed)*.

### How they interlock — 🟡 (computed against `lib/`)

- **The bead handle `AAGCAGTGGTATCAACGCAGAGTAC` (25 nt) = `rt.SMART_HANDLE` + `AC`.** It
  is the SMARTer/ISPCR handle, as on Drop-seq / SMART-seq oligo-dT, but here it is a
  PCR and sequencing-primer site only — there is no TSO carrying the same handle.
- **TPP1 (66 nt) = `illumina.P5` + `GCCTGTCCGCGG` + bead handle.** It anneals to the
  complement of the bead handle and adds P5.
- **Read 1 primer (37 nt) = TPP1 minus P5**, i.e. `GCCTGTCCGCGG` + bead handle. It ends
  exactly at the last base of the handle, so **cycle 1 of Read 1 is the first base
  after the handle** — barcode1 (or the variable offset, §5). The 12-nt `GCCTGTCCGCGG`
  extension is unexplained 🔴; it lengthens the primer and probably raises its Tm. It
  contains no `nextera.S5` — the library has **no s5 at all**.
- **DNA Adapter N7xx = `illumina.P7` + i7 + `nextera.S7`** (the standard Nextera
  i7 primer truncated after s7; 24 + 8 + 15 = 47 nt). It primes on the s7 that the
  transposome left at the far end of each fragment.
- **Read 2 primer = `nextera.READ2_PRIMER` = `nextera.ADAPTOR_S7`** (s7 + ME, 34 nt).
- **Index primer = `nextera.INDEX1_PRIMER` = revcomp(`ADAPTOR_S7`)**: the standard Nextera
  i7 index primer. The library carries **no i5** (P5 is fused directly to the custom
  read-1 region), so it is a **single-indexed** library.
- The bead oligo from the handle through `GAC` (without dT) is **87 nt**; from barcode1
  to the end of `GAC` is **62 nt** (6 + 15 + 6 + 15 + 6 + 3 + 8 + 3).

## 3. Step by step

Volumes and conditions 🟢 from the Reference Guide (agreeing with the checklist); the
molecular interpretation of each step 🟡 (upstream page and the oligo layout).

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
   replacement strand" 🟢; upstream calls it RNase H + DNA Pol I (Gubler–Hoffman nick
   translation) 🟡 — consistent, but the enzymes are not named by Illumina. Safe
   stopping point.
8. **cDNA clean-up**: 44 µL SPB (0.55× on 80 µL) 🟡 *(ratio computed)*, elute 11 µL,
   keep 10 µL. QC: Bioanalyzer HS, **>2 ng** cDNA needed (less over-tagments); typical
   cDNA 400–8000 bp.
9. **Tagment**: 30 µL Tagmentation Mix (TCB 88 + TCE 44 µL per cartridge) on the 10 µL
   cDNA, **55 °C 5 min** → 4 °C (remove within 6 min); stop with 10 µL Tagment Stop
   Buffer, 5 min room temperature. Transposome is the **"Nextera SureCell
   transposome"** 🟢; upstream infers it is loaded only with the **s7-ME** adaptor
   ("highly likely a Tn5 homodimer") 🟡 — Illumina does not say. Each fragment that
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

## 4. Final library — 🟡 (assembled from the oligos; agrees with the upstream drawing)

Top strand, 5′→3′ (the strand whose 5′ end is P5):

```
5'- P5 (29) · GCCTGTCCGCGG (12) · bead handle (25) · [offset 0-3?] · bc1 (6) · spacer1 (15) · bc2 (6) · spacer2 (15) · bc3 (6) · ACG · UMI (8) · GAC · (dT)n V · cDNA (sense, 3' end of the transcript) · ME' (19) · s7' (15) · i7' (8) · P7' (24) -3'
```

Computed checks 🟡: upstream's final-library 5′ part equals TPP1 followed by the bead
oligo after its handle, and its 3′ tail equals `nextera.ME_RC` + `nextera.S7_RC` + i7 +
`illumina.P7_RC`, i.e. the reverse complement of N7xx + ME; upstream's bottom-strand line
is the exact complement of its top-strand line. Fixed length outside cDNA and dT:
29 + 12 + 25 + 62 + 19 + 15 + 8 + 24 = **194 nt** (without the offset) 🟡 *(computed)*.

Strand bookkeeping 🟡: the bead oligo *is* the first (antisense) cDNA strand, and the
top strand above is that strand with P5 + the 12-nt extension added to its 5′ end by
TPP1. So the cDNA written on the top strand is **antisense** to the mRNA (the upstream
"cDNA" label notwithstanding). Read 1 is primed on the bottom strand and copies the top
strand: barcodes, UMI, dT. Read 2 is primed on s7 + ME at the far end and copies the
bottom strand: it is **sense** to the mRNA, starting at the random Tn5 cut and reading
towards the polyadenylation site.

## 5. Read layout and sequencing

- **Read 1** (custom SP primer, ends at the handle): **68 cycles** (upstream) — reads
  barcode1, spacer1, barcode2, spacer2, barcode3, `ACG`, UMI, `GAC` and the start of the
  dT. 🟡 The umis example reads are all **68 nt** 🟡, agreeing.
- **Index 1** (standard Nextera i7 primer): **8 cycles**, the N7xx sample index. Nine of
  the ten example read headers carry `TAAGGCGA` (one has `TAAGGAGA`, one mismatch,
  presumably a read error) 🟡 (N701 by name; not checked against an
  index table). No Index 2.
- **Read 2** (standard Nextera Read 2 primer, s7 + ME): **75 cycles** (upstream) — the
  cDNA, sense, starting at the Tn5 cut.
- Illumina gives **no read lengths** in any document read 🔴; it requires Illumina
  Experiment Manager ≥ 1.13 / bcl2fastq ≥ 2.18 for "the appropriate UMI settings" and
  says SP is compatible "with this library and PhiX only" 🟢. Loading: NextSeq
  2.7–3 pM, HiSeq 2500 18–20 pM, HiSeq 3000/4000 350–400 pM, MiSeq 25–32 pM 🟢.

**Disagreement with upstream — a variable offset before barcode1.** 🟡 *(computed from
`umis_SureCell_K562_R1.fastq`)* Locating spacer1 in the ten example reads puts barcode1
at read position 1 in four reads, position 2 in two (preceding base `A`), and position 4
in three (preceding `GCA`); one read has no spacer1. The `umis` regex starts with a
greedy `(.*)` before barcode1 for the same reason. So Read 1 does **not** always start
with barcode1, as the upstream drawing shows: the bead oligos carry a **0–3 nt (at
least) phasing insert** between the handle and barcode1, presumably to diversify the
first cycles of a low-complexity read. Its length set and sequences are not published
🔴. With a 3-nt offset only 3 nt of dT remain within 68 cycles. Other observations 🟡:
cycle 6 is `N` in every example read (a run-specific base-call dropout, not chemistry),
one read has `TAACTCTGAGCTGAA` for spacer2 (one mismatch), and one read lacks
spacer2–UMI entirely (a truncated bead oligo or synthesis error).

## 6. Agreement between sources

| Point | Upstream page | Illumina documents | umis data |
|---|---|---|---|
| RT in droplets from bead oligo-dT | yes | yes 🟢 | — |
| Second strand by RNase H / Pol I | yes | "removes RNA, synthesizes replacement strand", enzymes not named | — |
| Nextera SureCell transposome | s7-ME only (inferred) | name only | — |
| PCR with TPP1 + N7xx | yes | yes, 15 cycles 🟢 | — |
| Custom Read 1 primer only | yes | yes 🟢 | — |
| Bead layout bc1·sp1·bc2·sp2·bc3·ACG·UMI·GAC·dT | yes | — | yes (regex, reads) |
| Read 1 starts at barcode1 | yes | — | **no**: 0–3 nt offset |
| Read lengths 68 / 8 / 75 | yes | not given | R1 = 68 |

## 7. Open questions

- 🔴 Every oligo sequence rests on a secondary source; Illumina printed none.
- 🔴 The phasing insert before barcode1: allowed lengths and bases.
- 🔴 Whether the barcode oligos are released from the beads (the "unbound barcodes" in
  the first clean-up suggests so), how, and what the 37 °C 30 min step and the
  "Enhancer Enzyme" do.
- 🔴 The transposome's adaptor load (s7 only, or s5 + s7 with the s5 ends simply lost).
- 🔴 The `(dT)` length, and the purpose of `ACG` / `GAC` flanking the UMI and of the
  12-nt `GCCTGTCCGCGG` in TPP1 and the Read-1 primer.
- 🟡 That the library PCR relies on the un-gapped strand only (no gap fill) — inferred.
- 🟡 N7xx index sequences assumed to be the standard Nextera i7 set — unchecked.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` found only the upstream page (the catalogue row has no DOI: a
vendor kit). The four Illumina PDFs were found on the kit's support "documentation"
page and the three `umis` files on GitHub, all fetched by hand into the source
directory; `tools/doctext.py` made the text twins. `tools/scrape_primers.py` found
every sequence on the upstream page (and recognised `rt.SMART_HANDLE`, `illumina.P5`,
`nextera.ADAPTOR_S7`, `nextera.INDEX1_PRIMER`), and correctly found none in the Illumina
PDFs. It does not scan `.json` / `.fastq` files, so the umis evidence was parsed by hand.
