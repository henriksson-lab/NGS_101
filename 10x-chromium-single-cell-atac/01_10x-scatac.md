# 10x Chromium Single Cell ATAC — bulk Tn5 tagging, droplet barcoding by linear PCR

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`). A claim found **only**
> on the upstream scg_lib_structs page is 🟡 (secondary source), even when quoted
> verbatim from it.

**10x Genomics Chromium Single Cell ATAC** (kit v1; Chip E). A commercial kit: there is
**no defining publication** (catalogue row: "vendor/kit protocol: no defining
publication"). Primary source:

- 10x Genomics. *Chromium Single Cell ATAC Reagent Kits User Guide*, CG000168 **Rev A**
  (© 2018). Kits PN-1000110 / PN-1000111 (Library & Gel Bead), PN-1000082 / PN-1000086
  (Chip E), PN-1000084 (i7 Multiplex Kit N, Set A). No DOI.

Papers upstream cites for the concepts (not for the kit itself):

- Amini S, Pushkarev D, Christiansen L, *et al.*, Shendure J, Gunderson KL, Steemers FJ.
  "Haplotype-resolved whole-genome sequencing by contiguity-preserving transposition and
  combinatorial indexing." *Nat Genet* 2014;46:1343–1349.
  doi:[10.1038/ng.3119](https://doi.org/10.1038/ng.3119) (PMC4409979) — **Tn5 stays bound
  after transposition** (CPT-seq), the reason nuclei can be tagged in bulk and stay intact.
- Plessy C, Bertin N, Takahashi H, *et al.* "Linking promoters to functional transcripts
  in small samples with nanoCAGE and CAGEscan." *Nat Methods* 2010;7:528–534.
  doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470) (PMC2906222) — origin of the
  term **semi-suppressive PCR**.
- Secondary: Teichmann lab, scg_lib_structs, "10x Chromium Single Cell ATAC"
  (<https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html>).

Sources read (in `$CHEM_DATA/sources/10x-chromium-single-cell-atac/`, never committed):

| File | What | Used for |
|---|---|---|
| `upstream_10xChromium_scATAC.html.txt` | scg_lib_structs page (fetched by `get_sources.py`) | complete second account: oligos, steps, final library, read layout |
| `CG000168_ChromiumSingleCellATAC_ReagentKits_UserGuide_RevA.pdf(.txt)` | 10x User Guide Rev A — fetched **by hand** from the copy upstream links (`teichlab.github.io/scg_lib_structs/data/10X-Genomics/`) | **every oligo** (Appendix "Oligonucleotide Sequences", txt lines 3343–3441), all reaction conditions, sequencing parameters |
| `ng.3119_PMC4409979.html.txt` | Amini 2014 full text (`get_sources.py`) | Tn5-stays-bound concept; CPT-seq release by SDS |
| `ng.3119_41588_2014_BFng3119_MOESM1_ESM.pdf.txt`, `…MOESM2_ESM.xlsx.txt` | Amini 2014 supplement + oligo table | checked: CPT-seq's indexed transposons, not used by 10x |
| `nmeth.1470_PMC2906222.html.txt` | nanoCAGE full text — fetched **by hand** from PMC (`get_sources.py` did not fetch this "cited" row) | definition of semi-suppressive PCR |

Not obtained:

| What | Why it matters | Where |
|---|---|---|
| User guides for **v1.1** and **v2** | upstream says all three versions give the same library; only Rev A (v1) was read | 10x support site (manual) 🔴 |
| Composition of ATAC Buffer, ATAC Enzyme, Barcoding Reagent / Enzyme, Reducing Agent B | whether the transposome is pre-loaded Tn5 with both adaptors; what releases Tn5 in the GEM | proprietary 🔴 |
| The 4 sequences in each Chromium i7 Sample Index well | sample-sheet index sequences | 10x website index CSV (manual) 🔴 |
| Amini 2014 PMC author supplement `NIHMS681868-supplement-Figures___Tables.doc` | not needed (Springer copy fetched) | <https://pmc.ncbi.nlm.nih.gov/articles/instance/4409979/bin/NIHMS681868-supplement-Figures___Tables.doc> |

---

## 1. What it is

Droplet single-nucleus ATAC-seq. **Nuclei are tagmented in bulk first**, then partitioned
one per GEM (Gel Bead-in-emulsion) on Chip E, where a gel bead releases a barcoded
primer that adds a **16-nt 10x Barcode** to every Tn5 fragment by **linear (single
primer) amplification**. The droplets are broken, and a bulk PCR adds P7 and an 8-nt
**i7 sample index**. 🟢 (CG000168, "Stepwise Objectives")

There is **no UMI**: the 16-nt cell barcode is the only random/designed sequence added in
the droplet. Duplicates are resolved by fragment coordinates. 🟡 (no UMI appears in any
oligo; the dedup rule is Cell Ranger ATAC behaviour, not stated in the guide)

| | New thing here | Builds on |
|---|---|---|
| 1 | **Tag first, fragment later**: Tn5 inserts adaptors into chromatin in intact nuclei and stays bound, so a nucleus can be moved into a droplet carrying its own tagmented DNA | CPT-seq (Amini 2014); [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |
| 2 | Cell barcode added by a **bead primer that ends in s5 only** (14 nt), priming on the Tn5 adaptor rather than on any capture handle | the same s5 entry point as Nextera N5xx primers (`nextera.n5xx_primer`) |
| 3 | Barcode sits in the **i5 index position** (16 nt), so the library is read as a standard dual-index Nextera library with standard primers | sci-ATAC-seq / dsciATAC also use index reads for cell identity |

The library is not an RNA library: no reverse transcription, template switching, ligation
or circularization is involved.

## 2. Oligos

🟢 Verbatim from CG000168 Rev A, Appendix "Oligonucleotide Sequences" (txt lines
3343–3441); the hyphens are 10x's segment separators, not modifications. The guide gives
**no modifications** for any oligo (no phosphorylation, no linker chemistry for the bead
attachment). Upstream writes the bead attachment as `|--5'-`. 🟡

```
Read 1N primer sequence (Transposition Mix)     5'-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-3'
Read 2N primer sequence (Transposition Mix)     5'-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-3'
Gel Bead Oligo Primer PN-2000132                5'-AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTC-3'
SI-PCR Primer B PN-2000128 (Forward Primer)     5'-AATGATACGGCGACCACCGAGA-3'
i7 Sample Index Plate N, Set A PN-3000262       5'-CAAGCAGAAGACGGCATACGAGAT-NNNNNNNN-GTCTCGTGGGCTCGG-3'
  (Reverse Primer)
```

The guide labels the bead oligo segments "P5 · 10x Barcode · Partial Read 1N", SI-PCR
Primer B "Partial P5", and the i7 primer "P7 · Sample Index N · Partial Read 2N". 🟢 Each
i7 well is "a unique mix of 4 oligos", base-balanced; the sequences are not in the guide. 🔴

Upstream additionally lists the two index-read sequencing primers (not in CG000168) 🟡:

```
Sample index sequencing primer (index1)   CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Cell barcode sequencing primer (index2)   CTGTCTCTTATACACATCTGACGCTGCCGACGA
```

### How they interlock — 🟡 (computed against `lib/`)

- **Read 1N = `nextera.ADAPTOR_S5`** (s5 + ME, 33 nt) and **Read 2N = `nextera.ADAPTOR_S7`**
  (s7 + ME, 34 nt): the transposome carries the ordinary Nextera adaptors; the guide calls
  them "Read 1N/2N primer sequence" because the same sequences are the sequencing primers.
  They are therefore also `nextera.READ1_PRIMER` / `READ2_PRIMER`.
- **Gel Bead oligo = `illumina.P5` (29 nt) + 16 N + `nextera.S5` (14 nt) = 59 nt** —
  i.e. exactly `nextera.n5xx_primer(<16-nt barcode>)`: a Nextera i5 index primer with a
  16-nt instead of 8-nt index. It ends **before** the ME; the ME bases of the product come
  from the template.
- **SI-PCR Primer B = `illumina.P5[:22]`**, 22 nt; it lacks the last 7 nt of P5
  (`TCTACAC`), so it primes on the outermost P5 bases only and cannot add anything new.
- **i7 primer = `nextera.n7xx_primer(<8-nt index>)`** = `illumina.P7` + 8 N + `nextera.S7`,
  47 nt: a standard Nextera N7xx primer. Anneals over s7 (15 nt) only.
- Upstream's **index1 primer = `nextera.INDEX1_PRIMER` = revcomp(Read 2N)** and **index2
  primer = `nextera.INDEX2_PRIMER` = revcomp(Read 1N)**: standard Nextera index primers.
- The transposed-product strand ends given by 10x — top 3' end
  `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC` = revcomp(Read 2N), bottom strand written 3'→5'
  `AGCAGCCGTCGCAGTCTACACATATTCTCTGTC` = revcomp(Read 1N) — confirm the ordinary s5/s7
  tagmentation product.

## 3. Step by step

Conditions 🟢 from CG000168 Rev A unless marked.

1. **Nuclei** isolated by a separate Demonstrated Protocol (not read 🔴), resuspended in
   Diluted Nuclei Buffer (PN-2000153 diluted 1:20); the guide warns the buffer
   composition, "including Magnesium concentration", is optimised for the transposition
   and barcoding steps. Nuclei stock volume = Targeted Nuclei Recovery × 1.53.
2. **Transposition (bulk)**: 7 µl ATAC Buffer + 3 µl ATAC Enzyme + 5 µl nuclei/buffer =
   15 µl; **37 °C 60 min**, hold 4 °C (lid 50 °C). Tn5 inserts s5-ME and s7-ME adaptors
   into accessible chromatin with the usual 9-bp gap — see
   [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). Three end combinations arise
   (s5/s5, s7/s7, s5/s7); only s5/s7 becomes a sequenceable library. 🟡 (upstream draws
   all three; `nextera.TAGMENTATION_OUTCOMES`)
   - The DNA is **not released**: Tn5 stays bound after transposition until denatured
     (Amini 2014: dissociates only with SDS or protease 🟢 in Amini), so nuclei stay intact
     and keep their own fragments. 🟡 for the 10x case (the guide only says the transposase
     "fragments" DNA in open regions).
3. **GEM generation**: 65 µl Master Mix (61.5 µl Barcoding Reagent, 1.5 µl Reducing
   Agent B, 2 µl Barcoding Enzyme) added to the 15 µl transposed nuclei; 75 µl loaded,
   plus 40 µl Gel Beads and Partitioning Oil, on Chip E. ~90–99 % of GEMs contain no
   nucleus. The gel bead dissolves and releases the barcoded oligo. 🟢 Reducing Agent B
   in the mix is presumably what dissolves the (disulfide-crosslinked) bead. 🟡
4. **GEM incubation** (lid 105 °C, 125 µl setting): **72 °C 5 min**; 98 °C 30 s;
   **12 × (98 °C 10 s, 59 °C 30 s, 72 °C 1 min)**; 15 °C hold. 🟢
   - 72 °C 5 min = **gap fill-in** by the Barcoding Enzyme before any denaturation
     (`nextera.GAP_FILL_TEMP_C`), completing the s5/s7 adaptors on both strands. 🟡
   - Tn5 release in the droplet: by heat (72 → 98 °C) and/or a denaturant in the
     Barcoding Reagent. 🔴 (not stated; Amini 2014 used SDS at 55 °C, which 10x does not
     list)
   - 12 cycles with **one primer** (the bead oligo) — "Linear Amplification" in the
     guide 🟢. The bead primer's 3' `TCGTCGGCAGCGTC` anneals to the s5 complement at the
     3' end of each strand that carries an s5 adaptor, and extends through ME, insert and
     the far adaptor 🟡 (mechanism inferred). Product (10x's own drawing):
     `P5 · <cbc16> · Read 1N · insert · revcomp(Read 2N)`. 🟢
5. **Post-GEM cleanup**: 125 µl Recovery Agent, remove; Dynabeads MyOne SILANE cleanup
   (200 µl cleanup mix, 10 min RT, two 80 % ethanol washes), elute 40.5 µl Elution
   Solution I (Buffer EB + 10 % Tween 20 + Reducing Agent B); then **SPRIselect 1.2×** (48 µl to 40 µl),
   elute in 40.5 µl EB, keep 40 µl. The SPRI step removes unused bead primers. 🟢
6. **Sample Index PCR** (100 µl): 40 µl sample + 50 µl Amp Mix + 7.5 µl SI-PCR Primer B +
   2.5 µl one i7 Sample Index well. 98 °C 45 s; **N × (98 °C 20 s, 67 °C 30 s, 72 °C 20 s)**;
   72 °C 1 min. N = 10–12 (cell lines) or 11–13 (primary cells) by targeted recovery
   (500–2,000: 12/13; 2,001–6,000: 11/12; 6,001–10,000: 10/11). 🟢 The i7 primer
   anneals to s7 complement at the 3' end of the barcoded strand and adds
   `<i7 index>` + P7; SI-PCR Primer B then amplifies from P5. 🟡
7. **Double-sided size selection**: 40 µl SPRIselect to the 100 µl PCR (0.4×, keep
   130 µl supernatant), then +74 µl SPRIselect to the supernatant (bind), two 80 % ethanol
   washes, elute 20.5 µl EB, keep 20 µl. 🟢 (The guide gives volumes; the cumulative ratio
   is not stated.)
8. QC and **qPCR quantification** (KAPA Library Quantification Kit; 95 °C / 67 °C, 30 cycles); Bioanalyzer High Sensitivity or TapeStation. 🟢

### Why only s5/s7 fragments become library — 🟡

- **s7/s7**: no s5 complement anywhere, so the bead primer never extends on it; no
  barcode, no P5 → not sequenced.
- **s5/s5**: the bead primer can barcode **both** strands, giving
  `P5 · <cbc> · Read 1N · insert · revcomp(Read 1N) · revcomp(<cbc>) · revcomp(P5)` with
  inverted repeats at the two ends (same barcode, same GEM). Upstream calls these "not
  amplifiable due to semi-suppressive PCR" (the ends fold into a panhandle); in any case
  they carry no P7 and cannot form clusters.
- **s5/s7**: one barcoded strand per template per cycle, ending in revcomp(Read 2N) — the
  substrate for the i7 primer.

nanoCAGE (Plessy 2010) defines semi-suppressive PCR for linkers that are "similar (but not
identical)" at the two ends, with identical-end templates forming homo-duplexes. 🟢 in
Plessy 2010; its application to scATAC is upstream's interpretation. 🟡

## 4. Final library — 🟢 (CG000168 Appendix, "Sample Index PCR Product")

```
5'- P5 · <cbc16> · s5 · ME · <genomic insert> · ME' · s7' · <i7 index 8> · P7' -3'
```

As sequence (10x's own top strand, hyphens added by 10x):

```
AATGATACGGCGACCACCGAGATCTACAC-NNNNNNNNNNNNNNNN-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG---insert---CTGTCTCTTATACACATCTCCGAGCCCACGAGAC-NNNNNNNN-ATCTCGTATGCCGTCTTCTGCTTG
```

Segment lengths 🟡 (computed): `illumina.P5` 29 + barcode 16 + Read 1N 33 + insert +
revcomp(Read 2N) 34 + i7 8 + `illumina.P7_RC` 24 = **144 nt of adaptor** plus the insert.
The tail `ATCTCGTATGCCGTCTTCTGCTTG` = `illumina.P7_RC`. 🟡 (computed)

Upstream's final structure is identical, segment for segment and base for base. 🟡
(agreement checked against CG000168)

## 5. Sequencing

🟢 (CG000168 "Sequencing"): paired-end, dual indexing; **Read 1N 50 cycles, i7 index 8
cycles, i5 index 16 cycles, Read 2N 50 cycles**; 25,000 read pairs per nucleus. Read 1N
and Read 2N contain insert only (from opposite ends of the fragment); the 8-bp sample
index is the i7 read; the **16-bp 10x Barcode is the i5 read**. Verified on MiSeq,
NextSeq 500/550 HO, HiSeq 2500 RR, HiSeq 3000/4000, NovaSeq; ~1 % PhiX; "DO NOT pool
Single Cell ATAC libraries with other 10x Genomics libraries".

Upstream gives the same cycle numbers (50 / 8 / 16 / 50) 🟡 (agreement). Read primers are
standard: Read 1 = Read 1N = `nextera.READ1_PRIMER`, Read 2 = Read 2N =
`nextera.READ2_PRIMER`, index 1 = `nextera.INDEX1_PRIMER`, index 2 =
`nextera.INDEX2_PRIMER` 🟡 (computed; the guide names no custom primers). Because the
barcode is read as i5, its orientation in the FASTQ depends on the instrument's i5
workflow (forward vs reverse complement); Cell Ranger ATAC handles this. 🟡 (not in the
guide)

Because ME is part of the read primer, Read 1N and Read 2N start at the first genomic
base after the ME — the 9-bp target-site duplication is read at both ends. 🟡

## 6. Upstream vs 10x guide — agreements and disagreements

| Point | Upstream | CG000168 Rev A | Verdict |
|---|---|---|---|
| Bead oligo, SI-PCR Primer B, i7 primer, Read 1N/2N | as §2 | identical | agree 🟢 |
| Final library | as §4 | identical | agree 🟢 |
| Read lengths | 50/8/16/50 | 50/8/16/50 | agree 🟢 |
| Number of barcodes | 737,280 (`737K-cratac-v1.txt.gz`) | "~750,000 10x Barcodes" | consistent; exact count 🟡 (upstream only) |
| v1, v1.1, v2 have the same library | stated | Rev A covers v1 only | unchecked 🔴 |
| Tn5 release "using heat or denaturing agents such as SDS", attributed to Amini 2014 | stated | not stated | **partial disagreement**: Amini 2014 names SDS and protease, not heat 🟢 (Amini); heat release in the GEM is plausible but unsourced 🔴 |
| Index-read primer sequences | given | not given | 🟡 (standard Nextera, computed) |
| Gap fill-in as first step of the linear PCR at 72 °C | stated | 72 °C 5 min first step, purpose not stated | consistent 🟡 |

## 7. Open questions

- 🔴 What releases Tn5 from the DNA in the GEM (heat alone, or a detergent in the
  Barcoding Reagent)? What polymerase is the Barcoding Enzyme?
- 🔴 Whether v1.1 (Next GEM) and v2 changed any oligo, condition or barcode list (upstream
  says the library is the same).
- 🔴 Exact sequences of the i7 Sample Index Set A wells (4 oligos each).
- 🟡 Barcode-to-whitelist orientation across instruments — software, not chemistry.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py "10x-chromium-single-cell-atac"` fetched the upstream page and
Amini 2014 (PMC text + Springer supplements; the PMC author `.doc` went to `(manual)`).
It did **not** fetch the nanoCAGE paper (role "cited") nor the 10x user guide that the
upstream page links; both were fetched by hand (`curl`) and converted with
`tools/doctext.py`. `tools/scrape_primers.py --max-hits 80` found every 10x oligo in the
guide's Appendix, rejoining `-NNNN-` runs. Its upstream hits include one chimera: on
upstream line 67 it joined the two strands of a double-strand drawing into a 38-nt
"oligo" (`CTGTCTCTTATACACATCT` + `TCTACACATATTCTCTGTC`). `--find` reports 0 locations for
any query containing N (e.g. the bead oligo with its 16 N), so the N-containing oligos
were checked by their fixed segments with `--find` and by reading the Appendix lines
directly.
