# scDNase-seq (Pico-Seq) and scMNase-seq — nuclease digestion of one cell, then a plain ligation library

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note. They are not asserted in a self-test because the primary papers do
> not print the oligo sequences.
>
> **Important for these protocols:** neither defining paper prints a library oligo.
> Pico-Seq's cited library method identifies the exact historical Illumina kit, whose
> oligos are vendor-published and therefore 🟢 below. scMNase-seq only points to an
> older universal-adapter method; its reconstructed endpoint remains 🟡.

**scDNase-seq** (the paper calls it **Pico-Seq**) — Jin W, Tang Q, Wan M, Cui K, Zhang Y,
Ren G, Ni B, Sklar J, Przytycka TM, Childs R, Levens D, Zhao K. "Genome-wide detection of
DNase I hypersensitive sites in single cells and FFPE tissue samples." *Nature* 528,
142–146 (2015). doi:[10.1038/nature15740](https://doi.org/10.1038/nature15740) · PMC4697938.
Data: GEO GSE61844 (Pico-Seq and RNA-seq).

**scMNase-seq** — Lai B, Gao W, Cui K, Xie W, Tang Q, Jin W, Hu G, Ni B, Zhao K.
"Principles of nucleosome organization revealed by single-cell micrococcal nuclease
sequencing." *Nature* 562, 281–285 (2018).
doi:[10.1038/s41586-018-0567-3](https://doi.org/10.1038/s41586-018-0567-3) · PMC8353605.
Data: GEO GSE96688.

Both from Keji Zhao's lab (NHLBI). The library-construction step of each paper points to an
earlier lab paper instead of describing it: Pico-Seq to Hu *et al.* 2013 *Nat Immunol*
(doi:10.1038/ni.2712, ref. 29), scMNase-seq to Barski *et al.* 2007 *Cell*
(doi:10.1016/j.cell.2007.05.009, ref. 33). Hu *et al.* was fetched and identifies
Illumina Multiplexing Sample Prep Oligonucleotide Kit ref. 1005709. Barski's exact
adapter/index bases remain unavailable.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/scdnase-seq-scmnase-seq__10.1038+nature15740/`, never committed):

| File | What | Used for |
|---|---|---|
| `upstream_scDNase_scMNase.html.txt` | scg_lib_structs page (both methods) | **the only source of oligo sequences**, library diagram, read layout |
| `nature15740_PMC4697938.html.txt` | Pico-Seq paper, PMC author manuscript with Methods | DNase I digestion, carrier DNA, two-step PCR, size selection, sequencer |
| `nature15740_41586_2015_BFnature15740_MOESM35_ESM.pdf.txt` | Pico-Seq Supplementary Information | only the technical-repeat design (one cell's worth split in two); no chemistry |
| `nature15740_..._MOESM36..40_ESM.xls` | Pico-Seq Supplementary Tables | read via a LibreOffice CSV conversion in scratch: GO enrichment tables and a VCF — not chemistry |
| `nature15740_..._MOESM31..34_ESM.ppt` | Pico-Seq supplementary slides | `catppt` returns no text (image-only); not used |
| `s41586-018-0567-3_PMC8353605.html.txt` | scMNase-seq paper (main text; Methods say "see Supplementary Methods") | cell numbers, fragment-size classes |
| `s41586-018-0567-3_41586_2018_567_MOESM1_ESM.pdf.txt` | scMNase-seq Supplementary Information | **all scMNase-seq wet-lab methods** |
| `s41586-018-0567-3_..._MOESM3_ESM.xlsx.txt`, `MOESM4_ESM.xlsx.txt` | per-library QC stats, dataset accessions | library counts; not chemistry |
| `s41586-018-0567-3_..._MOESM2_ESM.pdf.txt` | Reporting Summary | not chemistry |
| Hu *et al.* 2013, doi:10.1038/ni.2712 | Pico-Seq's cited library method | identifies Illumina Multiplexing Sample Prep Oligonucleotide Kit ref. 1005709 |
| Illumina adapter-sequence guide | authoritative vendor document | obsolete Multiplexing-kit adapter and sequencing-primer bases |

Could not be fetched (listed `(manual)` in `MANIFEST.tsv`):

| File | URL | Matters? |
|---|---|---|
| Pico-Seq PMC supplement | https://pmc.ncbi.nlm.nih.gov/articles/instance/4697938/bin/NIHMS723562-supplement-1.pdf | maybe — may duplicate MOESM35, but could hold library details 🔴 |
| scMNase-seq Supplementary Methods (PMC copy) | https://pmc.ncbi.nlm.nih.gov/articles/instance/8353605/bin/NIHMS1727380-supplement-1727380_SupMethod.pdf | probably the same text as MOESM1 |
| scMNase-seq Sup. Tables 1–2 (PMC copy) | https://pmc.ncbi.nlm.nih.gov/articles/instance/8353605/bin/NIHMS1727380-supplement-1727380_SupTab1.xlsx , `..._SupTab2.xlsx` | probably the same as MOESM3/4 |
| Barski *et al.* 2007 (the cited library protocol for scMNase-seq) | https://doi.org/10.1016/j.cell.2007.05.009 | **yes**, same reason |

---

## 1. What it is

Not a barcoding method. One FACS-sorted cell per PCR tube; the chromatin of that cell is
cut by a nuclease **in the tube**, the DNA is purified, and a completely ordinary
end-repair / A-tail / T-overhang-adaptor ligation library is made from it. **The cell's
identity is the sample index on the PCR primer**, one library per cell. No cell barcode,
no UMI, no Tn5. 🟢 (Pico-Seq Methods; scMNase-seq Supplementary Methods) / 🟡 (that the
index is the only cell label — upstream says so; the papers just say "index primers").

The two methods differ only in the nuclease and what it reports 🟢:

| | scDNase-seq / Pico-Seq (2015) | scMNase-seq (2018) |
|---|---|---|
| Nuclease | DNase I — cuts accessible DNA (DHSs) | micrococcal nuclease — eats linker DNA, leaving nucleosome cores (140–180 bp) and sub-nucleosomal particles (≤ 80 bp, TF footprints / accessible sites) |
| Readout | DHS presence per cell | nucleosome positions **and** accessibility in one experiment |
| Carrier DNA | **30 ng circular plasmid** added with the stop buffer | **not mentioned** in the methods 🔴 |
| PCR | two-step: 6 cycles index primers → gel 160–300 bp → 22 cycles P5/P7 | one step: 23 cycles with indexed primers → gel 160–300 bp |
| Cells | 5 NIH3T3 + 14 mESC single cells (38 Pico-Seq libraries incl. pools and FFPE) | 48 NIH3T3, 198 mESC, 278 naïve CD4 T cells |

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Circular carrier DNA** (Pico-Seq): 3 × 10⁷-fold excess of plasmid protects the < 0.001 pg of DHS DNA from loss during extraction and ligation; being circular it has **no ends**, so it cannot take an adaptor and is (mostly) not amplified 🟢 | carrier RNA / glycogen in low-input kits; linear carrier would be ligated and sequenced |
| 2 | **Two-step PCR with gel cut in between** to favour short DHS fragments over carrier background 🟢 | — |
| 3 | Nuclease digestion of a **single cell in its own tube**, then a bulk-style library | later superseded by Tn5 (scATAC) and combinatorial indexing; see [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) for the end-joining alternative |

Even with the carrier trick, only **~2 % of single-cell Pico-Seq reads map** to the
genome — the rest is amplified carrier (1000-cell libraries: ~40 %). 🟢 (Pico-Seq Methods,
data analysis)

## 2. Oligos

For Pico-Seq, 🟢: Hu *et al.* identifies Multiplexing kit ref. 1005709 and Illumina
publishes the obsolete kit sequences. For scMNase-seq, the same layout remains 🟡:
it agrees with upstream scg_lib_structs, but the defining paper and cited Barski method
do not establish the precise indexed oligo set. Notation: `/Phos/` = 5' phosphate,
`*` = phosphorothioate bond, `[i7]` = sample index.

```
Illumina adaptor top              /Phos/ GATCGGAAGAGCACACGTCT
Illumina adaptor bottom           ACACTCTTTCCCTACACGACGCTCTTCCGATCT

Illumina PCR Primer 1.0           AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T
Illumina Multiplexing PCR Primer  CAAGCAGAAGACGGCATACGAGAT[i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T

Illumina P5 adapter               AATGATACGGCGACCACCGAGATCTACAC
Illumina P7 adapter               CAAGCAGAAGACGGCATACGAGAT

TruSeq Read 1 primer              ACACTCTTTCCCTACACGACGCTCTTCCGATCT
TruSeq Read 2 primer              GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
Sample index sequencing primer    GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

What the papers do say: Pico-Seq — "Illumina kits", first PCR with "index primers",
second PCR with "the P5 and P7 primers"; its cited Hu method supplies the missing kit
identity. scMNase-seq says only "Universal adaptors" and "indexed primers". The exact
scMNase-seq index catalogue and length are not given. 🔴

### How they interlock — 🟡 (computed with `lib/illumina`)

- **Adaptor bottom = `illumina.TRUSEQ_READ1`** exactly (33 nt), ending in the 3' **T**
  overhang.
- **Adaptor top = `illumina.INDEX1_PRIMER[:20]`** (20 nt), starting with
  `illumina.STEM` `GATCGGAAGAGC`.
- Annealed, the two oligos pair only over 12 bp: top[1–12] `GATCGGAAGAGC` with the
  bottom's 3' `GCTCTTCCGATC`, leaving a single **3'-T overhang** on the bottom, a
  20-nt unpaired 5' arm on the bottom (`ACACTCTTTCCCTACACGAC`) and an 8-nt unpaired 3' arm
  on the top (`ACACGTCT`). A **Y (forked) adaptor**; upstream draws it this way.
- After T/A ligation each insert strand reads
  `5'- TRUSEQ_READ1 · insert(+A) · adaptor top -3'`. Its 3' end `A·GATCGGAAGAGCACACGTCT`
  is complementary to the **last 21 nt of `TRUSEQ_READ2`** — that is where the
  Multiplexing PCR Primer anneals.
- **PCR Primer 1.0 = `illumina.P5` + `TRUSEQ_READ1[4:]`** (58 nt; P5 ends and Read 1
  starts with the same `ACAC`, so they overlap by 4) — identical to
  `illumina.NEBNEXT_UNIVERSAL_PRIMER`. It is the **same sense** as the ligated strand's
  5' end, so it cannot prime on the ligation product; it primes only on the strand made by
  the Multiplexing primer. (Upstream says the same: "no place to anneal" in round 1.)
- **Multiplexing PCR Primer = `illumina.P7` + i7 + `illumina.TRUSEQ_READ2`** — 58 nt plus
  the index (64 nt with a 6-nt index).
- **Sample index sequencing primer = `illumina.INDEX1_PRIMER` = revcomp(`TRUSEQ_READ2`)
  without its 5' `A`**, and it begins with the 20-nt adaptor top.
- **P5 / P7 adapter** = `illumina.P5` / `illumina.P7` — the second-round primers of
  Pico-Seq's two-step PCR.

## 3. Step by step

### scDNase-seq / Pico-Seq — 🟢 (Pico-Seq Methods)

1. **Sort** single live (DAPI-negative) cells into PCR strip-tubes holding 30 µL lysis
   buffer: 10 mM Tris-HCl pH 7.5, 10 mM NaCl, 3 mM MgCl₂, 0.1 % Triton X-100.
2. **DNase I**: 0.2–1 U (Roche), 37 °C 5 min.
3. **Stop + carrier**: 80 µL stop buffer (10 mM Tris pH 7.5, 10 mM NaCl, 0.15 % SDS,
   10 mM EDTA) carrying 1 µL proteinase K (20 mg/mL) and 5 µL of 6 ng/µL **circular carrier
   DNA** (30 ng). 55 °C 1 h.
4. **Purify**: phenol–chloroform, ethanol precipitation with 20 µg glycogen.
5. **Library** "using Illumina kits as described" (Hu 2013). The paper's Fig. 1a legend
   names end repair, adaptor ligation and PCR; **A-tailing is not stated** but is implied by
   the T-overhang adaptor 🟡 (upstream lists it).
6. **PCR 1, index primers**: 6 × (98 °C 10 s, 67 °C 30 s, 72 °C 30 s). The circular
   carrier is still present; it has no adaptor so is not exponentially amplified. 🟡
7. **Gel**: isolate 160–300 bp on a 2 % E-gel.
8. **PCR 2, P5 + P7 primers**: 22 × (98 °C 10 s, 68 °C 30 s, 72 °C 30 s).
9. **Gel** again, 160–300 bp; sequence on HiSeq 2500.

FFPE variant 🟢: cells scraped from the marked slide area, deparaffinised (Qiagen
solution, 56 °C 3 min), lysis buffer 37 °C 2 h, DNase I as above, then **crosslink reversal
65 °C overnight** before purification and library prep.

### scMNase-seq — 🟢 (Supplementary Methods, MOESM1)

1. **Sort** single live cells into PCR strip-tubes holding 32 µL lysis buffer: 10 mM
   Tris-HCl pH 7.5, 10 mM NaCl, **3 mM CaCl₂** (MNase is Ca²⁺-dependent 🟡), 0.1 % Triton
   X-100.
2. **MNase**: 8 µL MNase diluted 1:3000 in lysis buffer, 37 °C 5 min. (Titration
   experiment: three groups of 10 NIH3T3 single cells at different MNase amounts; the
   Supplementary Methods text names 0.6 and 2.4 U per million cells, and the titration
   sheet of Supplementary Table 1 labels a group "0.1U" — the spacing conclusions held in
   both conditions compared.)
3. **Stop**: 80 µL stop buffer with **10 mM EGTA** (Ca²⁺ chelator, in place of Pico-Seq's
   EDTA) and 1 µL proteinase K (20 mg/mL); 55 °C 1 h.
4. **Purify**: phenol–chloroform, ethanol precipitation with glycogen. **No carrier DNA
   is mentioned.** 🔴 (whether it was omitted or just not written)
5. **Ligate "Universal adaptors"** "as described" (Barski 2007).
6. **PCR, indexed primers**: 23 × (98 °C 10 s, 67 °C 30 s, 72 °C 30 s). No second PCR.
7. **Gel**: 160–300 bp on E-gel; **paired-end** sequencing on HiSeq 2500.

Note on the size window 🟡: a 160–300 bp library minus 122 nt of adaptor/primer
(computed for a 6-nt i7, §4) is an insert of roughly 40–180 bp — exactly the
sub-nucleosomal (≤ 80 bp) to mono-nucleosome (140–180 bp) range scMNase-seq analyses, and
it excludes di-nucleosomes.

## 4. Final library

For Pico-Seq this is computed from the identified vendor kit (🟢 source, 🟡
assembly); for scMNase-seq it remains an upstream-supported reconstruction (🟡):

```
5'- P5 · TruSeq Read 1 (minus its first 4 nt, shared with P5) · <insert> · A · INDEX1_PRIMER (= TruSeq Read 2') · <i7'> · P7' -3'
```

As `lib/` names, top strand:
`illumina.P5` + `illumina.TRUSEQ_READ1[4:]` + insert + `"A"` + `illumina.INDEX1_PRIMER` +
revcomp(i7) + `illumina.P7_RC`. Fixed part 122 nt with a 6-nt i7 (29 + 29 + 34 + 6 + 24).
The `A` before `GATCGG…` is the A-tail; strictly it pairs with the adaptor's 3' T, so the
insert's genomic sequence is bounded by the T on one side and the A on the other. 🟡

Agreement with upstream: identical, including the 6-nt i7 placeholder `NNNNNN`. The
six-base width is defined for the historical Multiplexing kit used by Pico-Seq; for
scMNase-seq it is only the upstream reconstruction. 🟡

## 5. Sequencing

- 🟢 HiSeq 2500 for both. scMNase-seq: **paired-end** (fragment length is the readout).
  Pico-Seq: reads handled "if pair-end sequencing was performed" — so at least some
  Pico-Seq libraries were paired-end, some apparently single-end. 🟢 Read lengths: not
  given. 🔴
- Pico-Seq: 🟢 historical Multiplexing-kit Read 1, Index 1 and Read 2 primers;
  one six-base i7 identifies the library, and there is no i5. scMNase-seq: the same
  three-primer placement is 🟡 because its exact adapter/index set is unresolved.
- 🟢 Mapping detail: unmapped single-cell and low-input Pico-Seq reads were iteratively
  trimmed by 5 bp and re-aligned, until shorter than 26 bp. 🟡 (our reading) This is
  consistent with short DHS inserts read through into the adaptor.

## 6. Agreements and disagreements with upstream

| Upstream claim | Papers | Verdict |
|---|---|---|
| Both methods use "almost the same procedures" | yes — same lysis/stop/purification/gel; nuclease, divalent cation and chelator differ | agrees 🟢 |
| Carrier plasmid is "the trick" for both | stated only for Pico-Seq; absent from scMNase-seq methods | **partly disagrees** 🔴 |
| End repair and A-tailing | Pico-Seq says end repair; A-tailing not stated in either | plausible 🟡 |
| Oligo sequences (TruSeq Y adaptor, PCR Primer 1.0, Multiplexing PCR Primer) | papers name only "Illumina kits", "index primers", "Universal adaptors", "P5 and P7 primers" | consistent, unverifiable 🟡 |
| One PCR with PCR Primer 1.0 + Multiplexing PCR Primer | Pico-Seq used **two** PCRs (6 cycles index primers, gel, 22 cycles P5/P7); scMNase-seq one PCR of 23 cycles | upstream fits scMNase-seq; for Pico-Seq it omits the second, short-primer PCR 🟡 |
| i7 is 6 nt and is the cell barcode | papers: one library per cell with index primers; length not given | plausible 🟡 |
| "Incubate nuclei" | cells are sorted into a Triton-containing lysis buffer and digested there | agrees in effect 🟡 |

## 7. Open questions

- 🔴 Which exact adapter and indexed-primer set scMNase-seq used; Barski 2007 does
  not publish enough sequence detail to establish it.
- 🔴 Was carrier DNA used in scMNase-seq? The 2018 methods do not mention it, while the
  2015 paper presents it as essential.
- 🔴 The identity of the circular carrier plasmid (not named), and whether linearised
  carrier molecules (nicked/broken plasmid) are the source of the ~98 % non-mapping
  single-cell reads.
- 🔴 Pico-Seq second-round "P5 and P7 primers": the short 29/24-nt P5/P7 oligos upstream
  lists, or full-length primers? The upstream final library is the same either way.
- 🔴 Read lengths and read configuration per library.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched both PMC full texts, the Springer supplements, and the
Hu *et al.* paper cited for Pico-Seq's library construction; the Pico-Seq PMC supplement
and the scMNase-seq PMC copies sat behind the download gate.
`tools/doctext.py` made no `.txt` twins for the legacy `.xls` (BIFF) and `.ppt` files;
the `.xls` were read via a LibreOffice CSV conversion in scratch (GO tables and a VCF, no
chemistry), the `.ppt` gave no text. `tools/scrape_primers.py` found the library oligos
in the upstream page (all recognised against `lib/illumina`); the biological papers
contain only unrelated Sanger, ChIP-qPCR and EMSA/reporter oligos. Pico-Seq's promotion
rests on the cited Hu method's exact kit reference plus Illumina's vendor sequence table,
not on sequence text in the paper. scMNase-seq remains 🟡 because that source chain
does not identify its indexed set.
