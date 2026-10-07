# HyDrop-ATAC — droplet scATAC-seq on dissolvable, split-pool-barcoded hydrogel beads

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note; they are not yet asserted in a self-test, because
> this protocol has no `tools/` module yet (status `notes`).

**HyDrop** — De Rop FV, Ismail JN, Bravo González-Blas C, Hulselmans GJ, Flerin CC,
Janssens J, Theunis K, Christiaens VM, Wouters J, Marcassa G, de Wit J, Poovathingal S,
Aerts S. "Hydrop enables droplet-based single-cell ATAC-seq and single-cell RNA-seq using
dissolvable hydrogel beads." *eLife* 11:e73971 (2022).
doi:[10.7554/eLife.73971](https://doi.org/10.7554/eLife.73971) · PMID 35195064 ·
PMC8993220.

Other records for the same method (from `catalogue/scg_lib_structs.tsv`):

- protocols.io "HyDrop-ATAC v1.0 v3",
  doi:[10.17504/protocols.io.b4xvqxn6](https://doi.org/10.17504/protocols.io.b4xvqxn6) —
  the bench protocol (not fetched, see below).
- Bead manufacturing protocol, doi:10.17504/protocols.io.b4cyqsxw (cited in the
  methods; not fetched).
- Cited by upstream for the "semi-suppressive PCR" idea: nanoCAGE / CAGEscan,
  doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470).
- The sister protocol from the same paper and bead chemistry:
  [HyDrop-RNA](../hydrop-rna__10.7554+eLife.73971/01_hydrop-rna.md).

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/hydrop-atac__10.7554+eLife.73971/`, never committed):

| File | What | Used for |
|---|---|---|
| `eLife.73971_PMC8993220.html.txt` | full text (PMC) | methods: bead barcoding, ATAC library prep, sequencing, barcode processing |
| `eLife.73971_supp_elife-73971-supp4.xlsx.txt` | Supplementary file 4, the oligonucleotide list (sheets `Protocol_primers`, `20200130_plate-1-96`, `-2-96`, `-3-96-ATACseq`, `-3-96-RNAseq`) | **every oligo** |
| `eLife.73971_supp_elife-73971-supp1.docx.txt` | Supplementary file 1, "Molecular sequence description of HyDrop bead barcoding" | bead oligo after each extension, one example barcode |
| `eLife.73971_supp_elife-73971-supp2.docx.txt` | Supplementary file 2, "Molecular sequence description of HyDrop-ATAC" | capture, linear PCR, bulk PCR, two sequencing schemes |
| `eLife.73971_supp_elife-73971-supp5.xlsx.txt` | Supplementary file 5, cost sheet | enzyme list and amounts only |
| `eLife.73971_supp_elife-73971-fig1-figsupp2..4-data1.docx.txt` | bead QC tables | not chemistry |
| `upstream_HyDrop_ATAC.html.txt` | scg_lib_structs page (Teichmann lab) | second source, checked below |

Not obtained:

| What | URL | Why it matters |
|---|---|---|
| protocols.io HyDrop-ATAC v1.0 | https://doi.org/10.17504/protocols.io.b4xvqxn6 (also https://www.protocols.io/view/hydrop-atac-v1-0-bxsbpnan) | **bulk PCR cycle number**, volumes; the upstream page also cites a newer oligo sheet `20210712_supp_methods_table_hydrop_oligonucleotide_list.xlsx` from there |
| protocols.io bead protocol | https://doi.org/10.17504/protocols.io.b4cyqsxw | bead synthesis details |
| `elife-73971-data1.zip`, `elife-73971-transrepform1.docx` | https://pmc.ncbi.nlm.nih.gov/articles/instance/8993220/bin/elife-73971-data1.zip | source data; not chemistry |

---

## 1. What it is

A droplet scATAC-seq method built on open-source hardware (a 3-inlet PDMS chip run on an
Onyx pump system). Nuclei are **tagmented in bulk** (standard Nextera Tn5), then
co-encapsulated with a **dissolvable hydrogel bead** carrying a cell-barcoded primer.
DTT in the PCR mix dissolves the bead and releases the primer, and the droplet is
thermocycled: Tn5 is denatured, the gaps are filled, and the bead primer **linearly**
copies each fragment's s7 end, writing the cell barcode onto it. The emulsion is broken
and a bulk PCR adds P5/P7 and a sample index. 🟢 (Results, "Implementation and accuracy
assessment of HyDrop-ATAC"; Methods; Supplementary file 2). That the cleavable bonds are
the BAC cross-linker and the `iThioMC6-D` disulfide on the oligo, and that the gap is
9 nt (standard Tn5), is our reading 🟡.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | Bead primer released by **DTT dissolution** of a BAC-cross-linked polyacrylamide bead (vs. UV release in inDrop) | 10x Gel Beads are also DTT-dissolvable |
| 2 | **Three 96-way split-pool extensions** on the bead (inDrop-style, 3 × 10 nt instead of 2 halves) → 96³ = 884,736 barcodes 🟢 (whitelist size in Methods; 96³ computed 🟡) | inDrop: 2 rounds of 384 |
| 3 | Bead capture site = **s7 (`GTCTCGTGGGCTCGG`)**, the 15-nt Nextera i7-side primer site, so the bead primes directly on tagmented DNA | 10x scATAC captures the s5 side (`nextera.ADAPTOR_S5`-type bead oligo) |
| 4 | Same bead backbone for ATAC and RNA: only the third barcode plate differs (s7 tail vs. UMI + dT) | [HyDrop-RNA](../hydrop-rna__10.7554+eLife.73971/01_hydrop-rna.md) |

Builds on [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md): every insert is flanked
by the 19-bp mosaic end (ME), with s5 or s7 on the outside and the 9-nt gap after
transposition.

## 2. Oligos

🟢 Verbatim from Supplementary file 4 (`supp4.xlsx`), sheet `Protocol_primers` unless
stated. IDT notation as written: `/5Acryd/` = 5' acrydite, `/iThioMC6-D/` = internal
disulfide (thiol-modifier C6 S-S), `*` = phosphorothioate bond, `/56-FAM/` = 5' FAM.

### Bead backbone and the three barcode plates

```
Acrydite_primer   /5Acryd//iThioMC6-D/TTTTTTTTAATACGACTCACTATAGGGAAGCAGTGGTATCAACGCAGAGTAC     (HPLC)

20200130_plate-1-96          GCAGTAGCTG <bc1:10> GTACTCTGCG          96 oligos, 30 nt
20200130_plate-2-96          AGGGTACTCG <bc2:10> GCAGTAGCTG          96 oligos, 30 nt
20200130_plate-3-96-ATACseq  CCGAGCCCACGAGAC <bc3:10> AGGGTACTCG     96 oligos, 35 nt
20200130_plate-3-96-RNAseq   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN <bc3:10> AGGGTACTCG   (RNA beads only, 58 nt)
```

(The plate sheets give only sequences; the `<bc>` layout and the constant flanks were
read off all 96 rows of each sheet 🟡 *(computed)*. First rows, verbatim:
`GCAGTAGCTGTGTAGCAAGTGTACTCTGCG`, `AGGGTACTCGTTAGTTGGACGCAGTAGCTG`,
`CCGAGCCCACGAGACTGACCGTACTAGGGTACTCG`.)

### Final ATAC bead oligo

🟢 Supplementary file 1, with one example barcode:

```
/5Acryd//iThioMC6-D/TTTTTTTTAATACGACTCACTATAGGGAAGCAGTGGTATCAACGCAGAGTACTTCCTGTGAGCAGCTACTGCTCGGACTTATCGAGTACCCTGGCTGAATTAGTCTCGTGGGCTCGG   117 bp
```

### Library PCR primers

```
HYi7_1_CGCTCAGTTC   CAAGCAGAAGACGGCATACGAGATCGCTCAGTTCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C   ... HYi7_16, 69 nt
HYi5_1_TCGTGGAGCG   AATGATACGGCGACCACCGAGATCTACACTCGTGGAGCGTCGTCGGCAGCGTCAGATGTG             ... HYi5_16, 60 nt
HYi5_ATAC_P_1_TCGTGGAGCG  AATGATACGGCGACCACCGAGATCTACACTCGTGGAGCGTCGTCGGCAGCGTCAGATG*T*G  (1-5; "phosphorylated version (not used currently)")
```

Index names 🟢: HYi7 1–16 = `CGCTCAGTTC TATCTGACCT ATATGAGACG CTTATGGAAT TAATCTCGTC
GCGCGATGTT AGAGCACTAG TGCCTTGATC CTACTCAGTC TCGTCTGACT GAACATACGG CCTATGACTC TAATGGCAAG
GTGCCGCTTC CGGCAATGGA GCCGTAACCG`; HYi5 1–16 = `TCGTGGAGCG CTACAAGATA TATAGTAGCT
TGCCTGGTGG ACATTATCCT GTCCACTTGT TGGAACAGTA CCTTGTTAAT GTTGATAGTG ACCAGCGACA CATACACTGT
GTGTGGCGCT ATCACGAAGG CGGCTCTACT GAATGCACGA AAGACTATAG`.

(The "ATAC_P" i5 set is described as "phosphorylated" but the written sequence carries
only phosphorothioates, no `/5Phos/`. 🟡 — a naming slip or an unstated modification.)

### Custom sequencing primers and QC probes

```
HyDrop_CustSeq_R2      CTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC     (HPLC; "custom read 2 read primer for hydrop-rna")
HyDrop_CustSeq_Short   GTACTCTGCGTTGATACCACTGCTTCCGCGGACAG     (HPLC; "custom i7 read primer for hydrop-rna")

Anti-Acrydite_FAM   /56-FAM/TTTTTGTACTCTGCGTTGATACCAC
Anti-ATAC_FAM       /56-FAM/AAAAAACCGAGCCCACGAGAC
Anti-BC1_FAM        /56-FAM/TTTTTCTATCCGTCAGTAC
Anti-BC2_FAM        /56-FAM/TTTTTACACGTTGTGGCAG
Anti-BC3_FAM        /56-FAM/TTTTTCTCCTATCATAGGG
```

The RNA-only oligos in the same sheet (`TSO` `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`, `TSO-P`,
`HYi5_TruSeq_*`) are covered in the HyDrop-RNA note.

### How they interlock — 🟡 (computed with `lib/`)

**Acrydite primer (52 nt after the modifiers)** = `T8` + `TAATACGACTCACTATAGGG` (the T7
promoter, inherited from inDrop; unused in ATAC) + **`rt.SMART_HANDLE`**
(`AAGCAGTGGTATCAACGCAGAGT`) + `AC`. The SMART handle starts at position 28 (1-based).

**Split-pool extensions.** Each plate oligo is a template whose 3' 10 nt pairs with the
current 3' end of the bead oligo; the polymerase copies the barcode and the 5' 10 nt,
which become the next anchor:

- plate 1 3' end `GTACTCTGCG` = revcomp of the bead's 3' `CGCAGAGTAC`. Its 5'
  `GCAGTAGCTG` is copied as linker 1 = **`CAGCTACTGC`**.
- plate 2 3' end `GCAGTAGCTG` pairs with linker 1; its 5' `AGGGTACTCG` becomes linker 2 =
  **`CGAGTACCCT`**.
- plate 3 (ATAC) 3' end `AGGGTACTCG` pairs with linker 2; its 5'
  `CCGAGCCCACGAGAC` = **`nextera.S7_RC`**, copied as **`nextera.S7`** =
  `GTCTCGTGGGCTCGG`, the capture site.

So the bead carries `T8 · T7 · SMART handle · AC · bc1' · L1 · bc2' · L2 · bc3' · s7`,
52 + 10 + 10 + 10 + 10 + 10 + 15 = **117 nt**, matching the length the paper prints. Each
barcode on the bead is the **reverse complement** of the 10 nt written in the plate
oligo. The BC3 sets of the ATAC and RNA third plates are identical (same 96 × 10 nt, same
order).

The example barcode in Supplementary file 1 (`TTCCTGTGAG` / `TCGGACTTAT` / `GGCTGAATTA`)
is **not** in plates 1/2/3 (reverse-complemented or not) — it is illustrative only.

**QC probes**: `Anti-ATAC_FAM` = A6 + `nextera.S7_RC` (hybridises to the capture site);
`Anti-Acrydite_FAM` = T5 + revcomp of the bead's last 20 nt before bc1; the BC1/BC2/BC3
probes match plate-1 row 3, plate-2 row 28 and plate-3 row 19 respectively (one barcode
out of 96, used to measure ~1/96 positive beads).

**HYi7** (all 16) = **`illumina.P7`** + 10-nt index **as named** + `CTGTCCGC` +
`GGAAGCAGTGGTATCAACGCAGAGTAC`. Its last 27 nt are identical to the bead oligo's
`GGAAGCAGTGGTATCAACGCAGAGTAC` (bead positions 26–52: the `GG` of the T7 promoter + SMART
handle + `AC`), so it primes on the complement of the bead backbone. The 8-nt
`CTGTCCGC` is a spacer between i7 and the handle; no reason given 🔴.

**HYi5** (all 16) = **`illumina.P5`** + 10-nt index **as named** +
**`nextera.ADAPTOR_S5[:21]`** (s5 + ME[:7] `AGATGTG`). Same design as the standard Nextera
i5 primers, with 10-nt indices.

**Custom primers**: `HyDrop_CustSeq_Short` is the exact reverse complement of
`HyDrop_CustSeq_R2`; `CustSeq_R2` is the 3' 35 nt of HYi7 (after the index). So `R2`
reads from the bead side into bc1, and `Short` reads the HYi7 index from the other
direction.

**Nextera Index 1 primer** (`nextera.INDEX1_PRIMER` = ME' + s7') ends exactly where bc3
starts, so the barcode read is `bc3 · L2 · bc2 · L1 · bc1`, i.e. **each barcode as written
in its plate oligo**, and 50 nt cover all five segments.

## 3. Step by step

🟢 conditions from Methods, unless marked.

**A. Beads** (made once; Methods "Barcoded hydrogel bead manufacturing")

1. 6 % acrylamide, 0.55 % bisacryloylcystamine (written "bisacryloylcystoylamine"),
   12 µM acrydite primer, 0.6 % APS, in TBSET; 50 µm droplets in HFE-7500 / EA-008,
   65 °C 14 h. Break emulsion, wash (PFO, SPAN-80/hexane, TBSET).
2. Three rounds of split–pool, 96 wells each (Hamilton STAR): 22.5 µL beads + 2.5 µL
   100 µM plate oligo + 25 µL KAPA HiFi HotStart master mix; 95 °C 3 min, 5 × (98 °C 20 s,
   38 °C 4 min, 72 °C 2 min), 98 °C 1 min, 38 °C 10 min, 72 °C 4 min, with vortexing at
   each annealing. Stop with STOP-25 (EDTA), pool, strip the template strand with
   150 mM NaOH / 85 mM Brij-35, neutralise. Round 3 uses the `-3-96-ATACseq` plate for
   ATAC beads.
3. Filter 70 µm, store in glycerol Bead Freezing Buffer at −80 °C.

**B. Library**

4. **Nuclei**: lysis with 0.1 % NP-40, 0.1 % Tween-20, 0.01 % digitonin and **70 µM
   Pitstop** (clathrin inhibitor, "to increase nucleus permeability to Tn5").
5. **Bulk tagmentation**: 50,000 nuclei in 50 µL (10 % DMF, Tris-HCl pH 7.4, 5 mM MgCl₂,
   5 ng/µL Tn5, 70 µM Pitstop, 0.1 % Tween-20, 0.01 % digitonin), 37 °C 1 h. Wash, resuspend
   in 40 µL 0.1 % BSA/PBS. Supplementary file 4 (`reagents` sheet) lists Illumina Tn5 20034198 ("you can use
   your own Tn5"); Supplementary file 5 prices it and notes the cost with in-house Tn5.
6. **Encapsulate** nuclei + 100 µL PCR mix (1.3× Phusion HF buffer, 15 % OptiPrep,
   1.3 mM dNTPs, **39 mM DTT**, 0.065 U/µL Phusion HF, 0.065 U/µL **Deep Vent**,
   0.013 U/µL **ET SSB**) with 35 µL HyDrop-ATAC beads. DTT dissolves the beads in the
   droplet. Why both Phusion and Deep Vent, and what the ET SSB is for, is not stated 🔴.
7. **In-droplet thermocycling** ("linear amplification program"): **72 °C 15 min**,
   98 °C 3 min, then **13 × (98 °C 10 s, 63 °C 30 s, 72 °C 1 min)**.
   - 72 °C: gap fill-in — the 3' ends of each tagmented strand are extended through the
     9-nt gap, ME' and the opposite adapter, so every fragment with an s7 end now has
     `s7'` (`CCGAGCCCACGAGAC`) at a 3' end. Gap fill 🟢 (Supplementary file 2, "Gap Fill
     step"); its assignment to the 72 °C step is upstream's 🟡.
   - 98 °C 3 min releases Tn5 🟡 (Results only says thermocycling denatures Tn5).
   - Cycles: the bead oligo's 3' `s7` anneals to `s7'` and is extended across ME, insert,
     ME' and s5', writing `<bead backbone + barcode> · s7 · ME · insert · ME' · s5'`.
     Only the bead primer is present, so amplification is **linear**. 🟢 (Supplementary
     file 2 shows the bead-primed copies)
8. **Recover**: per 50 µL emulsion, 125 µL 20 % PFO/HFE + 55 µL GITC buffer + 5 µL 1 M DTT,
   ice 5 min; Dynabeads silane clean-up, 1× AMPure.
9. **Bulk index PCR**: 1× KAPA HiFi, 1 µM HYi7 and 1 µM HYi5. HYi7 primes on the
   complement of the bead backbone (its 3' 27 nt), HYi5 on s5' (its 3' 21 nt). **Cycle
   number not given** in the paper 🔴 (protocols.io). Double-sided 0.4×–1.2× AMPure.

Fragment classes 🟡 (inferred; upstream calls the effect "semi-suppressive PCR"):
- s5…s7 fragments → bead-primed at the s7 end, then amplified by HYi5 + HYi7: the
  library.
- s7…s7 fragments → can be barcoded at both ends; HYi7 then primes both ends (identical
  ends → panhandle suppression, and no s5 for HYi5).
- s5…s5 fragments → no s7 site, never barcoded, no HYi7 site.

## 4. Final library — 🟡 (assembled from the oligos above; identical to upstream's drawing)

Top strand, 5'→3':

```
P5 · i5(10, as named) · s5 · ME · <insert> · ME' · s7' · bc3 · L2' · bc2 · L1' · bc1 · revcomp(GGAAGCAGTGGTATCAACGCAGAGTAC) · GCGGACAG · i7'(10) · P7'
```

with `L2'` = `AGGGTACTCG`, `L1'` = `GCAGTAGCTG`, and bc1/bc2/bc3 in the orientation
written in the plate oligos. The 3' tail after bc1 is `GTACTCTGCGTTGATACCACTGCTTCC`
followed by `GCGGACAG` (= revcomp of the HYi7 spacer region). The T8 and the rest of the
T7 promoter on the bead are **not** in the library: HYi7 primes inside the bead backbone.
Fixed sequence outside the insert: 225 nt 🟡 *(computed)*.

## 5. Sequencing

🟢 NextSeq 500 / NextSeq 2000: **Read 1 50** (ATAC mate 1), **Index 1 52** (cell barcode),
**Index 2 10** (sample index), **Read 2 50** (ATAC mate 2). Standard Nextera primers 🟡
(implied by the paper's "standard primers" scheme in Supplementary file 2; upstream
names them):

| Read | Primer | Reads |
|---|---|---|
| Read 1 | `nextera.READ1_PRIMER` (s5 + ME) | genomic insert, from the s5 end |
| Index 1 | `nextera.INDEX1_PRIMER` (ME' + s7') | `bc3 · L2' · bc2 · L1' · bc1` = 50 nt, then 2 nt of `GTACTCTGCG` 🟡 |
| Index 2 | standard i5 read | HYi5 sample index |
| Read 2 | `nextera.READ2_PRIMER` (s7 + ME) | genomic insert, from the s7 end |

The **HYi7 index is not read** in this scheme — only the i5 index demultiplexes samples
(upstream marks it "deprecated, not sequenced"; consistent with the paper's cycle list
🟡). The 52-cycle Index 1 needs a custom recipe on some instruments 🟢 (Supplementary
file 2). Barcode processing 🟢: the inter-barcode linkers are trimmed from the index read and the
result is matched against the 884,736-entry whitelist with ≤ 1 mismatch (the trimmed
barcode is 3 × 10 = 30 nt 🟡).

Supplementary file 2 also shows a **deprecated custom-primer scheme** 🟢: Read 1 with
the Nextera read 1 primer, a "HyDrop custom primer" Read 2 starting in the bead handle
(`CTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC` = `HyDrop_CustSeq_R2`, reads bc1 → bc3), the I7
index read with `GTACTCTGCGTTGATACCACTGCTTCCGCGGACAG` (= `HyDrop_CustSeq_Short`), and the i5
read. In that scheme the second ATAC mate is not sequenced (single-ended) 🟡.

## 6. Checked against upstream (scg_lib_structs)

Agreements 🟡 *(computed)*: the acrydite primer, the three plate-oligo layouts, the final
bead oligo layout, HYi5 and HYi7 layouts, the four sequencing primers, the final library
string (equal character-for-character after replacing indices/barcodes by N), and the
read configuration (50 / 52 / 10 / 50; i7 sample index not read).

Differences:
- Upstream writes HYi7 without the 3' phosphorothioates (`...CAGAGTAC`); Supplementary
  file 4 has `...CAGAGT*A*C`. 🟢 (both sources)
- Upstream draws the bead annealing to the s7' strand of a fragment drawn s5-top; the
  paper's Supplementary file 2 draws the s7-top strand. Same chemistry. 🟡
- Upstream names its oligo source as a later protocols.io sheet
  (`20210712_...oligonucleotide_list.xlsx`), not Supplementary file 4; not fetched, so a
  later sequence change cannot be excluded 🔴.
- Upstream does not mention the T7 promoter in the bead backbone, the
  `HYi5_ATAC_P` set, or the custom primers. 🟢

## 7. Open questions

- 🔴 Bulk PCR cycle number and annealing temperature (protocols.io only).
- 🔴 Why Phusion **and** Deep Vent, and ET SSB, in the droplet mix.
- 🔴 Purpose of the 8-nt `CTGTCCGC` spacer in HYi7, and why HYi7 carries an index that
  is not sequenced in the published scheme.
- 🟡 Whether the "ATAC_P" i5 primers carry an unstated 5' phosphate.
- 🟡 Supplementary file 1 writes the RNA third extension with 9 N, while the RNA plate-3
  oligos carry 8 N (`NNNNNNNN`) — relevant to HyDrop-RNA, not ATAC.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text and JATS XML; the Europe PMC
supplementary zip stopped at exactly 50,000,000 bytes and does not open as a zip; the
supplementary DOCX/XLSX files are on disk with `.txt` twins although `MANIFEST.tsv`
still lists them as `(manual)` (supp5 came with a reCAPTCHA page saved as
`recaptcha_page_supp5.html.bad`). `tools/scrape_primers.py --max-hits 80` put
Supplementary file 4 first and recognised the SMART handle, P5/P7 and s5/s7; the 80-hit
cap cut most of the 288 plate oligos, which were analysed directly from the sheet. The
reaction order came from reading the methods and Supplementary file 2.
