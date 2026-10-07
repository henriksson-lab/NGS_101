# VASA-seq — total RNA from single cells by fragmenting and poly(A)-tailing everything

> **Evidence marking.** 🟢 verbatim from the source (paper or its supplement) · 🟡 derived,
> inferred, or only in a secondary source (the upstream scg_lib_structs page) · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out with
> `lib/` while writing this note; they are not yet asserted in a self-test, because this
> protocol has no `tools/` module yet (status `notes`).

**VASA-seq** — Salmen F, De Jonghe J, Kaminski TS, Alemany A, Parada GE, Verity-Legg J,
Yanagida A, Kohler TN, Battich N, van den Brekel F, Ellermann AL, Martinez Arias A,
Nichols J, Hemberg M, Hollfelder F, van Oudenaarden A. "High-throughput total RNA
sequencing in single cells using VASA-seq." *Nature Biotechnology* 40 (2022) 1780–1793.
doi:[10.1038/s41587-022-01361-8](https://doi.org/10.1038/s41587-022-01361-8), PMID 35760914,
PMC9750877 (open access). Data: GEO GSE176588. The catalogue lists no other paper for this
method.

Two formats in one paper, both covered here: **VASA-plate** (FACS into 384-well plates,
CEL-seq2/SORT-seq barcoding and library) and **VASA-drop** (inDrop v3 hydrogel beads plus two
picoinjections, Nextera-style library).

## Sources

Fetched by `tools/get_sources.py` into
`_data/sources/vasa-seq__10.1038+s41587-022-01361-8/`, never committed:

| File | What | Used for |
|---|---|---|
| `s41587-022-01361-8_PMC9750877.html(.txt)` | paper, full text (PMC) | **all methods** (lines 1506–1760 of the `.txt`), read lengths, barcode lengths |
| `…_41587_2022_1361_MOESM13_ESM.xlsx(.txt)` | Supplementary Table 12, "oligonucleotides not found in cited publications" | **RA3, RTP, 46 + 46 VASA-drop index primers, 195 rRNA-depletion probes** |
| `…_41587_2022_1361_MOESM15_ESM.xlsx(.txt)` | Supplementary Table 14, reagents and equipment | which oligo comes from which publication; enzymes and kits |
| `…_41587_2022_1361_MOESM2_ESM.xlsx(.txt)` | Supplementary Table 1 | encapsulation / picoinjection efficiencies (86 %, 98 %) — not chemistry |
| `…_41587_2022_1361_MOESM14_ESM.xlsx(.txt)` | Supplementary Table 13 | cost per cell — not chemistry |
| `upstream_VASA-seq.html(.txt)` | scg_lib_structs "VASA-seq" page | the **CEL-seq2 RT primer, inDrop v3 bead oligo, TruSeq small-RNA PCR primers, VASA-drop ligation adapter and RT primer, read primers** — none of these sequences is in the paper |
| `…_MOESM1_ESM.pdf` | Reporting Summary (image-only PDF, 3 pages; text twin is empty) | not used |
| `…_MOESM3–8_ESM.xlsx`, `…MOESM9/10_ESM.tsv`, `…MOESM11/12` | Supplementary Tables 2–11 (markers, DE genes, splicing) | not chemistry |
| `…MOESM16/17` (videos), `…MOESM18_ESM.dxf` (device designs) | Supplementary Videos 1–2, Data 1 | not chemistry |

Not available to us:

| What | Where | Would settle |
|---|---|---|
| inDrop bead-barcoding protocol (Zilionis et al. 2017, *Nat Protoc*, doi:10.1038/nprot.2016.154) and inDrop v3 design (Briggs et al. 2018, *Science*, doi:10.1126/science.aar5780) | cited for the bead oligo and the IVT steps | the bead oligo and inDrop v3 barcode lists 🟢 (we have them only from upstream 🟡); the post-RT clean-up and IVT for VASA-drop |
| SORT-seq (Muraro et al. 2016, *Cell Syst*, doi:10.1016/j.cels.2016.09.002) | cited for the VASA-plate RT primers, pooling and IVT | the 384 CEL-seq2 primers (paper says the barcodes are in GSE176588) |
| CEL-seq (Hashimshony et al. 2012) | cited for VASA-plate RA3/RTP and PE1/PE2 | see `../cel-seq-family__10.1016+j.celrep.2012.08.003/01_cel-seq.md` |

---

## 1. What it is

VASA-seq turns *every* RNA in the lysate — mRNA, unspliced pre-mRNA, lncRNA, histone mRNA,
small RNA — into something oligo-dT can prime from. The single-cell lysate is **fragmented by
heat**, the fragment ends are **repaired with T4 PNK**, and **E. coli poly(A) polymerase** adds a
poly(A) tail to every fragment. From there on the chemistry is a standard 3'-barcoding,
T7-amplifying (CEL-seq-type) workflow: barcoded oligo-dT RT, second strand, IVT. Because the
amplified RNA is now mostly rRNA, **rRNA is removed after IVT** by probe hybridisation and
RNase H, and the library is finished by **3' adapter ligation, RT and PCR** as in CEL-seq.
🟢 (main text, Fig. 1a; Methods)

Each poly(A)-tailed fragment yields one molecule with a UMI, which the paper calls a
**unique fragment identifier (UFI)** — counting fragments rather than transcripts, strand
specific. 🟢

| | New thing here | Builds on |
|---|---|---|
| 1 | **Fragment, then poly(A)-tail all RNA** in the cell lysate, so random internal fragments and non-polyadenylated RNA get oligo-dT priming | earlier single-cell total-RNA methods (the paper benchmarks against Smart-seq-total, refs 16–18); here poly(A)-tailing is applied to the whole fragmented transcriptome |
| 2 | **rRNA depletion on the amplified RNA (aRNA)** with 195 DNA probes + thermostable RNase H, then DNase | probe/RNase H depletion is common for bulk RNA-seq; doing it after IVT means it is done once per pool, not per cell |
| 3 | **Two picoinjections** into each droplet (repair/poly(A), then RT) after bead–cell co-encapsulation (VASA-drop) | inDrop v3 beads and encapsulation; picoinjection (Abate 2010) |
| 4 | VASA-drop: a **custom i7 primer** that adds a second 8-nt index behind inDrop's barcode 1 so that both index reads can flag index hopping on the NovaSeq | inDrop v3 (where Index 1 reads only barcode 1) |

The "New thing" column follows the paper 🟢; the "Builds on" column is our framing 🟡.

Concept pages: [reverse transcription](../ref/concepts/reverse-transcription.md) (oligo-dT
priming), [small-RNA ligation](../ref/concepts/small-rna-ligation.md) (the 5'-adenylated,
3'-blocked adapter ligated to the aRNA 3' end by T4 RNA ligase 2 truncated). There is no
template switching and no tagmentation: the Nextera sequences in VASA-drop come from the bead
oligo and from the ligated adapter, not from Tn5.

## 2. Oligos

### 2a. From the paper (Supplementary Table 12, `MOESM13_ESM.xlsx`) — 🟢

Written exactly as the table gives them (IDT-style names without leading slash):

```
RTP oligonucleotide for VASA-drop   GCCTTGGCACCCGAGAATTCCA
RA3 oligonucleotide for VASA-drop   5rApp/TGGAATTCTCGGGTGCCAAGG/3SpC3

VD_ILMN8_i5_1    AATGATACGGCGACCACCGAGATCTACACAGCGCTAGTCGTCGGCAGCGTC     (… _i5_46; 51 nt)
VD_ILMN8_i5_46   AATGATACGGCGACCACCGAGATCTACACGGAACGTTTCGTCGGCAGCGTC
VD_ILMN8_i7_1    CAAGCAGAAGACGGCATACGAGATAACCGCGGGGGTGTCGGGTGCAG         (… _i7_46; 47 nt)
VD_ILMN8_i7_46   CAAGCAGAAGACGGCATACGAGATTGCATTGCGGGTGTCGGGTGCAG

template:
VD_ILMN8_i5_N    AATGATACGGCGACCACCGAGATCTACAC <i5> TCGTCGGCAGCGTC
VD_ILMN8_i7_N    CAAGCAGAAGACGGCATACGAGAT <i7> GGGTGTCGGGTGCAG
```

plus 195 rRNA-depletion probes `Ribo_Dep_aRNA_<rRNA>_<n>` (mostly 50 nt; 18S ×38, 28S ×101,
5.8S ×3, 5S ×3, `_12_` ×19 and `_16_` ×31 — presumably mitochondrial 12S/16S 🟡; counts and lengths 🟢 from the table, 188 of 195 are 50 nt), described as the
reverse complement of published probes (ref. 62). Example 🟢:
`Ribo_Dep_aRNA_18_1  GTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA`. Being the
reverse complement of probes designed against rRNA, they are **sense to the rRNA** — which
is what is needed, because the aRNA in both formats is antisense. 🟡

The **CEL-seq2/SORT-seq RT primers**, the **inDrop v3 bead oligo** and the **VASA-plate
PCR primers** are not in the paper: Supplementary Table 14 attributes the first two to Muraro 2016
and Briggs 2018; for the plate PE1/PE2 PCR primer mix the Methods cite refs 1 and 27
(Hashimshony 2012, Muraro 2016). 🟢 (attribution) / 🔴 (sequences)

### 2b. From the upstream page only — 🟡

VASA-plate (CEL-seq2 / TruSeq Small RNA; same oligos as in the CEL-seq note):

```
Barcoded RT primer   GCCGGTAATACGACTCACTATAGGGAGTTCTACAGTCCGACGATC[6-bp UMI][6-bp cell barcode]TTTTTTTTTTTTTTTTTTTTTTTTV
RA3 primer           /5rApp/ TGGAATTCTCGGGTGCCAAGG /3SpC3/
RTP primer           GCCTTGGCACCCGAGAATTCCA
RP1                  AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA
RPI                  CAAGCAGAAGACGGCATACGAGAT[6-bp RPI]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
Read 1 seq primer    GTTCAGAGTTCTACAGTCCGACGATC
Index read primer    TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC
Read 2 seq primer    GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
```

VASA-drop (the upstream page thanks J. De Jonghe, an author, for these):

```
Beads-oligo-dT19V    /5Acryd/iSpPC/CGATGACGTAATACGACTCACTATAGGGTGTCGGGTGCAG[8-bp barcode1]GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG[8-bp barcode2][6-bp UMI]TTTTTTTTTTTTTTTTTTV
Ligation adapter     /5rApp/ CTGTCTCTTATACACATCTGACGCTGCCGACGA /3SpC3/
aRNA RT primer       TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Library PCR P5       AATGATACGGCGACCACCGAGATCTACAC[8-bp i5 index]TCGTCGGCAGCGTC
Library PCR P7       CAAGCAGAAGACGGCATACGAGAT[8-bp i7 index]GGGTGTCGGGTGCAG
Read 1 seq primer    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Index 1 seq primer   CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Index 2 seq primer   CTGTCTCTTATACACATCTGACGCTGCCGACGA
Read 2 seq primer    GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
```

### 2c. Upstream vs paper — agreements and disagreements

- **Agree** 🟢/🟡: the VASA-drop library PCR primers. Upstream's two templates match all
  46 + 46 `VD_ILMN8` primers of Table 12 base for base outside the index (computed below).
- **Agree**: VASA-plate RA3 and RTP — upstream's sequences are identical to Table 12's.
- **Disagree — the VASA-drop ligation adapter and RT primer.** Table 12 (and Table 14) name
  the TruSeq small-RNA `RA3` / `RTP` as the oligos "for VASA-drop", and Table 14 lists
  separate "RA3/RTP for VASA-plate" from Hashimshony 2012. The upstream page instead gives a
  Nextera-type ligation adapter and RT primer for VASA-drop. **The paper's own VASA-drop PCR
  primers support upstream**: the i5 primers end in `TCGTCGGCAGCGTC` (Nextera s5), which has
  no match in RA3 or RTP on either strand, so an RTP-primed cDNA could not be amplified by
  them 🟡 *(computed)*. Most likely Table 12 printed the plate oligos under the drop name,
  and the true VASA-drop "RA3"/"RTP" are upstream's adapter and RT primer. 🔴 not confirmed
  in the paper.
- **Disagree — VASA-plate cell barcode length.** Upstream writes a 6-bp cell barcode (the
  original CEL-seq2 design). The paper says Read 1 has a 6-nt UFI/UMI followed by an
  **8-nt** cell barcode, 384 barcodes, one per well (SORT-seq set) 🟢; 6 + 8 + 12 = 26 is
  also the Read 1 length used 🟡. Upstream's `[6-bp cell barcode]` is wrong for VASA-plate.
- **Naming**: upstream calls the bead primer "oligo-dT19V" — it has 18 T plus V, 19 nt
  🟡 *(computed)*; the CEL-seq2 primer has T24V.
- **Order of UMI and cell barcode in the plate primer**: upstream writes `[UMI][cell barcode]`
  5'→3'; the paper's "Read 1 starts with the UFI, followed by the cell barcode" agrees. 🟢

## 3. How they interlock — 🟡 (computed with `lib/`)

**VASA-drop index primers** (all 92 checked):

- `VD_ILMN8_i5_N` = **`illumina.P5`** (29) + 8-nt i5 + **`nextera.S5`** (14 nt, =
  `nextera.ADAPTOR_S5[:14]`) — 51 nt, all 46.
- `VD_ILMN8_i7_N` = **`illumina.P7`** (24) + 8-nt i7 + `GGGTGTCGGGTGCAG` — 47 nt, all 46.
  The 15-nt tail is the **3' end of the T7 promoter plus the next 12 nt of the inDrop v3
  bead oligo** (`…CACTATA|GGGTGTCGGGTGCAG`, starting at the T7 +1 G). It anneals to the first
  15 nt transcribed by T7, i.e. the 5' end of the aRNA, and is the "15 bp common sequence"
  of Extended Data Fig. 1b.
- All 46 i5 and all 46 i7 indices are distinct within their set. The i7 indices `i7_1..22`
  are the same 8-mers as `i5_25..46`, in the same order (not reverse-complemented). The
  pairing of i5 and i7 primers per sample ("unique dual-indexed") is not given. 🔴
- None of the indices is in `illumina.NEBNEXT_I5_SET1/2` or `NEBNEXT_I7_SET1/2` (either
  orientation); the source set is not named. 🔴

**VASA-drop bead and library adapters** (upstream sequences):

- The bead linker between barcode 1 and barcode 2 is exactly **`nextera.ADAPTOR_S7`**
  (s7 + ME, 34 nt).
- The ligation adapter is exactly **revcomp(`nextera.ADAPTOR_S5`)** = `nextera.INDEX2_PRIMER`
  (33 nt); ligated to the aRNA 3' end, it gives the aRNA a 3' end that the aRNA RT primer
  (= `nextera.ADAPTOR_S5`) anneals to. The i5 PCR primer's s5 tail then primes on that
  cDNA's 5' end.
- Index 1 primer = `nextera.INDEX1_PRIMER` = revcomp(`ADAPTOR_S7`); Read 2 primer =
  `ADAPTOR_S7`; Read 1 primer = `ADAPTOR_S5`. All standard Nextera.

**VASA-plate** (CEL-seq / TruSeq Small RNA, as in the CEL-seq note):

- RT primer = 5 nt + **T7 promoter** (`TAATACGACTCACTATAGGG`, from position 5) + the last
  20 nt of RA5 (`AGTTCTACAGTCCGACGATC`) + UMI + barcode + T24V.
- RP1 = **`illumina.P5`** + `GTTCAGAGTTCTACAGTCCGA` (the first 21 nt of RA5, 50 nt); its
  3' 15 nt `AGTTCTACAGTCCGA` are in the RT primer, so it primes on the cDNA copy of the
  aRNA 5' end.
- RTP = `G` + revcomp(RA3) (22 nt).
- RPI tail = `GTGACTGGAGTT` (first 12 nt of `illumina.TRUSEQ_READ2`) + revcomp(RA3);
  the Read 2 primer is this tail and the Index read primer is its reverse complement.

## 4. Step by step

Volumes, concentrations and conditions 🟢 from the Methods unless marked.

### Shared front end (per cell)

1. **Lyse** with IGEPAL CA-630 + thermolabile proteinase K in First-Strand Buffer (drop mix
   also 3 mM MgCl₂): plate 25 °C 1 h, then 55 °C 10 min (heat-inactivates the thermolabile
   proteinase 🟡 role inferred); drop 20 min at room temperature (23 °C).
2. **Fragment RNA by heat** in that buffer: 85 °C (plate 3 min; drop 6 min 30 s), then onto ice.
   Divalent-cation-catalysed hydrolysis leaves 5'-OH and 2',3'-cyclic phosphate ends 🟡
   (standard chemistry; the paper does not say).
3. **Repair and poly(A)-tail** together, 37 °C: T4 PNK, E. coli poly(A) polymerase, ATP,
   DTT, RNaseOUT (plate 1 h; drop 25 min at 23 °C then 8 min at 37 °C). PNK removes the 3'
   phosphate so that the poly(A) polymerase has a 3'-OH to extend 🟡 (inferred role).
4. **RT** with SuperScript III and dNTPs, 50 °C (plate 1 h; drop 2 h, then 70 °C 20 min),
   primed by the **barcoded oligo-dT** (plate: CEL-seq2 primer already in the well; drop:
   bead oligo, released by UV photocleavage of the PC spacer for 7 min before fragmentation).
   First strand 🟡: `5'-[T7 promoter]-[adapter]-[UMI/barcodes]-T…-<cDNA of fragment> -3'`.

### VASA-plate specifics

- Cells FACS-sorted into 384-well plates holding 5 µL mineral oil and 50 nL CEL-seq2/SORT-seq
  primer (0.25 µM); all reagents dispensed at 50 nL (NanoDrop II), then 1,100 nL
  second-strand mix.
- **Second strand**: E. coli DNA Pol I + RNase H, 16 °C 2 h, 85 °C 20 min.
- Pool and **IVT** "as SORT-seq" (conditions not given here 🔴); then ExoSAP-IT 37 °C 15 min.

### VASA-drop specifics

- Co-encapsulation of cells (450/µL, 15 % OptiPrep), inDrop v3 beads and lysis mix in
  ~0.55 nL droplets; fluorinated oil with 5 % 008-FluoroSurfactant.
- **Picoinjection 1** (repair + poly(A) mix), **picoinjection 2** (RT mix), each by
  electro-coalescence (250 V, 10 kHz); collected in fractions of ~1,000 cells.
- De-emulsify with perfluoro-octanol. Downstream "up to and including IVT" as in inDrop
  (Zilionis 2017) 🟢 — the conditions are not in this paper 🔴. Table 14 lists ExoI,
  HinfI, NEBNext second-strand module and T7 kits (MEGAscript, HiScribe), which suggests the
  inDrop sequence of primer digestion, second strand and IVT. 🟡

### Shared back end (per pool)

5. **Clean up aRNA**: AMPure XP 1.8×; adjust to ≤ 100 ng/µL.
6. **rRNA depletion**: 6 µL aRNA + 2 µL probes (25 µM) + 2 µL hybridisation buffer
   (500 mM Tris-HCl pH 7.5, 1 M NaCl); 95 °C 2 min, ramp 0.1 °C/s to 45 °C; add thermostable
   RNase H + buffer (with 50 mM MgCl₂), 45 °C 30 min. Then RQ1 DNase + CaCl₂, 37 °C 30 min
   (removes the DNA probes 🟡 inferred role). AMPure 1.6×.
7. **3' ligation**: 1 µL "RA3" (20 µM, Table 12) to 5 µL aRNA, 70 °C 2 min, ice; T4 RNA
   ligase 2 truncated, 25 °C 1 h. The adapter is 5'-adenylated and 3'-blocked, so it joins
   only to aRNA 3' ends and cannot concatemerise ([small-RNA ligation](../ref/concepts/small-rna-ligation.md)) 🟡.
   (For VASA-drop, which adapter — see §2c.)
8. **RT of the aRNA**: "RTP" (20 µM, 2 µL) + dNTPs, 65 °C 5 min; SuperScript III, 50 °C 1 h,
   70 °C 15 min. **RNase A** 37 °C 30 min, AMPure 1×, elute 20 µL.
9. **PCR** on 10 µL:
   - plate: NEBNext High-Fidelity 2×, PE1/PE2 (RP1 + RPI, 🟡 identities from upstream);
     98 °C 30 s; 7–8 × (98 °C 10 s, 60 °C 30 s, 72 °C 30 s); 72 °C 10 min.
   - drop: Kapa HiFi HotStart, one `VD_ILMN8_i5` + one `VD_ILMN8_i7`; 98 °C 2 min;
     2 × (98 °C 20 s, **55 °C** 30 s, 72 °C 40 s); 5–6 × (98 °C 20 s, **65 °C** 30 s,
     72 °C 40 s); 72 °C 5 min. The two-step annealing fits primers that match the template
     with only their 3' 14–15 nt in the first cycles 🟡 (inferred).
   - Two 0.8× AMPure clean-ups.

## 5. Final libraries — 🟡 (assembled from the oligos above)

Top strand 5'→3' (the strand that starts with P5).

**VASA-plate**:

```
5'- P5 · RA5[:21] (GTTCAGAGTTCTACAGTCCGA) · CGATC · <UMI 6> · <cell barcode 8> · T24 · V · <insert> · RA3 (TGGAATTCTCGGGTGCCAAGG) · AACTCCAGTCAC · <rpi'> · P7' -3'
```

The top strand has the sense of the aRNA, so `<insert>` here is **antisense** to the
original RNA fragment: the P5 side carries UMI, barcode and poly(T), the P7 side the RA3
ligated at what was the fragment's 5' end. Read 2 (cDNA) therefore reads the RNA sense. 🟡

**VASA-drop**:

```
5'- P5 · <i5> · ADAPTOR_S5 (s5 + ME) · <insert> · A… · <UMI 6>' · <bc2 8>' · revcomp(ADAPTOR_S7) · <bc1 8>' · CTGCACCCGACACCC · <i7>' · P7' -3'
```

with `'` marking the reverse complement of the segment as written on the bead or primer.
The top strand is the cDNA made from the s5 RT primer on the aRNA, so `<insert>` has the
**sense** of the original RNA fragment, followed by the copy of its poly(A) tail. The
barcodes sit on the **P7 side**, the reverse of VASA-plate, and Read 1 (cDNA) reads the RNA
sense. 🟡 *(computed; assumes upstream's Nextera ligation adapter, §2c)*

## 6. Sequencing

| | VASA-plate | VASA-drop |
|---|---|---|
| Instrument | NextSeq 500, high-output 150 cycles | NovaSeq 6000 S2, 300 cycles |
| Read 1 | 26 cycles: UMI 6 + cell barcode 8 + 12 of the poly(T) | 247 cycles: cDNA |
| Index 1 | 8 cycles: RPI index (6 nt per upstream 🟡) | 31 cycles: barcode 1 (8) + 15-nt common sequence + i7 (8) |
| Index 2 | — | 8 cycles: i5 |
| Read 2 | 135 cycles: cDNA | 14 cycles: barcode 2 (8) + UMI (6) |

Instruments and cycle numbers 🟢 (Methods); read contents 🟡, except the VASA-plate Read 1 start
(UFI then 8-nt barcode 🟢) and the VASA-drop Index 1 layout (Extended Data Fig. 1b 🟢).
Read primers 🟡 (upstream): plate — small-RNA Read 1 (full RA5), RPI-tail Read 2 and its
reverse complement for the index; drop — standard Nextera Read 1 / Index 1 / Index 2 /
Read 2 (§3). Read-content arithmetic 🟡 *(computed)*: 6 + 8 + 12 = 26; 8 + 15 + 8 = 31;
8 + 6 = 14. Because Read 2 starts with `ADAPTOR_S7`, which is the bead's own linker, it
reads barcode 2 and the UMI **as written on the bead**; Index 1 primes on the other strand
and reads barcode 1 as its **reverse complement**. 🟡

The paper's pipeline merges the four VASA-drop reads into two, with a 16-nt cell barcode
(8 + 8) and a 6-nt UMI, and accepts barcodes within Hamming distance 2 (plate: 1). 🟢 Only
reads with a correct i5/i7 combination are kept (Extended Data Fig. 1b). 🟢

## 7. Open questions

- 🔴 **Which ligation adapter and RT primer VASA-drop really used.** Table 12 says RA3/RTP;
  the drop PCR primers and the upstream page say Nextera s5-based (§2c). Needs the authors
  or the GEO read structure (does VASA-drop Read 1 start with ME?).
- 🔴 The VASA-plate CEL-seq2/SORT-seq primer sequences and 384 barcodes (in GSE176588 /
  Muraro 2016, not fetched).
- 🔴 The inDrop v3 bead oligo and barcode lists (Zilionis 2017 / Briggs 2018); upstream only.
- 🔴 i5/i7 pairing scheme, and why i7_1..22 reuse i5_25..46 8-mers.
- 🔴 IVT conditions and post-RT clean-up for VASA-drop ("as inDrop").
- 🟡 The 5' end of the insert: heat fragments start with 5'-OH, and T4 PNK will also
  phosphorylate them; nothing is ligated to the RNA 5' end in this workflow, so it does not
  matter for the library.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC full text and **all 18 Springer supplements**; the PMC copies
of the supplements were saved as `…_supp_…` files but are ~21 KB reCAPTCHA HTML pages, not
the files (the scraper skips them as unreadable). The Reporting Summary PDF is image-only, so
`doctext.py` produced an empty twin. `tools/scrape_primers.py` put the upstream page first and
found the oligo table (`MOESM13`, Supplementary Table 12) second, with all `VD_ILMN8` primers
recognised as `illumina.P5`-containing; it did not flag the RA3/RTP mislabelling — that came
from reading Table 14 and computing against the i5 primers.
