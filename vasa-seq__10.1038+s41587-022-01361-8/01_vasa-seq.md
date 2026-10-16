# VASA-seq — total RNA from single cells by fragmenting and poly(A)-tailing everything

> **Evidence marking.** 🟢 verbatim from the source (paper or its supplement) · 🟡 derived
> or inferred · 🔴 not published / not available to us. Relationships marked 🟡 *(computed)*
> were worked out with `lib/`; the model is `tools/vasa.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length, kit or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence. An earlier revision of this note carried a whole
> section of sequences copied from that page — a CEL-seq2 barcoded RT primer, RP1/RPI,
> an inDrop v3 bead oligo, a Nextera ligation adapter and aRNA RT primer, and eight read
> primers — and reasoned from them. **None of those sequences has a source we hold**, so
> none of them is reproduced or relied on any more. What the paper itself prints
> (Supplementary Table 12) is marked 🟢; what it names a predecessor for but does not print
> is 🔴, and the model draws it as a bracketed `inferred=True` placeholder.

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
| `GSM5369495_bc2cellID_vasaplate.tsv.gz(.txt)` | Original NCBI GEO GSE176588 VASA-plate barcode-to-cell-ID mapping | **384 unique 8-base plate barcodes**, IDs 001–384; the scraper's default length omits these, so read the table directly |
| `…_41587_2022_1361_MOESM2_ESM.xlsx(.txt)` | Supplementary Table 1 | encapsulation / picoinjection efficiencies (86 %, 98 %) — not chemistry |
| `…_41587_2022_1361_MOESM14_ESM.xlsx(.txt)` | Supplementary Table 13 | cost per cell — not chemistry |
| `upstream_VASA-seq.html(.txt)` | scg_lib_structs "VASA-seq" page | **inspiration only, not a source** — nothing in this note or the model is cited from, or verified against, it |
| `…_MOESM1_ESM.pdf` | Reporting Summary (image-only PDF, 3 pages; text twin is empty) | not used |
| `…_MOESM3–8_ESM.xlsx`, `…MOESM9/10_ESM.tsv`, `…MOESM11/12` | Supplementary Tables 2–11 (markers, DE genes, splicing) | not chemistry |
| `…MOESM16/17` (videos), `…MOESM18_ESM.dxf` (device designs) | Supplementary Videos 1–2, Data 1 | not chemistry |

**The authors' own provenance map.** Supplementary Table 14, section "Oligonucleotides"
(`…MOESM15_ESM.xlsx.txt:79–87`), says exactly which oligo comes from where 🟢. This is the
only authority used for the regions the paper does not print:

| Oligo | Source the authors name | On disk? |
|---|---|---|
| inDrop v3 polyacrylamide bead barcoding | Briggs et al. 2018, doi:10.1126/science.aar5780 | Main article local; original detailed supplement gated 🔴 |
| inDrop bead barcoding reagents (and the IVT sequence) | Zilionis et al. 2017, doi:10.1038/nprot.2016.154 | Original oligo supplement local; main protocol gated 🟡 |
| Cel-seq2 primers for VASA-plate | Muraro et al. 2016, doi:10.1016/j.cels.2016.09.002 — "Indexes used available at GSE176588" | Muraro supplement and **384 exact GEO indexes now local** 🟢; full plate capture oligo still unresolved |
| RA3 and RTP **for VASA-plate** | Hashimshony et al. 2012, doi:10.1016/j.celrep.2012.08.003 | 🔴 no (but see `../cel-seq__10.1016+j.celrep.2012.08.003/01_cel-seq.md`, which traces the same CEL-seq arms to Illumina's TruSeq Small RNA adapter document) |
| RA3, RTP and the 92 dual-index PCR primers **for VASA-drop** | "This study", Table S12 | 🟢 yes — `…MOESM13_ESM.xlsx.txt` |

The GEO index mapping and two predecessor supplements were recovered after the initial
note. The 384 plate barcode members are published; `tools/vasa.py` still draws an
inferred placeholder for the unverified full capture arm. Briggs's v3 bead supplement
and the original Hashimshony 2012 plate adapter source remain unavailable. No
unverified predecessor bases are filled in from a reconstruction.

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

### 2b. What the paper does **not** print — 🔴

No sequence is given here for any of the following, and no reconstruction is substituted.
The role boundaries that *are* established by the paper (Methods, Extended Data Fig. 1b) are
listed because the model draws exactly them and nothing more:

| Region | Established by the paper | Not established 🔴 |
|---|---|---|
| VASA-plate capture oligo | it is a CEL-seq2/SORT-seq primer carrying a T7 promoter, a 6-nt UFI, an **8-nt** cell barcode from a 384-member set (one per well), then anchored poly(T), read in that 5'→3' order 🟢 | every base; the barcode list; the T7/adapter 5' region; whether the poly(T) is T24V |
| VASA-drop bead oligo | an inDrop v3 hydrogel bead oligo carrying barcode 1, a linker, barcode 2, a 6-nt UFI and poly(T); Index 1 spans barcode 1 (8) + a "15 bp common sequence" + i7 (8) 🟢 | every base, including the bead linker between the two barcodes, the barcode whitelists, and the photocleavable/acrydite chemistry |
| VASA-plate terminal PCR arms ("PE1/PE2") | the Methods call them PE1 and PE2 and cite refs 1 and 27 (Hashimshony 2012, Muraro 2016) 🟢 | the two sequences, and which index set the plate libraries used |
| all sequencing primers, both formats | read **lengths** and read **contents** (above) 🟢 | not one primer sequence is printed for either format |

`tools/vasa.py` represents the first two as `capture_primer(fmt)` — a bracketed
`[plate capture region]` / `[drop capture region]` token plus `N`-width UFI and cell-barcode
placeholders and an `inferred=True` oligo-dT — the third as `[plate P5-side arm]` /
`[plate P7-side arm]`, and the fourth by declaring **no** sequencing primer at all
(`SequencingPlan(limitation=…)` and `seqprimers.unavailable_diagram(...)` on the page).

### 2c. What the paper's own tables cannot settle — 🔴

One substantive gap is visible from first-party documents alone, without appeal to any
reconstruction:

> **The printed VASA-drop ligation/RT pair does not interlock with the printed VASA-drop PCR
> arms.** Supplementary Table 12 labels `/5rApp/TGGAATTCTCGGGTGCCAAGG/3SpC3/` and
> `GCCTTGGCACCCGAGAATTCCA` "RA3/RTP oligonucleotide **for VASA-drop**" 🟢, while
> Supplementary Table 14 lists a *separate* "RA3/RTP for VASA-**plate**" credited to
> Hashimshony 2012 🟢. If the drop aRNA really received that RA3 and were primed by that
> RTP, the resulting cDNA would carry `RTP` at its 5' end — and the printed
> `VD_ILMN8_i5` primers end in `TCGTCGGCAGCGTC`, which occurs nowhere in RA3 or RTP on
> either strand 🟡 *(computed with `lib/`)*. So the i5 primer would have no landing site and
> the drop library could not be amplified as described. Either Table 12 lists the plate
> oligos under the drop name, or a further drop-specific adapter/RT primer exists that the
> paper does not print. **We cannot tell which, and this note does not choose.** The model
> accordingly gives the drop library-completion step
> `PCRGeometry.unavailable(..., reason="The VASA-drop bead linker and complete
> sequencing-primer landing sites are not published")`, and the page's step caption says so.

Two smaller points that *are* settled, by the paper against itself:

- **VASA-plate cell barcode is 8 nt, not 6.** The Methods state a 6-nt UFI followed by an
  8-nt cell barcode drawn from 384 well-specific barcodes 🟢, and the 26-cycle Read 1
  closes exactly on 6 + 8 + 12 🟡 *(computed)*. The model uses
  `PLATE_BARCODE_LEN = 8`.
- **Order inside the plate capture oligo** is UFI then cell barcode, from "Read 1 starts
  with the UFI, followed by the cell barcode" 🟢; `vasa._require_order` asserts that order
  on the model rather than leaving it to prose.

## 3. How they interlock — 🟡 (computed with `lib/`)

**VASA-drop index primers** (all 92 checked):

- `VD_ILMN8_i5_N` = **`illumina.P5`** (29) + 8-nt i5 + **`nextera.S5`** (14 nt, =
  `nextera.ADAPTOR_S5[:14]`) — 51 nt, all 46.
- `VD_ILMN8_i7_N` = **`illumina.P7`** (24) + 8-nt i7 + `GGGTGTCGGGTGCAG` — 47 nt, all 46.
  That 15-nt tail is the **"15 bp common sequence"** that Extended Data Fig. 1b places
  between barcode 1 and the i7 inside the 31-cycle Index 1 read 🟢. For the primer to work
  it must anneal at the 5' end of the aRNA, i.e. to the first 15 transcribed bases 🟡
  *(computed)*; its leading `GGG` is the canonical T7 +1…+3 that any class-III T7 promoter
  produces 🟡. The remaining 12 bases are bead-oligo sequence immediately downstream of the
  promoter, and **the bead oligo is not published** 🔴 — so the model keeps them as the
  verbatim `DROP_COMMON` string from Table 12 and does not claim to know the bead context
  they come from.
- All 46 i5 and all 46 i7 indices are distinct within their set. The i7 indices `i7_1..22`
  are the same 8-mers as `i5_25..46`, in the same order (not reverse-complemented). The
  pairing of i5 and i7 primers per sample ("unique dual-indexed") is not given. 🔴
- None of the indices is in `illumina.NEBNEXT_I5_SET1/2` or `NEBNEXT_I7_SET1/2` (either
  orientation); the source set is not named. 🔴

**VASA-drop bead linker and ligation adapter** — 🔴 **not published.** The paper gives the
*order* of the bead's segments (barcode 1 · linker · barcode 2 · UFI · poly(T), via the
31-cycle Index 1 and 14-cycle Read 2 layouts) but no base of the linker, and §2c shows that
the drop ligation/RT chemistry as printed does not close. `tools/vasa.py` therefore draws
`[proprietary bead linker]` as an `inferred=True` bracketed token inside
`drop_library()` and asserts nothing about the ligation junction; the i5 arm it draws is
only what Table 12 prints (`illumina.P5` + i5 + `nextera.S5`). An earlier revision of this
note identified the linker as `nextera.ADAPTOR_S7` and the ligation adapter as
`nextera.INDEX2_PRIMER`; that identification came from the drawn page and has no source, so
it is withdrawn.

**RA3 / RTP interlock** — the one adapter pair the paper does print, checked rather than
eyeballed:

- `RA3` = `illumina.TRUSEQ_SMALL_RNA_RA3` (21 nt) — Table 12's bases are the canonical
  TruSeq Small RNA 3' adapter, so `vasa.py` imports the constant instead of retyping it 🟢.
- `RTP` = `G` + revcomp(`RA3`) (22 nt) 🟡 *(computed)*. `vasa.rtp_segments()` raises if that
  relation ever stops holding, so the 21-nt annealing core is derived from RA3 rather than
  transcribed twice.
- The 5'-rApp / 3'-SpC3 modifications are recorded in
  `vasa.ra3_ligation_geometry()` as a `single-strand`, `5-prime` ligation end with an
  explicit `[3SpC3]` inferred block — not as a duplex overhang, because Table 12 prints a
  single strand 🟢.

**VASA-plate terminal PCR arms** — 🔴 **not published**, cited to Hashimshony 2012 and
Muraro 2016 (neither on disk). `plate_library()` draws `[plate P5-side arm]` and
`[plate P7-side arm]`. The CEL-seq lineage that supplies them is documented in
[the CEL-seq note](../cel-seq__10.1016+j.celrep.2012.08.003/01_cel-seq.md), which traces the
same arms to Illumina's *TruSeq Small RNA* adapter document; that chain is 🟡 for VASA-seq
and is deliberately **not** used to fill in bases on the VASA-seq page.

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
   - plate: NEBNext High-Fidelity 2×, the primer pair the Methods call **PE1/PE2**, cited to
     Hashimshony 2012 and Muraro 2016 and not printed 🔴;
     98 °C 30 s; 7–8 × (98 °C 10 s, 60 °C 30 s, 72 °C 30 s); 72 °C 10 min.
   - drop: Kapa HiFi HotStart, one `VD_ILMN8_i5` + one `VD_ILMN8_i7`; 98 °C 2 min;
     2 × (98 °C 20 s, **55 °C** 30 s, 72 °C 40 s); 5–6 × (98 °C 20 s, **65 °C** 30 s,
     72 °C 40 s); 72 °C 5 min. The two-step annealing fits primers that match the template
     with only their 3' 14–15 nt in the first cycles 🟡 (inferred).
   - Two 0.8× AMPure clean-ups.

## 5. Final libraries — as far as the paper supports, and no further

These are the two constructs `tools/vasa.py` builds. Bracketed regions are
`inferred=True` placeholders: the page draws them dotted and inside `<inf>`, and the
preamble caveat says why. The **order** of the named segments is asserted by
`vasa._require_order` and so cannot silently drift.

**VASA-plate** (`plate_library()`):

```
5'- [plate P5-side arm] · <UFI 6> · <cell barcode 8> · T12 · <fragment cDNA> · RA3 (TGGAATTCTCGGGTGCCAAGG) · [plate P7-side arm] -3'
```

Established: the UFI → cell barcode → poly(T) → cDNA → RA3 order and the 6/8 widths 🟢;
RA3's bases 🟢 (Table 12). Not established: both terminal arms 🔴, the index set 🔴, the
poly(T) length (12 is a drawing width, not a measurement) 🔴. The top strand has the sense
of the aRNA, so `<fragment cDNA>` is **antisense** to the original RNA fragment and Read 2
reads the RNA sense 🟡 *(computed from the T7/oligo-dT geometry)*.

**VASA-drop** (`drop_library()`):

```
5'- P5 · <i5 8> · s5 (TCGTCGGCAGCGTC) · <fragment cDNA> · A… · <UFI 6>' · <bc2 8>' · [proprietary bead linker] · <bc1 8>' · CTGCACCCGACACCC · <i7 8>' · P7' -3'
```

with `'` marking the reverse complement of the segment as written on the bead or primer.
Established: both terminal arms verbatim from the 92 printed `VD_ILMN8` primers 🟢; the
barcode-2 / UFI and barcode-1 / common-sequence / i7 groupings and their widths, from the
14- and 31-cycle read layouts of Extended Data Fig. 1b 🟢. Not established: the bead
linker 🔴, and what lies between the s5 arm and the cDNA — see §2c, where the printed drop
RA3/RTP pair and the printed i5 arm are shown not to interlock 🔴. The barcodes sit on the
**P7 side**, the reverse of VASA-plate, so Read 1 reads the fragment and Read 2 reads
barcode 2 + UFI 🟡 *(computed)*.

## 6. Sequencing

| | VASA-plate | VASA-drop |
|---|---|---|
| Instrument | NextSeq 500, high-output 150 cycles | NovaSeq 6000 S2, 300 cycles |
| Read 1 | 26 cycles: UMI 6 + cell barcode 8 + 12 of the poly(T) | 247 cycles: cDNA |
| Index 1 | 8 cycles: the plate sample index; its width and set are not published 🔴 | 31 cycles: barcode 1 (8) + 15-nt common sequence + i7 (8) |
| Index 2 | — | 8 cycles: i5 |
| Read 2 | 135 cycles: cDNA | 14 cycles: barcode 2 (8) + UMI (6) |

Instruments and cycle numbers 🟢 (Methods); read contents 🟡, except the VASA-plate Read 1 start
(UFI then 8-nt barcode 🟢) and the VASA-drop Index 1 layout (Extended Data Fig. 1b 🟢).
**Read primers: 🔴 the paper prints none, for either format.** No primer is declared in
`vasa.protocol_graph()` and none is drawn: both `SequencingPlan`s carry a `limitation`
string instead, and the page renders `seqprimers.unavailable_diagram(...)` for each
library. Read-content arithmetic 🟡 *(computed)*: 6 + 8 + 12 = 26; 8 + 15 + 8 = 31;
8 + 6 = 14 — the declared cycle counts close exactly on the declared segment widths, which
is what fixes the read contents without a primer sequence. Which strand each read enters
from follows from the index positions in §5 🟡; barcode 1 is therefore read as its reverse
complement relative to the bead 🟡 *(computed)*.

The paper's pipeline merges the four VASA-drop reads into two, with a 16-nt cell barcode
(8 + 8) and a 6-nt UMI, and accepts barcodes within Hamming distance 2 (plate: 1). 🟢 Only
reads with a correct i5/i7 combination are kept (Extended Data Fig. 1b). 🟢

## 7. Open questions

- 🔴 **Which ligation adapter and RT primer VASA-drop really used.** Table 12 says RA3/RTP,
  but the printed i5 arm has no site on an RTP-primed cDNA (§2c). Needs the authors, the
  cited inDrop protocols, or the GEO read structure.
- 🔴 The VASA-plate CEL-seq2/SORT-seq primer sequences and 384 barcodes — cited to Muraro
  2016 with the barcodes in GSE176588; neither fetched (requested in `extra_sources.tsv`).
- 🔴 The inDrop v3 bead oligo and barcode lists — cited to Briggs 2018 and Zilionis 2017;
  neither fetched (requested in `extra_sources.tsv`).
- 🔴 The VASA-plate PE1/PE2 PCR primers — cited to Hashimshony 2012; not fetched
  (requested in `extra_sources.tsv`).
- 🔴 i5/i7 pairing scheme, and why i7_1..22 reuse i5_25..46 8-mers.
- 🔴 IVT conditions and post-RT clean-up for VASA-drop ("as inDrop").
- 🟡 The 5' end of the insert: heat fragments start with 5'-OH, and T4 PNK will also
  phosphorylate them; nothing is ligated to the RNA 5' end in this workflow, so it does not
  matter for the library.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` got the PMC full text and **all 18 Springer supplements**; the PMC copies
of the supplements were saved as `…_supp_…` files but are ~21 KB reCAPTCHA HTML pages, not
the files (the scraper skips them as unreadable). The Reporting Summary PDF is image-only, so
`doctext.py` produced an empty twin. `tools/scrape_primers.py` put the drawn page first and
found the oligo table (`MOESM13`, Supplementary Table 12) second, with all `VD_ILMN8` primers
recognised as `illumina.P5`-containing; it did not flag the RA3/RTP inconsistency — that came
from reading Table 14 and computing against the i5 primers.

`tools/selftest_source_oligos.py vasa-seq__10.1038+s41587-022-01361-8` now checks the two
printed oligo inputs (RA3, RTP) against the extracted text of Supplementary Table 12 and
reports `0 inputs without exact states`: every oligo input in the graph either has a
source-verbatim sequence or is explicitly inferred.
