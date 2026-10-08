# CH-ATAC-seq — combinatorial-hybridization single-cell ATAC-seq

> **Evidence marking.** 🟢 verbatim from a source read for this note · 🟡 derived,
> inferred, or taken only from the secondary upstream page · 🔴 not published / not
> available to us. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test (no `tools/` module yet).
>
> **Source boundary:** the subscription article text remains unavailable, but its original
> Supplementary Table 1 is present as an XLSX mirror. It directly confirms every
> CH-ATAC oligo used in the schematic. Detailed reaction conditions remain 🔴.

**CH-ATAC-seq** — Zhang G, Fu Y, Yang L, Ye F, Zhang P, Zhang S, Ma L, Li J, Wu H, Han X,
Wang J, Guo G. "Construction of single-cell cross-species chromatin accessibility
landscapes with combinatorial-hybridization-based ATAC-seq." *Developmental Cell*
59:793–811.e8 (2024). doi:[10.1016/j.devcel.2024.01.015](https://doi.org/10.1016/j.devcel.2024.01.015)
· PMID 38330939. Not open access (Europe PMC: "subscription required"; no PMC copy, no
preprint found).

Other papers:

- **CH-RNA-seq** (the RNA parent, same lab) — Ye F, Zhang G, E W, … Han X, Guo G.
  "Construction of the axolotl cell landscape using combinatorial hybridization sequencing
  at single-cell resolution." *Nat Commun* 13:4228 (2022).
  doi:[10.1038/s41467-022-31879-z](https://doi.org/10.1038/s41467-022-31879-z), PMC9307617
  (CC-BY 4.0). Not listed in the catalogue for this protocol; found by searching for the
  method name. Used here as the primary source for the shared oligos.
- Cited by upstream for the "semi-suppressive PCR" argument: nanoCAGE / CAGEscan,
  Plessy et al., *Nat Methods* 2010, doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470).
  Concept only; not fetched.

Sources read (in `_data/sources/ch-atac-seq__10.1016+j.devcel.2024.01.015/`, never committed):

| File | What | Used for |
|---|---|---|
| `scg_CH-ATAC-seq_SupplementaryTable1.xlsx.txt` | defining paper's original Supplementary Table 1, mirrored by scg_lib_structs | 🟢 common primers, 384 Tn5 barcode primers, 768 HY primers and 96 MGI P7 index primers |
| `upstream_CH-ATAC-seq.html.txt` | scg_lib_structs page for CH-ATAC-seq (secondary; says its oligos come from the paper's Supplementary Table 1) | **all CH-ATAC-specific oligos**, step order, final library, read layout |
| `related_CH-RNA-seq_PMC9307617.html.txt` | CH-RNA-seq full text (Europe PMC JATS XML, fetched by hand) | hybridization / blocking conditions, MGI sequencing with dark cycles |
| `related_CH-RNA-seq_MOESM4_ESM.xlsx.txt` | CH-RNA-seq Supplementary Data 1 = "Table S1. List of all oligonucleotide sequences used in CH-RNA-seq" (fetched by hand) | 🟢 HY head, 768 barcoded HY oligos, block tail, PCR P5 (with modification), 96 MGI P7 index primers |

Could not be fetched:

| What | URL | Why it matters |
|---|---|---|
| CH-ATAC-seq paper (STAR Methods) | https://doi.org/10.1016/j.devcel.2024.01.015 | every reaction condition; Tn5 loading; whether a ligase is used; cycle numbers; sequencing recipe |

---

## 1. What it is

A **three-level combinatorial-indexing scATAC-seq** for the BGI/MGI (DNBSEQ) platform, in
the line of sci-ATAC-seq3 🟡 (upstream):

1. **Barcoded Tn5** (one barcode per well, 384) tags nuclei in bulk-per-well.
2. Pool → split; a **barcoded hybridization (HY) oligo** (768) is annealed onto a handle
   carried by the Tn5 adapter — no barcode-specific ligation of a short splint as in
   SPLiT-seq; the HY oligo is itself a pre-annealed partial duplex.
3. Pool → split; **i7 index PCR** (96 MGI P7 index primers).

384 × 768 × 96 barcode combinations 🟡 (computed from the three primary-source sheets).

| | New thing here | Builds on |
|---|---|---|
| 1 | The second barcode is added by **hybridizing a pre-annealed duplex with a 23-nt single-stranded overhang** onto the Tn5-adapter's 5' handle | CH-RNA-seq does the same onto its RT primer; SPLiT-seq / sci-RNA-seq3 for split-pool |
| 2 | The Tn5 barcoded adapter carries **only ME + barcode + SMART handle** — no s5. s5 (the read-1 primer site) arrives with the HY oligo | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) with a non-standard "A" adapter |
| 3 | Fragments with the barcoded adapter at both ends get HY/s5 at both ends, and are **suppressed** in PCR (inverted repeats) rather than removed | semi-suppressive PCR (nanoCAGE, Plessy 2010) |
| 4 | MGI/DNBSEQ adapters, read 1 with **dark cycles** over the constant linker and the ME | CH-RNA-seq's MGI recipe |

## 2. Oligos

### 2a. CH-ATAC-specific oligos — 🟡 (upstream page only; paper's Table 1 not read)

As upstream writes them (`[…]` are upstream's barcode placeholders):

```
Tn5 ME bottom              /Phos/ CTGTCTCTTATACACATCT
Tn5_barcode_primer_{1..384}  AAGCAGTGGTATCAACGCAGAGT [10-bp Tn5 barcode] AGATGTGTATAAGAGACAG
Tn5_Primer_C_oligo         GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
MGI_P5                     GAACGACATGGCTACGATCCGACTT
MGI_P7                     TGTGAGCCAAGGAGTTGTTGTCTTC
Read 1 sequencing primer   TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Index i7 sequencing primer CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Read 2 sequencing primer   GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
```

The 384 Tn5 barcode sequences are not on the upstream page and were not obtained. 🔴

### 2b. Shared with CH-RNA-seq — 🟢 (CH-RNA-seq Supplementary Data 1)

Verbatim from `related_CH-RNA-seq_MOESM4_ESM.xlsx`, sheets "common primers",
"768 2nd hybridization primers", "96 MGI P7 indexed primers"; the modification column as
given:

```
HY_head_oligo            TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Block_tail_primer_oligo  AAGCAGTGGTATCAACGCAGAGT
PCR P5 primer            GAACGACATGGCTACGATCCGACTTTCGTCGGCAGCGTC     5'Phosphorylation
Barcoded_HY_oligo_1      ACTCTGCGTTGATACCACTGCTTTTCTCGCATGCTGTCTCTTATACACATCTGACGCTGCCGACGA
  … _1.._768             ACTCTGCGTTGATACCACTGCTT <10-nt HY bc> CTGTCTCTTATACACATCTGACGCTGCCGACGA   66 nt
MGI_P7_index_1           TGTGAGCCAAGGAGTTGTTGTCTTCTAGGTCCGATGTCTCGTGGGCTCGG
  … _1.._96              TGTGAGCCAAGGAGTTGTTGTCTTC <10-nt i7> GTCTCGTGGGCTCGG              50 nt
```

**Agreement with upstream:** HY_head, Block_tail, PCR P5 and the layout of the HY and
MGI P7 index oligos are identical to what upstream lists for CH-ATAC-seq (all 768 HY oligos
and all 96 P7 oligos have exactly the constant parts upstream gives — 🟡 computed). So
CH-ATAC-seq reuses the CH-RNA-seq second- and third-round oligo plates, or at least their
design. Whether the barcode lists are the same is 🔴 (CH-ATAC Table 1 not read).

**Disagreement with upstream:** the PCR P5 primer is **5'-phosphorylated** in CH-RNA-seq
🟢; upstream's CH-ATAC page shows no modification. A 5' phosphate fits the MGI workflow,
where the PCR product is circularized into single-stranded circles (§5). Whether the
CH-ATAC P5 is also phosphorylated is 🔴.

For comparison, the CH-RNA-seq first-round RT primer (not used in CH-ATAC) is 🟢
`A*AGCAGTGGTATCAACGCAGAGTNNNNNNNN<10-nt RT bc>TTTTTTTTTTTTTTTTTTTTTTTTTVN` (`*` =
phosphorothioate): the **same 23-nt handle** in front of the barcode. CH-ATAC-seq replaces
"handle + UMI + barcode + oligo-dT" with "handle + barcode + ME" on Tn5. 🟡

### How the oligos interlock — 🟡 (computed)

- **Tn5 ME bottom** = `nextera.ME_RC` (19 nt), 5'-phosphorylated: the standard
  non-transferred strand.
- **Tn5_barcode_primer** = `rt.SMART_HANDLE` (23) + 10-nt barcode + `nextera.ME` (19)
  = 52 nt. The 5' handle is exactly the SMART / ISPCR handle `AAGCAGTGGTATCAACGCAGAGT`.
  It has **no s5**.
- **Tn5_Primer_C_oligo** = `nextera.ADAPTOR_S7` (s7 + ME, 34 nt) — the ordinary Nextera
  "B" adapter.
- **HY_head_oligo** = `nextera.ADAPTOR_S5` (s5 + ME, 33 nt) — also identical to
  `nextera.READ1_PRIMER`, and to upstream's "Read 1 sequencing primer".
- **Barcoded_HY_oligo** (66 nt) = revcomp(`rt.SMART_HANDLE`) (23) + 10-nt barcode +
  revcomp(HY_head) (33; = `nextera.INDEX2_PRIMER`).
  - After pre-annealing, HY_head pairs with the 3' 33 nt of the HY oligo along its whole
    length → a blunt-ended duplex at the s5 end with a **23-nt 5' single-stranded
    overhang** = revcomp of the SMART handle, plus the 10-nt barcode left single-stranded
    between them.
  - That overhang is complementary to the 5' end of every Tn5_barcode_primer that was
    transferred into genomic DNA → this is the hybridization step.
- **Block_tail_primer_oligo** = `rt.SMART_HANDLE` = revcomp of the HY oligo's overhang.
  Added in excess after hybridization it occupies the overhang of free HY duplexes so they
  cannot capture handles after the next pooling.
- **PCR P5 primer** (39 nt) = MGI_P5 (25) + `nextera.S5` (14). Its 3' 14 nt are the 5' 14
  nt of HY_head, so it primes on the complement of s5 brought in by the HY oligo.
- **MGI_P7_index** (50 nt) = MGI_P7 (25) + 10-nt i7 + `nextera.S7` (15). Primes on the
  s7 brought in by Primer C.
- **Sequencing primers**: Read 1 = `nextera.READ1_PRIMER`, Index = `nextera.INDEX1_PRIMER`
  (ME_RC + s7_RC), Read 2 = `nextera.READ2_PRIMER` — i.e. all three are the standard
  Nextera ones, used on DNBSEQ (CH-RNA-seq calls them custom "TM (Tn5 modified)" primers 🟢).
- MGI_P5 and MGI_P7 are not in `lib/`; they are taken as given.

## 3. Step by step

Order and products from upstream 🟡; CH-ATAC conditions 🔴 (paper not read). Where the
shared hybridization step has published conditions in CH-RNA-seq, they are given and marked
as such.

1. **Assemble indexed Tn5.** Anneal ME bottom with Tn5_barcode_primer_N, and ME bottom
   with Primer C; load Tn5 with both. Each well's Tn5 carries one barcoded "A" adapter and
   the shared s7 "B" adapter. 🟡 (upstream). Loading ratio, enzyme source 🔴.
2. **Tagment nuclei** in 384 wells, one Tn5 barcode per well. Three products, each with
   the 9-nt gap ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)) 🟡:
   - barcode adapter at both ends → later gets HY/s5 at both ends, suppressed in PCR;
   - Primer C at both ends → no handle, never gets an HY oligo, cannot get P5 (only s7
     both ends);
   - one of each → the only productive fragment.
3. **Pool, split** into 768 wells (8 × 96 in CH-RNA-seq). 🟡
4. **Hybridize** the pre-annealed HY duplex (HY_head + Barcoded_HY_oligo_N) to the handle.
   CH-RNA-seq conditions 🟢 (not CH-ATAC): HY_head and barcoded HY oligo 50 µM each, mixed
   equally, 95 °C 2 min, ramp −0.1 °C/s to 25 °C; 2 µL of 25 µM duplex into 3 µL cells in
   hybridization buffer (50 mM Tris-HCl, 10 mM MgCl₂, 10 mM DTT, 0.1 % Triton X-100,
   10 % PEG8000, RNase inhibitor); 37 °C 90 min.
5. **Block**: Block_tail_primer_oligo (CH-RNA-seq: 0.5 µL of 100 µM, 37 °C 30 min 🟢), then
   pool.
6. **Join / fill.** Upstream says "ligate the HY barcode" and then "gap fill-in" 🟡. After
   hybridization the top strand has a **10-nt gap** (opposite the HY barcode) between the
   3' end of HY_head and the 5' end of the Tn5 barcode oligo; the bottom strand has a
   second **10-nt gap** opposite the Tn5 barcode (between the 3' end of ME bottom and the
   5' end of the HY oligo, as upstream's step 3.2 drawing also shows), besides the 9-nt Tn5
   gap on the insert side. 🟡 (computed from the oligo layout). CH-RNA-seq explicitly says its hybridization round works **without T4 ligase**
   🟢, and follows it with T4 PNK (37 °C 30 min, NEB) and a second-strand synthesis
   enzyme mix (16 °C 3 h; vendor not named in the text) 🟢. If that mix is the usual
   *E. coli* polymerase + ligase type, the gaps could be filled and sealed there 🟡
   (inference, not stated in the paper). Which enzymes CH-ATAC-seq uses, and whether the Tn5 barcode
   oligo carries the 5' phosphate a ligation would need, is 🔴.
7. **Pool, split** into 96 wells; lyse; **index PCR** with PCR P5 primer + one
   MGI_P7_index_N per well. A 72 °C step before denaturation presumably completes any remaining
   gap fill-in 🟡 (as in CH-RNA-seq's PCR
   program: 72 °C 5 min, then cycling 🟢 for CH-RNA-seq; 🔴 for CH-ATAC).
8. **MGI conversion** (CH-RNA-seq 🟢): circularize to ssDNA, make DNA nanoballs, sequence
   on DNBSEQ-T7. Presumably the same for CH-ATAC 🟡.

## 4. Final library — 🟡 (computed; equals upstream's drawn final library exactly)

Top strand, 5'→3':

```
5'- MGI_P5(25) · s5(14) · ME(19) · <HY bc'>(10) · SMART handle(23) · <Tn5 bc>(10) · ME(19)
    · <genomic insert> · ME'(19) · s7'(15) · <i7'>(10) · MGI_P7'(25) -3'
```

`'` = reverse complement of the sequence as written in the oligo. Assembling this from
`lib/` constants and the upstream oligos reproduces upstream's final-library top strand
character for character, and upstream's bottom strand is its exact complement 🟡
(computed). Fixed adapter length outside the insert: 25+14+19+10+23+10+19 = 120 nt on the
P5 side, 19+15+10+25 = 69 nt on the P7 side.

Barcode orientation 🟡 (computed):

- **HY barcode**: the top strand carries the **reverse complement** of the 10 nt written in
  Barcoded_HY_oligo (the HY oligo is the bottom strand). E.g. Barcoded_HY_oligo_1 has
  `TTCTCGCATG`; read 1 sees `CATGCGAGAA`.
- **Tn5 barcode**: as written in Tn5_barcode_primer (it is the top strand).
- **i7**: the index read (primer = ME_RC + s7_RC, which has top-strand sense and copies
  the bottom strand, so it reports top-strand sequence) reads the **reverse complement** of the 10 nt in MGI_P7_index. E.g. MGI_P7_index_1 has
  `TAGGTCCGAT`; the index read is `ATCGGACCTA`.

## 5. Sequencing

Upstream says it was guessing at MGI sequencing; its read layout 🟡:

| Read | Primer | Length | Contents |
|---|---|---|---|
| Read 1 | `nextera.READ1_PRIMER` (s5+ME) | 100 cycles, **dark 11–33 and 44–62** | 1–10 HY bc', 11–33 SMART handle (dark), 34–43 Tn5 bc, 44–62 ME (dark), 63–100 genomic |
| Index (i7) | `nextera.INDEX1_PRIMER` | 10 | i7' |
| Read 2 | `nextera.READ2_PRIMER` (s7+ME) | 100 | genomic |

The dark-cycle windows match the computed segment lengths exactly (10 + 23 = 33; 33 + 10 =
43; 43 + 19 = 62) 🟡 (computed). CH-RNA-seq used the same idea 🟢: DNBSEQ-T7, custom
"TM (Tn5 modified)" sequencing primers, read 1 = 51 cycles with dark 11–33, read 2 = 100
cycles. The CH-ATAC cycle numbers themselves are not confirmed 🔴.

## 6. Open questions

- 🔴 Every CH-ATAC reaction condition (Tn5 loading, tagmentation, nuclei numbers per well,
  hybridization, PCR cycles, size selection) — STAR Methods not read.
- 🔴 Is the HY barcode joined by ligation (upstream's word) or only by fill-in / extension
  (CH-RNA-seq states its hybridization round needs no T4 ligase)? Is the Tn5 barcode oligo 5'-phosphorylated?
- 🟢 The defining table supplies all 384 Tn5 barcodes, all 768 HY oligos and all 96 MGI
  P7 primers. Its PCR P5 entry explicitly carries a 5′ phosphate.
- 🟡 Product with the barcoded adapter at both ends: upstream calls it non-amplifiable by
  semi-suppressive PCR. With HY at both ends it carries s5 at both ends and only P5 primes
  it; whether it is truly negligible in the data is not stated.
- 🟡 CH-RNA-seq RT primer 146 is one base longer than the other 383 (69 vs 68 nt without
  `*`) — irrelevant to CH-ATAC, noted only because the table was checked.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` found no open copy of the Developmental Cell article text, but the
source catalogue retrieved the defining Supplementary Table 1 mirror and converted it with
`tools/doctext.py`. The related open-access CH-RNA-seq material independently agrees with
the shared hybridization system. `tools/scrape_primers.py` recognises `rt.SMART_HANDLE`,
`nextera.ME` and `ADAPTOR_S5/S7`; MGI_P5 / MGI_P7 remain protocol-local constants.
