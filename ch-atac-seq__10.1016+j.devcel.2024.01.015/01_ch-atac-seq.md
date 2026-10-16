# CH-ATAC-seq — combinatorial-hybridization single-cell ATAC-seq

> **Evidence marking.** 🟢 verbatim from a source read for this note · 🟡 derived or
> inferred · 🔴 not published / not available to us. Relationships marked 🟡 *(computed)*
> were worked out with `lib/`; the model is `tools/ch_atac.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length, step or read layout in this note or in the model is taken from it, and
> agreement with it is not evidence. An earlier revision of this note was built around
> agreement with that page (§2 "Agreement/Disagreement with upstream", a §4 that said the
> final library "equals upstream's drawn final library exactly", a §3 whose step order was
> "from upstream", and a §5 built on "its read layout"). All of that is withdrawn.
>
> **Source boundary.** The defining paper's own **Supplementary Table 1** is on disk and is
> the authority for every CH-ATAC base drawn: generic Illumina/Nextera bases come from
> `lib/nextera.py` and the genomic insert is a placeholder. The subscription STAR Methods
> are **not** available, so no reaction condition, enzyme choice, well count per round or
> sequencing cycle number is asserted as CH-ATAC's own 🔴; where the open-access sibling
> **CH-RNA-seq** publishes the same shared chemistry, it is cited as that, and the transfer
> to CH-ATAC is marked 🟡.
>
> **On the `scg_` filename.** `scg_CH-ATAC-seq_SupplementaryTable1.xlsx` is a *mirror of the
> publisher's supplement*, not scg-authored content: its first row reads "related to figure
> 1, oligo primers used in this research", it has Cell Press's `related to figure`
> convention and an IDT-style `Modification` column, and its four sheets ("common primers",
> "384 1st Tn5 Barcode primers", "768 2nd hybridization primers", "96 MGI P7 indexed
> primers") are the ordered oligo inventory. A filename prefix is a download route, not
> provenance; this is an ordinary original and is cited as the paper's Supplementary Table 1
> throughout.

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
- `catalogue/scg_lib_structs.tsv` carries a `cited` row for nanoCAGE / CAGEscan, Plessy
  et al., *Nat Methods* 2010, doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470),
  scraped from the drawn page's own citation list. Background reading on suppression PCR;
  not fetched, and nothing here rests on it 🟡.

Sources read (in `_data/sources/ch-atac-seq__10.1016+j.devcel.2024.01.015/`, never committed):

| File | What | Used for |
|---|---|---|
| `scg_CH-ATAC-seq_SupplementaryTable1.xlsx.txt` | **the defining paper's own Supplementary Table 1** — a publisher original that happens to have been downloaded through the scg_lib_structs `data/` mirror (see the note above) | 🟢 **every base in the model**: common primers, 384 Tn5 barcode primers, 768 HY primers and 96 MGI P7 index primers |
| `upstream_CH-ATAC-seq.html.txt` | scg_lib_structs page for CH-ATAC-seq | **inspiration only, not a source** — nothing in this note or the model is cited from, or verified against, it |
| `related_CH-RNA-seq_PMC9307617.html.txt` | CH-RNA-seq full text (Europe PMC JATS XML, fetched by hand) | hybridization / blocking conditions, MGI sequencing with dark cycles |
| `related_CH-RNA-seq_MOESM4_ESM.xlsx.txt` | CH-RNA-seq Supplementary Data 1 = "Table S1. List of all oligonucleotide sequences used in CH-RNA-seq" (fetched by hand) | 🟢 HY head, 768 barcoded HY oligos, block tail, PCR P5 (with modification), 96 MGI P7 index primers |

Could not be fetched:

| What | URL | Why it matters |
|---|---|---|
| CH-ATAC-seq paper (STAR Methods) | https://doi.org/10.1016/j.devcel.2024.01.015 | every reaction condition; Tn5 loading; whether a ligase is used; cycle numbers; sequencing recipe |
| Research Square preprint / Protocol Exchange record | https://doi.org/10.21203/rs.3.rs-960203/v1 | `MANIFEST.tsv` already carries this as a `(manual)` row (Research Square sits behind Cloudflare). It is the most likely open route to the same STAR Methods, and would settle most of §6 — worth a hand download. |

---

## 1. What it is

A **three-level combinatorial-indexing scATAC-seq** for the BGI/MGI (DNBSEQ) platform. The
three rounds, and their order, come from Supplementary Table 1's own sheet names and
Description column 🟢 — "384 **1st** Tn5 Barcode primers", "768 **2nd** hybridization
primers", "96 MGI P7 indexed primers", with "Anneal with 2nd hybridization primers" against
`HY_head_oligo`, "Block free annealing hybridization primers" against
`Block_tail_primer_oligo`, "Anneal with Tn5 ME oligo as Transposase P5/P7 Adaptor" against
the two Tn5 arms, and "Library preparation PCR P5 primer" against the P5 oligo:

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

### 2a. CH-ATAC-specific oligos — 🟢 (defining Supplementary Table 1)

Verbatim from `scg_CH-ATAC-seq_SupplementaryTable1.xlsx.txt`, a mirror of the
defining paper's original Supplementary Table 1. The representative sequences are
the first row of each indexed family; `*` is a phosphorothioate bond, as the table
states. The table gives all 384 Tn5, 768 HY and 96 MGI P7 index oligos.

```
Tn5 ME oligo              /5Phos/CTGTCTCTTATACACATCT
Tn5_barcode_primer_1     A*AGCAGTGGTATCAACGCAGAGTTCCTACCAGTAGATGTGTATAAGAGACAG
Tn5 primer C oligo       GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
HY_head_oligo            TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Block_tail_primer_oligo  AAGCAGTGGTATCAACGCAGAGT
PCR P5 primer            /5Phos/GAACGACATGGCTACGATCCGACTTTCGTCGGCAGCGTC
Barcoded_HY_oligo_1      ACTCTGCGTTGATACCACTGCTTTTCTCGCATGCTGTCTCTTATACACATCTGACGCTGCCGACGA
MGI_P7_index_1          TGTGAGCCAAGGAGTTGTTGTCTTCTAGGTCCGATGTCTCGTGGGCTCGG
```

**Reading `A*AGCAGTGG…`.** The table's Modification column for that sheet reads
`*:Phosphorthioate`, so the `*` marks a phosphorothioate **bond** between the first and
second base — the same convention the table uses for the HY and P5 entries. The oligo is
therefore **52 bases** beginning `AAGCAGTGG…`, which is exactly `rt.SMART_HANDLE` (23 nt)
🟢. Two checks inside the same table confirm it: `Block_tail_primer_oligo` is printed as the
bare 23-nt handle `AAGCAGTGGTATCAACGCAGAGT` 🟢, and every `Barcoded_HY_oligo` opens with
its 23-nt reverse complement `ACTCTGCGTTGATACCACTGCTT` 🟢. There is **no extra 5′ A**: an
earlier revision of this note and of `tools/ch_atac.py` read `A*` as a capping base, drew a
53-nt oligo, and consequently marked the HY junction inferred. That is corrected — see the
junction paragraph in §4. 5′ phosphates and phosphorothioate bonds are recorded here but
are not drawn as bases.

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

**Two independent first-party tables agree.** CH-ATAC Supplementary Table 1 and CH-RNA-seq
Supplementary Data 1 list the *same* `HY_head_oligo`, `Block_tail_primer_oligo`,
`PCR P5 primer`, the same 768 `Barcoded_HY_oligo` sequences and the same 96
`MGI_P7_index` primers, base for base, including the 5′-phosphorylation on the PCR P5
primer 🟢. That is a genuine two-source confirmation of the shared half of the chemistry —
two supplements from the same lab, published in two journals, not a drawn page and a
supplement. What is CH-ATAC-only is the 384 `Tn5_barcode_primer` family and the two Tn5
arms; those rest on CH-ATAC Supplementary Table 1 alone 🟢.

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
  `nextera.READ1_PRIMER`, which is why the standard Nextera Read 1 primer lands on this
  library at all 🟡 *(computed; `seqprimers.verify` enforces it)*.
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

**Step order** from Supplementary Table 1's sheet names and Description column (§1) 🟢;
**products** computed from the oligos with `lib/` 🟡; **CH-ATAC conditions** 🔴 (STAR
Methods paywalled). Where the shared hybridization step has published conditions in
CH-RNA-seq, they are quoted and marked as CH-RNA-seq's.

1. **Assemble indexed Tn5.** Anneal ME bottom with Tn5_barcode_primer_N, and ME bottom
   with Primer C; load Tn5 with both. Each well's Tn5 carries one barcoded "A" adapter and
   the shared s7 "B" adapter. 🟢 (Supplementary Table 1 labels the two arms "Transposase
   P5 Adaptor" and "Transposase P7 Adaptor", both annealed to the same `Tn5 ME oligo`).
   Loading ratio, enzyme source 🔴.
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
6. **Join / fill.** 🔴 CH-ATAC does not say which enzymes do this. After
   hybridization the top strand has a **10-nt gap** (opposite the HY barcode) between the
   3' end of HY_head and the 5' end of the Tn5 barcode oligo; the bottom strand has a
   second **10-nt gap** opposite the Tn5 barcode (between the 3' end of ME bottom and the
   5' end of the HY oligo), besides the 9-nt Tn5
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

## 4. Final library — 🟡 (computed from Supplementary Table 1 with `lib/`)

Top strand, 5'→3':

```
5'- MGI_P5(25) · s5(14) · ME(19) · <HY bc'>(10) · SMART handle(23) · <Tn5 bc>(10) · ME(19)
    · <genomic insert> · ME'(19) · s7'(15) · <i7'>(10) · MGI_P7'(25) -3'
```

`'` = reverse complement of the sequence as written in the oligo. Fixed adapter length
outside the insert: 25+14+19+10+23+10+19 = 120 nt on the P5 side, 19+15+10+25 = 69 nt on
the P7 side 🟡 *(computed)*.

**The HY/Tn5 junction is determined, not inferred.** `revcomp(Barcoded_HY_oligo)` reads
`s5 · ME · HY-barcode' · SMART handle` and the indexed Tn5 transfer strand reads
`SMART handle · Tn5 barcode · ME`. `ch_atac.hy_tn5_junction()` searches every possible
overlap between the two and finds **exactly one**, 23 nt long and equal to the SMART
handle; it raises if that ever stops being true. So the order
`… ME · HY-bc' · handle · Tn5-bc · ME …` follows from the two printed oligos alone, with
no appeal to the unavailable methods text and none to a drawing 🟡 *(computed, enforced by
construction)*. The earlier "inferred junction" marking on this segment, and the page
caveat that went with it, existed only because the model had mis-read `A*AGCAGTGG…` as
carrying an extra 5′ A (§2a) and so saw a 24-nt rather than a 23-nt boundary.

Barcode orientation 🟡 (computed):

- **HY barcode**: the top strand carries the **reverse complement** of the 10 nt written in
  Barcoded_HY_oligo (the HY oligo is the bottom strand). E.g. Barcoded_HY_oligo_1 has
  `TTCTCGCATG`; read 1 sees `CATGCGAGAA`.
- **Tn5 barcode**: as written in Tn5_barcode_primer (it is the top strand).
- **i7**: the index read (primer = ME_RC + s7_RC, which has top-strand sense and copies
  the bottom strand, so it reports top-strand sequence) reads the **reverse complement** of the 10 nt in MGI_P7_index. E.g. MGI_P7_index_1 has
  `TAGGTCCGAT`; the index read is `ATCGGACCTA`.

## 5. Sequencing

**The CH-ATAC cycle counts remain unverified** 🔴 — the defining STAR Methods and
linked preprint have not been read. What *is* established: all three primers are the standard Nextera ones,
located on the modelled library by `seqprimers.verify` (the build fails if any fails to
land) 🟡 *(computed)*; the index read is 10 cycles because the i7 is 10 nt 🟡 *(computed)*;
and Read 1 must exceed 62 cycles to reach genomic DNA 🟡 *(computed)*. The paired-100
figures below are the model's assumption, declared as such in the graph's
`SequencingPlan(limitation=…)` and in the page's sequencing preamble, and supported only by
the sibling CH-RNA-seq paper's published recipe (below).

| Read | Primer | Length (assumed) | Contents |
|---|---|---|---|
| Read 1 | `nextera.READ1_PRIMER` (s5+ME) | 100 cycles, **dark 11–33 and 44–62** | 1–10 HY bc', 11–33 SMART handle (dark), 34–43 Tn5 bc, 44–62 ME (dark), 63–100 genomic |
| Index (i7) | `nextera.INDEX1_PRIMER` | 10 | i7' |
| Read 2 | `nextera.READ2_PRIMER` (s7+ME) | 100 | genomic |

The dark-cycle windows are **computed** from the segment lengths, not transcribed:
10 + 23 = 33; 33 + 10 = 43; 43 + 19 = 62 🟡. The first window is independently published
for the sibling protocol — CH-RNA-seq states verbatim 🟢
(`related_CH-RNA-seq_PMC9307617.html.txt:53`):

> "Amplified DNA nanoballs (DNBs) were sequenced with custom TM (Tn5 modified) sequencing
> primers on MGI DNBSEQ-T7 with a dark reaction model (CH-RNA-seq: 51 cycles of read1 with
> a dark reaction from 11 to 33, 100 cycles of read2)."

CH-RNA-seq's `11 to 33` is the same 10-nt HY barcode followed by the same 23-nt shared
handle, so that window transfers to CH-ATAC by construction 🟡. The **second** window
(44–62, over the Tn5 barcode and the mosaic end) is CH-ATAC-only and is computed, not
verified for CH-ATAC 🔴 — CH-RNA-seq has no Tn5 barcode at that position. The CH-ATAC cycle numbers
remain unresolved until its original Methods are checked 🔴.

## 6. Open questions

- 🔴 Every CH-ATAC reaction condition (Tn5 loading, tagmentation, nuclei numbers per well,
  hybridization, PCR cycles, size selection) — STAR Methods not read.
- 🔴 Whether the HY barcode is joined by ligation or only by fill-in / extension. The
  model says nothing either way: `hy_fill_{a,b}` is a plain `StepKind.EXTENSION` with an
  unnamed polymerase and no ligase input. CH-RNA-seq states its hybridization round needs
  no T4 ligase 🟢, which is suggestive but is a different protocol 🟡. Supplementary Table 1
  marks the `Tn5 ME oligo` and the `PCR P5 primer` 5′-phosphorylated but gives no
  modification for the Tn5 barcode primer 🟢, so if a ligation happened there it would need
  a phosphate the table does not list 🟡.
- 🟢 The defining table supplies all 384 Tn5 barcodes, all 768 HY oligos and all 96 MGI
  P7 primers. Its PCR P5 entry explicitly carries a 5′ phosphate.
- 🟡 Product with the barcoded adapter at both ends: it carries s5 at both ends, so only
  the P5 primer can prime it and the resulting inverted repeat is suppressed in PCR
  *(computed from the two arms)*; the paper does not quantify how negligible it is 🔴, and
  the model therefore routes it to `graph.discard(tagment_waste_*)` with that reason rather
  than claiming a mechanism.
- 🟡 CH-RNA-seq RT primer 146 is one base longer than the other 383 (69 vs 68 nt without
  `*`) — irrelevant to CH-ATAC, noted only because the table was checked.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` found no open copy of the Developmental Cell article text, but the
source catalogue retrieved the defining Supplementary Table 1 mirror and converted it with
`tools/doctext.py`. The related open-access CH-RNA-seq material independently agrees with
the shared hybridization system. `tools/scrape_primers.py` recognises `rt.SMART_HANDLE`,
`nextera.ME` and `ADAPTOR_S5/S7`; MGI_P5 / MGI_P7 remain protocol-local constants.

What the tooling did **not** catch was the `A*` mis-read (§2a): the indexed Tn5 transposome
is an `ADAPTER` input built from two annealed strands, so
`tools/selftest_source_oligos.py` has no single whole-strand sequence to compare and
reported `inputs without exact states` rather than a failure. The error surfaced only on
re-reading the table's Modification column and counting the oligo against
`Block_tail_primer_oligo`. The two MGI PCR primers now carry exact graph states, so the
guard does check those two against the table text.
