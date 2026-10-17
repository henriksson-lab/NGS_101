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
> `lib/nextera.py` and the genomic insert is a placeholder. The defining article's
> **STAR Methods** were obtained on 2026-10-17 (publisher PDF with supplemental
> information, `mmc3`) and are the authority for reaction order, enzymes, conditions and
> the sequencing recipe 🟢. The sibling **CH-RNA-seq** paper is now background only; where
> it is still quoted, it is cited as CH-RNA-seq's.
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
· PMID 38330939. No PMC copy (Europe PMC: "subscription required"). Crossref
(checked 2026-10-17) lists an Elsevier open-archive VOR licence from 2025-03-25
(PII S1534580724000352). The article PDF with supplemental information (`mmc3`) and the
supplements `mmc1` (figures) and `mmc2` (tables) were hand-fetched on 2026-10-17 from the
Elsevier publisher CDN (`ars.els-cdn.com`) and read.

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

Sources read (in `/husky/henriksson/NGS101_papers/sources/ch-atac-seq__10.1016+j.devcel.2024.01.015/`, never committed):

| File | What | Used for |
|---|---|---|
| `scg_CH-ATAC-seq_SupplementaryTable1.xlsx.txt` | **the defining paper's own Supplementary Table 1** — a publisher original that happens to have been downloaded through the scg_lib_structs `data/` mirror (see the note above) | 🟢 **every base in the model**: common primers, 384 Tn5 barcode primers, 768 HY primers and 96 MGI P7 index primers |
| `j.devcel.2024.01.015_1-s2.0-S1534580724000352-mmc3.pdf.txt` (and the PDF, p. e4 rendered) | **the defining article** with STAR Methods; the text twin renders µ as "m" (checked against the rendered page) | 🟢 Results lines 172-192; Key resources lines 1470-1500, 1645-1692; "Library construction of CH-ATAC-seq" lines 1834-1863 (p. e4); "MGI library preparation and sequencing" lines 1884-1888 (p. e5) |
| `j.devcel.2024.01.015_1-s2.0-S1534580724000352-mmc2.zip` | published Supplementary Tables 1, 3, 5, 6, 7 | Table S1 has the same md5 (`712d8483…`) as the `scg_` mirror above; Tables 3/5/6/7 are cell-lineage, annotation, TF and regulon tables (sheet names checked), no chemistry |
| `j.devcel.2024.01.015_1-s2.0-S1534580724000352-mmc1.pdf.txt` | supplemental figures SFig. 1-7 and Supplementary Tables 2 and 4 | legends checked: all relate to Figures 2-7 (atlases, cCREs, TFs, GRNs), none to library design |
| `upstream_CH-ATAC-seq.html.txt` | scg_lib_structs page for CH-ATAC-seq | **inspiration only, not a source** — nothing in this note or the model is cited from, or verified against, it |
| `related_CH-RNA-seq_PMC9307617.html.txt` | CH-RNA-seq full text (Europe PMC JATS XML, fetched by hand) | hybridization / blocking conditions, MGI sequencing with dark cycles |
| `related_CH-RNA-seq_MOESM4_ESM.xlsx.txt` | CH-RNA-seq Supplementary Data 1 = "Table S1. List of all oligonucleotide sequences used in CH-RNA-seq" (fetched by hand) | 🟢 HY head, 768 barcoded HY oligos, block tail, PCR P5 (with modification), 96 MGI P7 index primers |

Could not be fetched:

| What | URL | Why it matters |
|---|---|---|
| Research Square record doi:10.21203/rs.3.rs-960203/v1 | https://doi.org/10.21203/rs.3.rs-960203/v1 | **Not a CH-ATAC source.** Crossref metadata (checked 2026-10-17) titles it "The single-cell transcriptional landscape of the axolotl" (posted 2021-10-13) with `is-preprint-of` 10.1038/s41467-022-31879-z, i.e. the CH-RNA-seq preprint, superseded by PMC9307617. It cannot supply the CH-ATAC STAR Methods. |

History: before the hand download, these routes were tried for the Dev Cell article on
2026-10-17 (round deferred-followup-02), all without full text: OpenAlex/Unpaywall (bronze, publisher location only, no repository
copy); Europe PMC (MED 38330939, "Subscription required", not in PMC); cell.com fulltext
and PDF URLs and ScienceDirect `pdfft` (Cloudflare 403); Elsevier article API
(metadata confirms `openArchiveArticle=true`, full text needs an API key); no CH-ATAC
preprint found (Europe PMC and Crossref title searches). GEO SubSeries GSE232043 (zebrafish),
GSE232271 (Drosophila), GSE232272 (earthworm) under SuperSeries GSE232162 were found;
GSM7324821's record (`geo_GSM7324821_self_brief.soft.txt` lines 19-22, 33) gives nuclei
lysis buffer, "Library preparations were performed with the CH ATAC-seq protocol",
two FASTQ files per library with the HY (768) and Tn5 (384) indices "in specific base
positions" of one of them, and DNBSEQ-T7 🟢 — no reaction conditions or read lengths.
The PDF was then obtained from the publisher CDN (above).

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
3. Pool → T4 PNK → split; gap fill and seal; lysis; **i7 index PCR** (96 MGI P7 index
   primers).

384 × 768 × 96 barcode combinations 🟡 (computed from the three primary-source sheets).

**Published discrepancy — number of PCR indices.** The Results (`mmc3.pdf.txt:183-185`,
p. 795) say "a third round of indices were introduced by 384 indexed P7 adapters" 🟢. Supplementary
Table 1 has a sheet named "96 MGI P7 indexed primers" with 96 sequences
(`scg_CH-ATAC-seq_SupplementaryTable1.xlsx.txt:1165-1262`) 🟢, the Key resources list
"3rd MGI P7 indexed primers (Table S1)", and the STAR Methods split nuclei into "one or
more 96-well plates" 🟢. The model draws the 96 printed primers. How 384 PCR indices would
be made from them is not stated 🔴. The page shows this discrepancy.

| | New thing here | Builds on |
|---|---|---|
| 1 | The second barcode is added by **hybridizing a pre-annealed duplex with a 23-nt single-stranded overhang** onto the Tn5-adapter's 5' handle | CH-RNA-seq does the same onto its RT primer; SPLiT-seq / sci-RNA-seq3 for split-pool |
| 2 | The Tn5 barcoded adapter carries **only ME + barcode + SMART handle** — no s5. s5 (the read-1 primer site) arrives with the HY oligo | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) with a non-standard "A" adapter |
| 3 | Fragments with the barcoded adapter at both ends get HY/s5 at both ends, and are **suppressed** in PCR (inverted repeats) rather than removed | semi-suppressive PCR (nanoCAGE, Plessy 2010) |
| 4 | MGI/DNBSEQ adapters, read 1 with **dark cycles** over the constant linker and the ME | CH-RNA-seq's MGI recipe (one dark window); CH-ATAC prints two |

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

**Step order and CH-ATAC conditions** from the defining STAR Methods, "Library
construction of CH-ATAC-seq" (`mmc3.pdf.txt:1834-1863`, p. e4; volumes in µL/µM as on
the rendered page) 🟢; **products** computed from the oligos with `lib/` 🟡.

1. **Assemble indexed Tn5.** Per 96-well plate, Tn5 primer C oligo (25 pM), Tn5 ME oligo
   (50 pM) and barcoded Tn5 primer oligos (25 pM) are mixed, heated to 95 °C 2 min and
   cooled to 10 °C at 0.1 °C/s; Tn5 (Vazyme S111, TruePrep Tagment Enzyme) is mixed with
   the annealed primers, 30 °C 60 min, final indexed Tn5 40 ng/µL 🟢. (The "pM" units are
   as printed.) Each well's Tn5 carries one barcoded "A" adapter and the shared s7 "B"
   adapter; Supplementary Table 1 labels the two arms "Transposase P5 Adaptor" and
   "Transposase P7 Adaptor", both annealed to the same `Tn5 ME oligo` 🟢.
2. **Tagment nuclei** in four 96-well plates (384 wells): 23 µL nuclei (10,000-15,000) +
   2 µL indexed Tn5, 55 °C 30 min, then ice 5 min 🟢. Three products, each with the 9-nt
   gap ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)) 🟡:
   - barcode adapter at both ends → later gets HY/s5 at both ends, suppressed in PCR;
   - Primer C at both ends → no handle, never gets an HY oligo, cannot get P5 (only s7
     both ends);
   - one of each → the only productive fragment.
3. **Pool, wash, split** into eight 96-well plates (768 wells, 8 µL per well) in
   hybridization buffer (50 mM Tris-HCl, 10 mM MgCl₂, 10 mM DTT, 0.1 % Triton X-100,
   10 % PEG8000) 🟢.
4. **Hybridize** 2 µL of 25 µM pre-annealed hybridization primers (HY head oligo and
   barcoded HY primer oligo, 50 µM each, mixed equally, 95 °C 2 min, cooled to 25 °C at
   0.1 °C/s) per well, 37 °C 90 min 🟢.
5. **Block**: 0.5 µL of 100 µM block tail primer oligo per well "to block any redundant
   hybridization primers", 30 min at 37 °C, then pool 🟢. Supplementary Table 1 describes
   `Block_tail_primer_oligo` as "Block free annealing hybridization primers" 🟢; the model
   adds it in each HY well before pooling (graph steps `block_a/b`, exact oligo state
   checked against the table).
6. **Phosphorylate in the pool.** Wash, resuspend in 40 µL; add 60 µL PNK mix (10 µL 10X
   PNK buffer, 20 µL T4 polynucleotide kinase, 10 µL 10 mM ATP, NEB), 37 °C 30 min; stop
   with 1 mL ice-cold wash buffer, filter, pellet 🟢 (graph step `pnk`).
7. **Re-split and gap-fill.** Nuclei are split equally into one or more 96-well plates
   (8 µL, 3,000-10,000 nuclei per well), and 2 µL gap-fill mix is added per well: 1 µL 10X
   Ampligase DNA ligase buffer, 0.3 µL Ampligase DNA ligase, 0.5 µL 10 mM dNTP, 0.2 µL T4
   DNA polymerase; 37 °C 60 min 🟢 (graph steps `hy_fill_a/b`; placement printed, so no
   longer inferred). The Methods text names the ligase vendor "Lucigen"; the Key resources
   list "Ampligase Thermostable DNA Ligase Kit, Biosearch, Cat# A8101" and "T4 DNA
   Polymerase, NEB, Cat# M0203S" 🟢. After hybridization the top strand has a **10-nt
   gap** (opposite the HY barcode) between the 3' end of HY_head and the 5' end of the Tn5
   barcode oligo; the bottom strand has a second **10-nt gap** opposite the Tn5 barcode
   (between the 3' end of ME bottom and the 5' end of the HY oligo), besides the 9-nt Tn5
   gap on the insert side 🟡 *(computed from the oligo layout)*. That T4 DNA polymerase
   fills these gaps and Ampligase seals the resulting nicks (after PNK has supplied 5'
   phosphates) is the natural reading of the printed mix 🟡; the paper does not draw the
   junctions, so the fill products carry no exact state in the graph.
8. **Lyse and index PCR** in the same well: 2 µL nuclei lysis buffer (Tris-HCl,
   proteinase K, SDS), 55 °C 15 min; 12 µL 10 % Tween quench; 28 µL PCR mix (26 µL 2X
   KAPA HiFi HotStart ReadyMix, 1 µL 10 µM P5 primer, 1 µL 10 µM indexed P7 primer).
   Program: 72 °C 5 min; 98 °C 3 min; **12 cycles** of 98 °C 10 s, 65 °C 15 s, 72 °C 1 min;
   72 °C 5 min 🟢. Pool, purify with 0.5-1.0X VAHTS DNA Clean Beads 🟢. No other size
   selection is printed.
   The quench is printed as "12 µL 10% Tween-10" 🟢 (also "Tween-10" in the bulk
   benchmark, `mmc3.pdf.txt:1859, 1880`); the detergent identity is not modelled.

**Reduced-index variant (not drawn).** For the fresh and formalin-fixed HEK 293T/NIH 3T3
mixtures, "we only introduced 96 indices in tagmentation, 384 indices in hybridization and
32 indices in PCR reaction here" 🟢 (`mmc3.pdf.txt:1883`, p. e5). The bulk ATAC
benchmark tagments three wells with indexed Tn5 and goes straight to lysis and indexed PCR
(no HY step) 🟢 (`mmc3.pdf.txt:1878-1882`). The model draws only the full CH-ATAC design.

**Published discrepancy — nuclei per final well.** The Results (`mmc3.pdf.txt:182`,
p. 795) collect "3,000–6,000 nuclei for each well of 96-well plates" after PNK; the STAR
Methods (`mmc3.pdf.txt:1856`, p. e4) load "3,000-10,000 nuclei" per well 🟢. Loading
is not part of the molecular model, so the drawn chemistry is unaffected; the page states
the discrepancy and which range applies to which experiment remains unstated 🔴.
9. **MGI conversion** ("MGI library preparation and sequencing", `mmc3.pdf.txt:1884-1888`,
   p. e5): circularize to ssDNA with a VAHTS Circularization Kit for MGI (Vazyme), make
   DNBs with a DNBSEQ DNB preparation kit (MGI), sequence on DNBSEQ-T7 🟢. The model ends
   at the linear PCR library and states the conversion in its sequencing limitation; the
   circularization chemistry is not drawn here.

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

The STAR Methods print the recipe verbatim 🟢 (`mmc3.pdf.txt:1886-1888`, p. e5):

> "Amplified DNA nanoballs (DNBs) were sequenced with custom TM (Tn5 modified) sequencing
> primers on MGI DNBSEQ-T7 with a dark reaction model (100 cycles of read1 with dark
> reaction from 11-33, 44-62, 100 cycles of read2)."

The graph's `SequencingPlan` declares Read 1 = 100 and Read 2 = 100 with that source
statement in its limitation. The custom TM primer bases are **not printed** 🔴, so the
model places the standard Nextera primers by computation (`seqprimers.verify`; the build
fails if any fails to land) 🟡. The index read is 10 cycles because the i7 is 10 nt 🟡
*(computed)*; the paper prints no index-read length 🔴.

| Read | Primer (modelled) | Length | Contents |
|---|---|---|---|
| Read 1 | `nextera.READ1_PRIMER` (s5+ME) | 100 cycles 🟢, **dark 11–33 and 44–62** 🟢 | 1–10 HY bc', 11–33 SMART handle (dark), 34–43 Tn5 bc, 44–62 ME (dark), 63–100 genomic |
| Index (i7) | `nextera.INDEX1_PRIMER` | 10 🟡 | i7' |
| Read 2 | `nextera.READ2_PRIMER` (s7+ME) | 100 🟢 | genomic |

The modelled read-1 segments give the same windows as the printed ones: 10 + 23 = 33;
33 + 10 = 43; 43 + 19 = 62 🟡 *(computed)*. That agreement supports a read-1 primer ending
in s5+ME, which is what the model places. The same first window (`11 to 33`) is printed
for CH-RNA-seq, which shares the HY barcode and handle (`related_CH-RNA-seq_PMC9307617.html.txt:53`).
GEO GSM7324821 says one of two FASTQ files carries both the HY (768) and Tn5 (384) indices,
and the STAR Methods "Preprocessing" paragraph says the same 🟢.

## 6. Open questions

- 🟢 Resolved 2026-10-17 from the STAR Methods: Tn5 loading, tagmentation, well counts,
  hybridization, blocking, PNK, gap fill (T4 DNA polymerase + Ampligase in the final
  wells), lysis, the 12-cycle KAPA HiFi PCR, bead cleanup and the DNBSEQ-T7 read recipe
  (§3, §5). No size selection other than the 0.5-1.0X bead cleanup is printed.
- 🟡 The ligation is printed (Ampligase in the gap-fill mix, after T4 PNK); which nicks
  it seals is computed from the oligo layout. Supplementary Table 1 gives no
  5′-phosphate for the Tn5 barcode primer or the barcoded HY oligo 🟢; the pooled T4 PNK
  step would supply one 🟡.
- 🔴 Bases of the "custom TM (Tn5 modified) sequencing primers" and the index-read length.
- 🔴 How the Results' "384 indexed P7 adapters" relate to the 96 printed MGI P7 primers
  (§1).
- 🟢 The defining table supplies all 384 Tn5 barcodes, all 768 HY oligos and all 96 MGI
  P7 primers. Its PCR P5 entry explicitly carries a 5′ phosphate.
- 🟡 Product with the barcoded adapter at both ends: it carries s5 at both ends, so only
  the P5 primer can prime it and the resulting inverted repeat is suppressed in PCR
  *(computed from the two arms)*; the paper does not quantify how negligible it is 🔴. The
  model routes same-ended products to `graph.discard(tagment_waste_*)` with a reason stating
  that no physical separation is represented and the fate is computed, not observed.
- 🟢 Resolved 2026-10-17: the Research Square record doi:10.21203/rs.3.rs-960203/v1 in
  `MANIFEST.tsv` is the CH-RNA-seq (axolotl) preprint per its Crossref `is-preprint-of`
  relation, not a CH-ATAC record; it is superseded by the published CH-RNA-seq article and
  gives no route to the CH-ATAC STAR Methods.
- 🟢 GEO GSM7324821 confirms DNBSEQ-T7 and that one of the two FASTQ files carries both
  the HY (768) and Tn5 (384) indices; read lengths are not stated there 🔴.
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
