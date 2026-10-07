# scTHS-seq — single-nucleus chromatin accessibility by T7 transposons and in vitro transcription

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`).

**scTHS-seq** — Lake BB, Chen S, Sos BC, Fan J, Kaeser GE, Yung YC, Duong TE, Gao D,
Chun J, Kharchenko PV, Zhang K. "Integrative single-cell analysis of transcriptional and
epigenetic states in the human adult brain." *Nature Biotechnology* 36:70–80 (2018).
doi:[10.1038/nbt.4038](https://doi.org/10.1038/nbt.4038) · PMC5951394 · PMID 29227469.
Data: GEO GSE97942. (The same paper introduces snDrop-seq; this note covers only scTHS-seq.)

Builds on **THS-seq** (bulk) — Sos BC, Fung H-L, Gao DR, Osothprarop TF, Kia A, He MM,
Zhang K. "Characterization of chromatin accessibility with a transposome hypersensitive
sites sequencing (THS-seq) assay." *Genome Biology* 17:20 (2016).
doi:[10.1186/s13059-016-0882-7](https://doi.org/10.1186/s13059-016-0882-7) · PMC4743176.
The engineered enzyme Tn5059 is from Kia et al., *BMC Biotechnol.* 17:6 (2017) (not fetched).

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/scths-seq__10.1038+nbt.4038/`, never committed):

| File | What | Used for |
|---|---|---|
| `nbt.4038_PMC5951394.html.txt` | full text (PMC author manuscript) | **all scTHS-seq methods**, sequencing |
| `nbt.4038_41587_2018_BFnbt4038_MOESM11_ESM.xlsx.txt` | Supplementary Tables 1–13; **Table S12** | **every oligo** (430 rows: 384 r5 transposons, 16 i5, 24 i7, adaptors, read primers) |
| `nbt.4038_41587_2018_BFnbt4038_MOESM9_ESM.pdf.txt` | Supplementary Figures + legends | Fig. S2 (transposon design, read layout in words) |
| `nbt.4038_41587_2018_BFnbt4038_MOESM10_ESM.pdf.txt` | Reporting summary | nothing chemical |
| `s13059-016-0882-7_PMC4743176.html.txt` | THS-seq full text | bulk ancestor: rationale, IVT and RNA-seq conditions |
| `s13059-016-0882-7_supp_13059_2016_882_MOESM1_ESM.pdf.txt` | THS-seq Additional file 1 | bulk THS-seq oligos (Table S5) |
| `upstream_scTHS-seq.html.txt` | scg_lib_structs page | second account, checked below |

Not fetched (PMC author-manuscript supplements behind a download gate; probably the same
content as the Springer ESM files, which were obtained):

| File | URL |
|---|---|
| NIHMS920808-supplement-1.doc | https://pmc.ncbi.nlm.nih.gov/articles/instance/5951394/bin/NIHMS920808-supplement-1.doc |
| NIHMS920808-supplement-2.docx | https://pmc.ncbi.nlm.nih.gov/articles/instance/5951394/bin/NIHMS920808-supplement-2.docx |
| NIHMS920808-supplement-3.pdf | https://pmc.ncbi.nlm.nih.gov/articles/instance/5951394/bin/NIHMS920808-supplement-3.pdf |

---

## 1. What it is

An ATAC-like accessibility assay in which the Tn5 adaptor carries a **T7 promoter**, so
that every insertion is amplified **linearly by in vitro transcription (IVT)** rather than
exponentially by PCR between two insertions; single-cell resolution comes from
**two-round combinatorial indexing** (sci-ATAC-style): a 384-plex barcoded transposome
round, FACS redistribution of ~100 nuclei per well, then an i5 × i7 PCR round (8 × 12 per 96-well plate; 16 i5 and 24 i7 primers in
Table S12). 🟢
(main text and Methods; Fig. S2)

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **T7 promoter inside the Tn5 transposon** → IVT amplifies each insertion end, independent of a second insertion nearby | bulk THS-seq (Sos 2016); IVT amplification is the CEL-seq / MARS-seq idea applied to genomic DNA |
| 2 | **Homodimer transposome**: one adaptor type only, so no 50 % loss to s5/s5 and s7/s7 ends | standard Nextera loads s5 and s7 at random — [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |
| 3 | After IVT, a **second, s7-only Tn5** fragments the re-made dsDNA and adds the 3' adaptor | the "tagment the cDNA" step of Smart-seq2 / CEL-seq2-type RNA-seq |
| 4 | The **cell barcode sits between the T7/i5 connector and s5**, read in the i5 index read | sci-ATAC-seq puts its Tn5 barcode next to s5/s7 as well, but reads it differently |
| 5 | Engineered hyperactive **Tn5059** (Illumina) | Kia et al. 2017 |

The rationale (THS-seq paper) is that in ATAC-seq random adaptor orientation leaves only
half of the molecules amplifiable, and regions where two neighbouring cuts are too far
apart cannot be amplified by PCR; making
every single insertion amplifiable by IVT recovers small and distal accessible sites. 🟢
(paraphrased)

## 2. Oligos

🟢 Verbatim from **Supplementary Table S12** ("DNA sequences for the customized
transposons and related adaptors"; `MOESM11_ESM.xlsx`, sheet "Table S12"). All written
5'→3'; `/5Phos/` = 5' phosphate (IDT notation as given).

```
scTHS-seq RNA-seq and i7 transposon oligos
sss_scnXTv2                    GGGAGATCCACGCGC
nXTv2_i7_top                   GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
P-T7tspn-bot, P-nXTv2_i7_bot   /5Phos/CTGTCTCTTATACACATCT

Sequencing primers
nXTv2_read1                    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
nXTv2_i7_index_read            CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
nXTv2_i5_index_read            AATGATACGGCGACCACCGAGATCTACAC

Barcoded r5 transposon oligos (384; first, second and last shown)
scT7_r5001_i5_top   AATTAATACGACTCACTATAGGGAGATCCACGCGCTCTAATTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
scT7_r5002_i5_top   AATTAATACGACTCACTATAGGGAGATCCACGCGCACTCGTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
scT7_r5384_i5_top   AATTAATACGACTCACTATAGGGAGATCCACGCGCGTTCAGTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
  general form:     AATTAATACGACTCACTATAGGGAGATCCACGCGC <r5 bc, 6 nt> TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG

Modified i5 barcoded primers (16; first and last shown)
scT7_S502   AATGATACGGCGACCACCGAGATCTACACCTCTCTATGGGAGATCCACGCGC
scT7_S522   AATGATACGGCGACCACCGAGATCTACACTTATGCGAGGGAGATCCACGCGC
  general form: AATGATACGGCGACCACCGAGATCTACAC <i5, 8 nt> GGGAGATCCACGCGC

Used original i7 Nextera XT V2 barcoded primers (24; first and last shown)
N701   CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG
N729   CAAGCAGAAGACGGCATACGAGATGACGTCGAGTCTCGTGGGCTCGG
  general form: CAAGCAGAAGACGGCATACGAGAT <i7 as written, 8 nt> GTCTCGTGGGCTCGG
```

The 16 i5 primers are S502, S503, S505–S508, S510, S511, S513, S515–S518, S520–S522; the 24
i7 primers N701–N707, N710–N712, N714–N716, N718–N724, N726–N729. 🟢 Random hexamers for RT
are used (20 µM) but **their sequence is not given** — plain N6 or with a tail is unknown. 🔴

### How they interlock — 🟡 (computed with `lib/`, from all 430 rows of Table S12)

**The r5 transposon (74 nt, all 384)** = `AAT` + T7 promoter (`TAATACGACTCACTATAG`,
positions 4–21, the last G being the +1 transcription start) + `GGAGATCCACGCGC` + 6-nt r5 barcode +
**`nextera.READ1_PRIMER`** (= `nextera.S5 + nextera.ME`, 33 nt). (3 + 18 + 14 + 6 + 33 = 74.) Checked: every one of the
384 has exactly this prefix and suffix; 384 distinct barcodes, minimum pairwise Hamming
distance 2 (the paper says "minimum edit distance of 2" 🟢 — consistent).

**The "i5 connector" is the start of the transcript.** Positions 21–35 of the transposon,
i.e. the T7 +1 G through the base before the barcode, are exactly **`sss_scnXTv2`
(`GGGAGATCCACGCGC`, 15 nt)**. So the 5' end of every IVT product is `GGGAGAUCCACGCGC…`, the
second-strand primer is identical to that 5' end, and the **3' end of every scT7_S5xx
primer is that same 15 nt**. One sequence is thus (a) the RNA 5' end, (b) the second-strand
primer, (c) the PCR priming site, and (d) the first 15 cycles after the i5 in the index
read.

**scT7_S5xx (52 nt)** = **`illumina.P5`** + 8-nt i5 + `sss_scnXTv2` — a Nextera XT v2
S5xx primer with s5 replaced by the connector. The i5 is written **as named** (forward),
identical to Illumina's S5xx set for all 16 (compared with Illumina's published S5xx
indices, which are not in the sources). 🟡

**N7xx (47 nt)** = **`nextera.n7xx_primer(i7)`** (`illumina.P7` + i7 + `nextera.S7`) for all
24 — unmodified Nextera XT v2. The i7 is written as the **reverse complement** of the
Illumina index name (`N701` carries `TCGCCTTA` = revcomp `TAAGGCGA`), all 24. 🟡

**nXTv2_i7_top** = **`nextera.ADAPTOR_S7`** exactly (34 nt); **P-T7tspn-bot /
P-nXTv2_i7_bot** = **`nextera.ME_RC`** with a 5' phosphate — **one bottom oligo serves both
transposomes**. 🟡

**Sequencing primers.** `nXTv2_read1` = `nextera.READ1_PRIMER`; `nXTv2_i7_index_read` =
`nextera.INDEX1_PRIMER`; `nXTv2_i5_index_read` = **`illumina.P5`** itself — not the Nextera
`INDEX2_PRIMER` (ME_RC + S5_RC). It anneals to the P5' end of the bottom strand and extends toward the insert:
i5 (8) → connector (15) → r5 barcode (6) → `TCG` of s5 (3) = **32 cycles**, matching the 32
index-2 cycles in the methods and Fig. S2C ("sequencing progressing through the linker
region and ending with the r5 barcode"). 🟡

**Bulk THS-seq ancestor** (Sos 2016, Additional file 1 Table S5) uses a different transposon,
`CATGAGATTAATACGACTCACTATAGGGAGATCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG` 🟢 (69 nt): T7
promoter at position 9, then a non-Nextera read-primer sequence and ME, with unbarcoded
one unbarcoded P5-tailed PCR primer plus barcoded P7 index primers (`Nxta_Ind…`, 8
variable nt in the primer) and 6-nt single-index sequencing (50 × 6). scTHS-seq re-plumbed this onto Nextera
XT v2 s5/s7 so standard Illumina read and i7 primers work. 🟡 (PDF-to-text separated the names column from the sequences, so only the transposon is
quoted)

## 3. Step by step

Conditions 🟢 from the Methods ("scTHS-seq transposon generation" through "library
validation, pooling and sequencing").

1. **Nuclei**: FACS-sorted nuclei, 5 min in lysis buffer (10 mM Tris pH 7.4, 10 mM NaCl,
   3 mM MgCl₂, 0.1 % NP-40, 2 % BSA, protease inhibitor, in PBS), then resuspended in
   1.5× tagmentation buffer (1×: 33 mM Tris-OAc pH 7.8, 66 mM K-OAc, 10 mM Mg-OAc, 16 % DMF).
   Human and mouse nuclei always mixed 1:1 for collision QC; ~2.4 × 10⁵ nuclei/mL.
2. **Anneal transposons**: each scT7_r5 top oligo + the phosphorylated ME bottom, 50 µM each,
   95 °C 2 min, cool to 14 °C at 0.1 °C/s; dilute to 8.4 µM in 50 % glycerol.
3. **Load Tn5059** fresh each run (activity decays over weeks): 1 µL 4.2 µM Tn5059 + 1 µL
   8.4 µM annealed r5 transposon per well of a 384-well plate, RT 30 min. Every well holds
   a **homodimer**: both arms carry the same T7–connector–r5bc–s5–ME transposon. 🟡
4. **Round 1, tagment**: 4 µL nuclei (~960, optimally ~2000 per well) into each of 384
   wells, 0.7 µM complex, 37 °C 30 min. Stop: 4 µL 50 mM EDTA, 37 °C 15 min; −20 °C
   overnight. Each insertion leaves the **74-nt top strand joined to the genomic 5' end**,
   with only the 19-nt ME bottom strand paired and the usual 9-nt gap — see
   [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). Both fragment ends carry the
   same transposon. 🟡
5. **Pool and redistribute**: add cold 2× FACS buffer, pool, spin, resuspend, PI stain, sort
   **100 nuclei per well** into 96-well plates (10 µL PBS); doublets gated out.
6. **Release DNA**: 11 µL guanidine hydrochloride per well, 1.8× AMPure, ethanol washes;
   beads stay in.
7. **Gap fill / end fill-in**: 10 µL 1× NEB Taq, 72 °C 3 min. This copies the transposon
   onto the bottom strand, which is what makes the T7 promoter **double-stranded** and
   therefore transcribable. 🟡 (the paper says only "end fill in")
8. **IVT**: NEB HiScribe T7 High Yield kit (10× buffer + ATP/CTP/GTP/UTP, 2 µL each) added
   directly, **37 °C 19 h**; checked on a TBU gel. Transcription starts at the +1 G of each
   end's promoter and runs through the insert to the far end (run-off), so each
   fragment yields aRNA from both ends, each copy beginning `GGGAGAUCCACGCGC` + r5 bc +
   s5 + ME + genomic sequence. 🟡 (computed transcript start; direction is the T7
   mechanism, not stated). Cleanup 2.0× SPRI buffer (20 % PEG 8000, 2.5 M NaCl), elute 9 µL.
9. **RT**: 2.5 µL 20 µM random hexamers, 70 °C 3 min, ice; Clontech SMART MMLV RT kit
   (first-strand buffer, dNTP, DTT), 22 °C 10 min, 42 °C 60 min, 70 °C 10 min. Random
   priming means cDNA 3' ends lie anywhere along the aRNA; the useful cDNA copies to the
   aRNA 5' end and so ends in the complement of the connector. 🟡 See
   [reverse transcription](../ref/concepts/reverse-transcription.md).
10. **RNase H** 0.5 U, 37 °C 20 min.
11. **Second strand**: 2.5 µL 20 µM **sss_scnXTv2**, 65 °C 2 min, ice; 5.9 µL NEB Taq 5×
    master mix, 72 °C 8 min. 2.0× SPRI, elute 7 µL. The primer anneals to the 3' end of the
    cDNA (complement of the aRNA 5' end) — product: ds `connector · r5bc · s5 · ME · gDNA…`.
    🟡
12. **Round 2 Tn5, s7 only**: pre-made custom **nXTv2_i7 Tn5059** homodimer (7 µM Tn5059 +
    10 µM annealed nXTv2_i7_top / ME-bottom, RT 30 min, diluted to 0.7 µM); 2 µL 5×
    tagmentation buffer + 2 µL complex (0.14 µM final), **55 °C 6 min**. Stop with 19 µL
    6.32 M guanidine-HCl (4 M final), 2.0× SPRI, elute 16 µL.
13. **Round 3, index PCR (qPCR)**: KAPA SYBR Fast 20 µL + 2 µL 10 µM **scT7_S5xx** + 2 µL
    10 µM **N7xx** ("nXTv2_i7XX" in the methods), 40 µL; 72 °C 3 min (fills the second Tn5 gap), 95 °C 30 s, then
    (95 °C 10 s, 63 °C 30 s, 72 °C 1 min) until saturation, typically 9–12 cycles; one
    8 × 12 i5 × i7 pair per 96-well plate.
14. **Pool and size-select** "as described" in the THS-seq paper; TBE gel check, Qubit.

Which round-2 products amplify 🟡 *(reasoned from the primers)*: only molecules with the
**connector at one end and s7 at the other**. Fragments with s7 at both ends lack a
scT7_S5xx site; a fragment with the connector end alone and no Tn5 cut has no N7xx site.
Each library molecule therefore carries **one original round-1 insertion** (its connector
end) and a round-2 cut somewhere downstream — the read is anchored at the genomic Tn5 site,
like a DNase-seq tag, not a paired fragment.

## 4. Final library — 🟡 (assembled from the oligos above)

```
5'- P5 · i5 (8) · connector GGGAGATCCACGCGC (15) · r5 bc (6) · s5 · ME · <gDNA, starting at the round-1 Tn5 site> · ME' · s7' · i7' · P7' -3'
```

As sequence (top strand; N = i5, B = r5 barcode, X = genomic, I = i7 as read):

```
AATGATACGGCGACCACCGAGATCTACAC NNNNNNNN GGGAGATCCACGCGC BBBBBB TCGTCGGCAGCGTC AGATGTGTATAAGAGACAG XXX…XXX CTGTCTCTTATACACATCT CCGAGCCCACGAGAC IIIIIIII ATCTCGTATGCCGTCTTCTGCTTG
```

Cell identity = **r5 barcode (384) × i5 (16) × i7 (24)**: round 1 in the 384 tagmentation
wells, round 2 in the sort plate (i5 × i7 identify the well of 100 nuclei). No UMI; clonal
reads were removed with samtools 1.3.1 🟢 (data processing section) — presumably by
mapping position, since there is nothing else to collapse on 🟡.

## 5. Sequencing

🟢 **Single-end**: MiSeq **50 + 32 + 32** for validation, HiSeq 2500 **50 + 8 + 32** for
data (Read 1 + Index 1 + Index 2). Read 1 = genomic DNA from the round-1 Tn5 site
(primer `nXTv2_read1` = Nextera Read 1); Index 1 = i7 (`nXTv2_i7_index_read` = Nextera
Index 1 primer), 8 cycles; Index 2 = i5 + connector + r5 barcode in one 32-cycle read
primed by `nXTv2_i5_index_read` = P5 (§2). Demultiplexing with bcl2fastq + deindexer,
0 mismatches. There is **no Read 2**.

## 6. Check against upstream scg_lib_structs

| Point | Upstream page | Paper / Table S12 | Verdict |
|---|---|---|---|
| Oligo sequences (r5 transposon form, ME bottom, nXTv2_i7_top, sss_scnXTv2, S5xx/N7xx forms, read1, i7 index read) | as in §2 | as in §2 | **agree** 🟢 |
| r5 barcodes / i5 primers | 384 × 6 nt / 16 | 384 / 16 (+ 24 N7xx) | agree |
| Nuclei per well | ~960, then ~100 | same | agree |
| aRNA 5' end | drawn as `GAGAUCCACGCGC…` | T7 +1 is the G of `TATAG`, so `GGGAGAUCC…` (= sss_scnXTv2) | **upstream drops `GG`** 🟡 (computed) |
| Read layout | paired-end, with a Read 2 step | single-end 50 + 8 + 32, no Read 2 | **disagree** — no Read 2 in the paper |
| i5 read primer | drawn annealing at ME/s5 and reading toward P5 (Nextera Index 2 style) | Table S12 lists `nXTv2_i5_index_read` = P5, reading from P5 toward the insert; 32 cycles | **disagree on primer**, agree on cycle count and content (8 + 15 + 6 + 3) |
| RT primer | "random hexamer, maybe with an overhang" | random hexamers, no sequence | both unknown 🔴 |
| Leftover heading | "Step-by-step generation of oligo-dT19V" | — | copy-paste slip in upstream; harmless |

## 7. Open questions

- 🔴 Exact sequence of the **random hexamer** RT primer (plain N6 or tailed).
- 🔴 Why **MiSeq** used 32 cycles for *Index 1* (50 + 32 + 32) when the i7 is 8 nt —
  possibly just reading on into s7/ME; not explained.
- 🔴 What the 5' `AAT` before the T7 promoter (`TAATACG…`) does (likely a landing pad so T7 RNAP
  binds a promoter not at the very DNA end 🟡).
- 🟡 Whether IVT run-off from the two ends of one short fragment interferes (upstream
  speculates one strand dominates) — not addressed by the paper.
- 🟡 Size selection conditions are only "as described" in Sos 2016; not re-read here.
- 🔴 PMC-only supplements (NIHMS920808-supplement-1..3) not fetched; assumed to duplicate
  the Springer ESM.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` had already fetched the paper (PMC HTML/XML), the Springer ESM
files and the THS-seq paper with its supplement; the three PMC author-manuscript
supplements were behind a download gate (`(manual)` in the manifest). `tools/scrape_primers.py`
ranked Table S12 first and recognised `nextera.ADAPTOR_S7`, `ME_RC`, `READ1_PRIMER`. It
found nothing in the THS-seq HTML (the oligos are only in Additional file 1, where
PDF-to-text separated names from sequences). Relationships across all 430 rows were
computed with a short script over the Table S12 text and `lib/` (`nextera`, `illumina`,
`chemdraw.revcomp`); every 🟢 sequence above was confirmed with
`scrape_primers.py --find`.
