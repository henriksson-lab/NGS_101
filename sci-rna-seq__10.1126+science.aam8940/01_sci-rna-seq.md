# sci-RNA-seq (2017): two-level combinatorial indexing of fixed cells or nuclei

> **Evidence marking.** 🟢 read in an original source · 🟡 derived or inferred.
> An inaccessible source is an open review obligation, not evidence of publication absence. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note.
> Upstream drawings are orientation aids, never source evidence. The 2017 public
> schematic is scoped to the two-level preprint implementation; the linked 2019
> successor supplies no 2017 oligo or reagent choices.

**sci-RNA-seq**: Cao J, Packer JS, Ramani V, Cusanovich DA, Huynh C, Daza R, Qiu X,
Lee C, Furlan SN, Steemers FJ, Adey A, Waterston RH, Trapnell C, Shendure J.
"Comprehensive single-cell transcriptional profiling of a multicellular organism."
*Science* 357:661–667 (2017). doi:[10.1126/science.aam8940](https://doi.org/10.1126/science.aam8940),
PMC5894354. Preprint: *bioRxiv* doi:[10.1101/104844](https://doi.org/10.1101/104844)
(v1, 2 Feb 2017; CC-BY-NC-ND 4.0).

**sci-RNA-seq3**: Cao J, Spielmann M, Qiu X, Huang X, Ibrahim DM, Hill AJ, Zhang F,
Mundlos S, Christiansen L, Steemers FJ, Trapnell C, Shendure J. "The single-cell
transcriptional landscape of mammalian organogenesis." *Nature* 566:496–502 (2019).
doi:[10.1038/s41586-019-0969-x](https://doi.org/10.1038/s41586-019-0969-x), PMC6434952.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/sci-rna-seq__10.1126+science.aam8940/`, never committed):

| File | What | Used for |
|---|---|---|
| `104844_v1.full.pdf` (.txt) | sci-RNA-seq preprint, with Materials & Methods | selected preprint two-level RT/PCR templates, reagents, order and read lengths |
| `science.aam8940_PMC5894354.html` (.txt) | Science author manuscript on PMC | workflow summary, Fig. 1 legend and the three-level indexed-Tn5 variant; the detailed selected implementation is printed in the preprint |
| `cao-combinatorial-indexing.pdf` (.txt) | Author-hosted final Science version of record | final Fig. 1, two-level/optional indexed-Tn5 branch, original DC1 supplementary inventory |

Could not be fetched (get by hand and drop in the same directory):

| What | URL | Why it matters |
|---|---|---|
| Science supplement (NIHMS956226, Supplementary Materials, 2.2 MB PDF) | https://pmc.ncbi.nlm.nih.gov/articles/instance/5894354/bin/NIHMS956226-supplement-supplement_1.pdf | the actual final Supplementary Materials, named by final-article DC1 inventory; its contents and index assignments remain to inspect, rather than being assumed from the filename (unresolved here) |
| Science supplementary landing | https://www.science.org/doi/suppl/10.1126/science.aam8940 | publisher route returns 403; historical DC1 also403; supplementary workbook DOI discovered in later bibliography returns404 and is not a validated recovered original; final main is recovered from the authors |

---

## 1. What it is

No compartments, no beads: **the cell (or nucleus) itself is the reaction vessel**, and its
identity is the combination of well-specific barcodes it collected while being split,
pooled and split again. 🟡 (paraphrase: the Science main text contrasts it with methods
that isolate cells in physical compartments; workflow per Fig. 1A)

- **sci-RNA-seq** (2017), two levels: methanol-fixed cells or nuclei are split into
  96/384 wells for **in situ RT with a barcoded oligo-dT** (barcode 1 + UMI), pooled,
  **FACS-sorted** 10–100 per well into a second plate, and given second-strand synthesis,
  Tn5 tagmentation and **indexed PCR** (barcode 2 = i5 + i7 pair). 🟢 An optional third
  level puts a barcode on **indexed Tn5** between RT and PCR (16 × 6 × 16 demonstrated). 🟢
  (Science main text; indexed-Tn5 original supplementary details remain inaccessible)

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Combinatorial indexing of whole fixed cells**: the first barcode is put on in situ, the cell stays intact through pooling | sci-ATAC-seq, sci-CAR, SPLiT-seq / SHARE-seq (ligation rounds), sci-Plex |
| 2 | **Read 1 handle on the RT primer, s7 from Tn5 at the other end**: only the 3' end of each cDNA, between the RT primer and the nearest s7 insertion, amplifies | the same 3'-tag logic as CEL-seq / Drop-seq, built from TruSeq + Nextera halves |

Concepts: [reverse transcription](../ref/concepts/reverse-transcription.md) (oligo-dT
priming, RNase H / Pol I second strand),
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). No template switching: the 5' end
of the transcript is never captured.

## 2. Oligos

### sci-RNA-seq (2017) 🟢

Verbatim from the preprint Materials & Methods (`104844_v1.full.pdf.txt` lines 899–900
and 922–931; the PDF wraps them, re-joined here). Brackets are the paper's placeholders.

```
anchored oligo-dT (RT)  5′-ACGACGCTCTTCCGATCTNNNNNNNN[10bp index]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′   (IDT)
P5 primer               5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′
P7 primer               5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′
Tn5                     Nextera DNA Sample Preparation kit, TDE1 (commercial intermediate structure unassigned)
```

"N" is any base, "V" is A, C or G. 🟢 The templates are positively printed in
preprint Methods; the GEO annotations now positively supply 96/384/96 RT barcode memberships and
cell-source groups for three experiments. Individual PCR-well i5/i7 assignments and
comparison with the original final Science supplementary material remain open. This is
not a claim that the authors failed to publish them. The same RT-primer template is
printed for the bulk Tn5-RNA-seq comparison (preprint Methods).

The source-selected two-level model keeps the RT barcode and the PCR-well i5/i7
pair as **parts of one combinatorial cell identifier**, rather than interpreting i5/i7
as ordinary sample indexes. The model preserves both V and N of the RT-primer VN
anchor. The canonical Nextera mosaic-end bases are explicitly inferred; the printed
PCR-primer arms themselves are source-positive. 🟡 The exact pre-PCR TDE1 substrate,
gaps and end chemistry are not assigned by this declaration.

The graph now includes the source's DNA-binding-buffer stop and AMPure purification
before library PCR. Carrier and other tagmentation products remain visible with
unresolved structures/fates: the source says carrier is **not appreciably amplified**,
which does not establish an immediate physical discard during tagmentation.

### Current 2017 source scope (2026-10-16)

The full bundle oligo scraper was rerun after retrieval (79 hits across 13 files);
2019 successor sources in that bundle were excluded from 2017 chemistry assignments.
The author-hosted final PDF, pp. 1–2/Fig. 1 and final supplementary inventory, and
preprint Methods lines 899–940/991–992 were read directly. Final supplementary
inventory names Materials and Methods, Figs. S1–S24, Tables S1–S14 and refs. 57–140
at `www.sciencemag.org/content/357/6352/661/suppl/DC1`.

The registered publisher supplemental landing and Europe PMC original NIHMS956226
mirror return 403; the original PMC route is a preparing-download gate. PMC Cloud
lists only three accessible manuscript objects, not the required supplementary PDF.
The actual author publication page links the final article PDF and Science DOI;
the inspected GitHub `pdfs/papers` directory inventory names that final PDF, not an
identified 2017 supplement. The original bioRxiv v1 landing exposes its main PDF and
article tabs; no supplementary tab was found among the inspected links. These are
bounded retrieval findings, not evidence of publication absence.

The preprint positively identifies a NEB mRNA Second Strand Synthesis module, Nextera
DNA Sample Preparation/TDE1, SSIV and NEBNext HiFi PCR. Exact relevant vendor-original
versions and any further predecessor chemistry in the inaccessible final supplement
remain unchecked linked-source obligations. GSE98561 is an actual final-article data
link. Its study metadata, annotation files and ten archive experiment XMLs were
retrieved and read in this pass; no dataset-wide index-map absence is asserted. The review remains blocked by
actual final supplementary access, and unassigned model values remain questions.

### How the 2017 primers interlock — computed

The RT handle is `illumina.TRUSEQ_READ1[-18:]`. The P5 PCR primer contains the
complete TruSeq Read 1 sequence and anneals through its terminal RT handle; PCR
therefore restores the full Read 1 site. The P7 PCR primer ends in `nextera.S7`.
The representative amplifiable strand must have the RT handle at one end and the
opposing s7 site at the other. This sequence relationship is computed locally;
it does not prove exact commercial transposome intermediates or every byproduct.

The later [sci-RNA-seq3 note](../sci-rna-seq3__10.1038+s41586-019-0969-x/01_sci-rna-seq3.md)
owns the distinct 2019 hairpin-ligation chemistry and oligo table. Those sources
are not evidence for a 2017 index member or transposome choice.

## 3. Step by step

### sci-RNA-seq (preprint methods) 🟢 unless marked

1. **Fix / isolate**: cells in ice-cold 100 % methanol, −20 °C 10 min, DEPC and
   SUPERase-In / BSA washes; or nuclei by IGEPAL lysis. 5000 cells or nuclei per µL.
2. **In situ RT, barcode 1**: per well 1,000–10,000 cells (2 µL) + 1 µL 25 µM barcoded
   oligo-dT + dNTPs, 55 °C 5 min, ice; SuperScript IV mix; **55 °C 10 min**; stopped with
   2× stop solution (40 mM EDTA, 1 mM spermidine). First strand (🟡):
   `5'- handle · UMI · RT barcode · T30 · V · N · <cDNA, antisense> -3'`.
3. **Pool, stain DAPI, FACS** 10–100 cells per well (an example range; per-experiment numbers are in Table S1 (unresolved here))
   into 5 µL EB; DAPI gating removes doublets.
4. **Second strand**: NEB mRNA Second Strand Synthesis module, 16 °C 150 min, 75 °C 20 min. The declaration does not assign a separate printed random primer.
5. **Tagment**: + **5 ng human genomic DNA as carrier** (against over-tagmentation and
   loss), Nextera TD buffer, 0.5 µL TDE1, 55 °C 5 min; Zymo DNA binding buffer, AMPure
   (36 µL), elute 16 µL. The methods state that carrier DNA is not appreciably amplified
   because the PCR primers are specific to the RT products 🟢; the mechanism (carrier has
   s5/s7 ends only, and the P5 primer needs the handle) is 🟡.
6. **Indexed PCR, barcode 2**: P5[i5] + P7[i7], NEBNext HiFi 2× master mix;
   **75 °C 3 min** (gap fill of the 9-nt Tn5 gap 🟡), 98 °C 30 s, **18–22 ×** (98 °C 10 s,
   66 °C 30 s, 72 °C 1 min), 72 °C 5 min. Pool, AMPure 0.8×, Qubit, 6 % TBE-PAGE.

   Which fragments amplify (🟡 *(computed)* from the two PCR-primer sites): only the one
   from the **RT-primer end to the nearest s7 insertion**. Handle-to-s5 fragments, s5–s5,
   s7–s7 and s5–s7 internal fragments (and all carrier DNA) lack one of the two primer
   sites, since the P5 primer needs the TruSeq handle and does not bind s5.
7. **Three-level variant** (Science only): after RT and second-strand synthesis, pool,
   split to 6 wells for **indexed Tn5**, pool, sort to 16 wells for PCR. 🟢 Adaptor
   sequences, where the Tn5 barcode lands in the read, and the sequencing recipe are in
   the Science supplement (unresolved here). The paper notes fewer UMIs per cell than two-level at equal
   depth, possibly from lower efficiency of indexed Tn5. 🟢

## 4. Final library and read layout

The representative top strand is reconstructed as:

```
5′-P5 · i5 · TruSeq R1 · UMI · RT barcode · T30 · V · N · antisense cDNA · inferred ME_RC · S7_RC · i7_RC · P7_RC-3′
```

The exact available RT/PCR templates are source-positive; the intervening canonical
Nextera mosaic end is inferred. Index placeholders express unassigned members,
not a claim that all ordered index combinations were used. cDNA length is an
illustrative display value, not a source-measured insert length.

The preprint explicitly names NextSeq 500 V2, Read 1=18, Index 1=10, Index 2=10,
Read 2=52 (actual extracted lines 991–992). Read 1 covers UMI8 then RT barcode10.
Index reads supply the PCR-well pair; Read 2 supplies transcript sequence. The
model locates each declared primer computationally. This layout is the preprint
implementation, with no claim that it proves the optional indexed-Tn5 branch's run.

## 5. Positive GEO source information

The original Series metadata names three small barcode annotation files. Its
`GSE98561_README_supplementary_files.txt` explicitly identifies their `barcode`
column as the **first-round RT barcode** and `source` as the associated cell group.
These are source annotations, not gene-count matrices or raw sequencing reads.

| Original file | RT members | Actual source scope |
|---|---|---|
| `GSE98561_celegans.barcode.csv.gz` | 96 unique 10-nt barcodes | C. elegans cells |
| `GSE98561_human.mouse.mix.1.barcode.csv.gz` | 384 unique 10-nt barcodes | pure and mixed cell groups |
| `GSE98561_human.mouse.mix.2.barcode.csv.gz` | 96 unique 10-nt barcodes | cell groups plus mixed nuclei |

Both readers inspected all rows and source-group fields, rather than treating the
scraper's lack of hits for short 10-nt sequences as absence. The model deliberately
retains the source's variable index windows; it does not assign a particular RT/PCR
combination to every experiment.

All ten GEO sample library/experiment/instrument metadata blocks directly cite the
2017 preprint protocol. Nine samples name NextSeq 500; the explicitly three-level
sample GSM2599705 names MiSeq. The ten linked archive experiment XMLs confirm these
instrument assignments and paired transcriptomic cDNA libraries. They do not assign
run cycle counts or PCR primer members in the inspected fields. Thus the modeled
18/10/10/52 run is the **preprint two-level implementation**, not an assertion about
every final study sample or the MiSeq third-level experiment.

The source's other README maps 1,717 run accessions to paired FASTQ filenames. It
is a file-access mapping, not the ordered PCR-primer table. The archive XMLs link
sample/submission records and file reports; these data references have not been
recursively inspected, and no universal archive absence is inferred.

The recovered Illumina reference guide 15027987 v01, January 2016, positively
describes the Nextera DNA engineered-transposome and index-PCR family. Its standard
indexes are 8 nt, whereas the sci-RNA-seq custom primer/run declarations have 10-nt
index windows. The guide therefore supplies family context, not the authors'
ordered index members or the selected commercial adapter's exact intermediate.

## 6. Remaining source and model questions

- The final Science supplement is an identified, unrecovered original. Its Methods,
  Figs. S1–S24, Tables S1–S14 and refs. 57–140 remain to inspect. The actual historical
  DC1 link and candidate direct publisher file paths return403; the original NIHMS
  route is gated and the canonical article/bin alternative returns404. A workbook
  DOI discovered through later bibliography returns404. These route failures do
  not prove unpublished information or validate a guessed filename.
- The final-article Waterston data-directory link is problematic: live HTTP and
  HTTPS downloads return an unrelated EPIC/Lineagomics RNAi page. Both wrong-identity
  artifacts were quarantined. The browser's cached primary directory lists a 2017
  vignette, sample description, count data and Monocle software, but the live
  sample-description CSV request also returns HTML. This is a retrieval/identity
  barrier, not a claim that the original directory lacked information.
- Exact historical NEB second-strand and NEBNext PCR vendor originals remain
  unchecked. The preprint positively names those reagents; no current vendor default
  fills their unassigned chemistry. The recovered Nextera guide's additional
  protocol/checklist/low-plex references remain named source-scope questions.
- The selected endpoint's canonical ME and exact gap-filled intermediate remain
  inferred. Tagmentation carrier/byproduct structures and physical fates are
  unresolved, rather than being rendered as a fictitious physical discard.

The source review remains blocked. Two independent fresh readings confirm the
preprint's full templates, reagent order and run, the final article's optional
indexed-Tn5 distinction, and the three GEO RT-membership annotations. This does not
close the inaccessible final supplement or the remaining linked-source inventory,
and accepts no claim of publication absence.
