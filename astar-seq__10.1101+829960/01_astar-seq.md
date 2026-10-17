# ASTAR-seq — chromatin accessibility and transcriptome from the same single cell

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/`; the model
> is `tools/astar.py` and its self-test (`tools/selftest.py`) pins the transcribed oligos.
>
> **Provenance.** Every oligo below is the authors' own, printed identically in the
> preprint's Supplementary Table 5 and in the journal version's Supplemental Table 5.
> `tools/selftest_source_oligos.py` re-finds the three C1 oligos in the locally extracted
> text of those tables, so the model's copies are machine-checked against the documents
> rather than against this note. Every other concrete base in the model is a canonical
> `lib/illumina.py` or `lib/nextera.py` constant; the i5/i7 windows are placeholders.
> **The scg_lib_structs drawn page is not a source here** — it is inspiration only, no
> copy of it was fetched for this protocol (the sources directory contains no
> `upstream_*` file), nothing is taken from it, and agreement with it would not be
> evidence.

**ASTAR-seq** — Xing QR, Farran CE, Yi Y, Warrier T, Gautam P, Collins J, Xu J, Li H,
Zhang L-F, Loh Y-H. "Parallel bimodal single-cell sequencing of transcriptome and
chromatin accessibility." *bioRxiv* 2019. doi:[10.1101/829960](https://doi.org/10.1101/829960)
(v1, 4 Nov 2019; CC-BY-NC 4.0). Data: GEO GSE113418.

**Published version (authoritative):** Xing QR *et al.* "Diversification of reprogramming trajectories revealed by
parallel single-cell transcriptome and chromatin accessibility sequencing." *Genome
Research* 30:1027–1039 (2020).
doi:[10.1101/gr.257840.119](https://doi.org/10.1101/gr.257840.119) · PMC7397874.
The published article and its Supplemental Table 5 are the authority. Its main
Methods (JATS lines 162–166) give only the on-/off-IFC step list and sequencing; the
detailed Supplemental Methods PDF it cites could not be retrieved, so exact reaction
conditions below come from the bioRxiv v1 supplement (`media-6`) and are labelled
as preprint values. 🟡 (explicit preprint use for omitted detail)

Sources read (fetched by `tools/get_sources.py astar-seq`, into
`$CHEM_DATA/sources/astar-seq__10.1101+829960/`, never committed):

| File | What | Used for |
|---|---|---|
| `829960_v1.full.pdf` | the preprint | workflow, the failed prototype |
| `829960_media-6.pdf` | Supplementary Information | **all methods** (the main text has none) |
| `829960_media-5.xlsx` | Supplementary Table 5 | **every oligo** |
| `829960_media-1..4.xlsx` | QC tables, motif enrichment | not chemistry |
| `gr.257840.119_PMC7397874.xml` | the journal version's full text (PMC Cloud JATS) | sequencing parameters, published wording |
| `gr.257840.119_supp_..._Supplemental_Table_5.xlsx` | the journal version's **Supplemental Table 5** | independent second printing of **every oligo** — checked row by row against the preprint table |
| `gr.257840.119_supp_..._Supplemental_Material_.docx`, `..._Protocol_Worksheet_ASTAR*.pdf` | supplemental text and the C1 script worksheets | no oligo sequences |

---

## 1. What it is

A plate-free, **Fluidigm C1** method: each cell sits in its own microfluidic chamber, so
there is **no cell barcode and no UMI** anywhere in the chemistry. Cell identity is the
chamber, carried into the library by a per-cell **Nextera dual index** at the very end.
One cell gives two libraries, ATAC and RNA, separated by **biotin** before either is
amplified to completion. 🟢 (Supplementary Information, "ASTAR-seq METHOD")

The order is the point of the method 🟢: **tagment first, then reverse transcribe.** The
earlier prototype did RT first and failed, because Tn5 cut the single-stranded cDNA too,
leaving cDNA fragments with ATAC adapters that could not be separated from the ATAC-DNA
(main text; Supplementary Fig. 1a–d).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | Two modalities from one cell separated by **biotin pull-down**, not by priming | G&T-seq separates by oligo-dT beads; SHARE-seq / SNARE-seq by barcoding both |
| 2 | **Tn5 before RT**, with EDTA to kill Tn5 and MgCl₂ to quench the EDTA | the Mg²⁺ bookkeeping recurs in every one-pot multi-enzyme protocol |
| 3 | A **custom C1 handle** in place of the SMART handle on oligo-dT, TSO and PCR primer | SMART-seq2 uses the ISPCR / SMART handle (`rt.SMART_HANDLE`) |

## 2. Oligos

🟢 Verbatim from Supplementary Table 5 (`media-5.xlsx`, sheet "Primers"), and identical in
the journal version's Supplemental Table 5. IDT notation as given: `/5BiosG/` = 5' biotin,
`rN` = ribonucleotide. Names are the table's own.

```
C1-P2-T31      /5BiosG/ GGCGACAACACCGATTGATCACG TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C1-P2-RNA-TSO  rGrGrCrGrArCrArArCrArCrCrGrArUrUrGrArUrCrArGrGrG
C1-P2-PCR-2    /5BiosG/ GGCGACAACACCGATTGATCA

ATAC adaptor 1 (forward)  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
ATAC adaptor 2 (reverse)  GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG

v2_Ad1.N_<i5>  AATGATACGGCGACCACCGAGATCTACAC <i5> TCGTCGGCAGCGTCAGATGTGTAT     N = 1..16, 61 nt
v2_Ad2.N_<i7>  CAAGCAGAAGACGGCATACGAGAT <revcomp i7> GTCTCGTGGGCTCGGAGATGTG     N = 1..12, 54 nt
```

(The methods call the first and third "Bio-C1-P2-T31" and "Bio-C1-P2-PCR-2"; the table
names them without "Bio-". Same oligos — the table shows the biotin. 🟡)

### How the three C1 oligos interlock — 🟡 (computed)

All three carry one 21-nt handle, **`GGCGACAACACCGATTGATCA`**:

- **C1-P2-PCR-2 *is* the handle**, biotinylated — 21 nt.
- **C1-P2-T31 = handle + `CG` + T₃₁** — 54 nt. No anchoring `VN`: it is a plain T31, so
  priming is not pinned to the poly(A) junction. The two-base `CG` spacer is not explained. 🔴
- **C1-P2-RNA-TSO = handle + `GGG`**, entirely RNA — 24 nt. Unlike the SMART-seq2 TSO
  (DNA with 3' `rGrG+G`), every base is ribo. 🟢 sequence; reason not given 🔴.

So after RT and template switching, **both ends of every cDNA carry the same handle**,
and a single primer (PCR-2) amplifies it — the SMART-seq2 single-primer design with a
different handle. Because that primer is 5'-biotinylated, **every amplified cDNA strand is
biotinylated at its 5' end**, which is what the streptavidin step sorts on.

### The ATAC indexing primers — 🟡 (computed against `lib/`)

- `v2_Ad1.N` = **`illumina.P5`** + 8-nt i5 + **`nextera.ADAPTOR_S5[:24]`** (s5 +
  the first 10 nt of the mosaic end). The i5 is written **as named** —
  `v2_Ad1.1_TAGATCGC` contains `TAGATCGC` — for all 16.
- `v2_Ad2.N` = **`illumina.P7`** + 8-nt i7 + **`nextera.ADAPTOR_S7[:22]`** (s7 + ME[:7]).
  The i7 is written as the **reverse complement** of its name — `v2_Ad2.1_TAAGGCGA`
  contains `TCGCCTTA` — for all 12. 16 × 12 = 192 index pairs: the 192 libraries per
  lane in §5. 🟡
- "ATAC adaptor 1 (forward) / 2 (reverse)" are exactly **`nextera.ADAPTOR_S5` /
  `ADAPTOR_S7`** — used here only as qPCR primers, to test whether cDNA had been tagmented
  (Supplementary Fig. 1c–d, which refers to "the ATAC adaptors (Supplementary Table 5)" 🟢).

The `Ad1.x` / `Ad2.x` naming and layout are those of the Buenrostro ATAC-seq indexing
primers; whether the 8-nt indices are identical to that set was not checked. 🟡

## 3. Step by step

All on the C1 IFC unless stated. The publisher's detailed
`Supplemental_Material.pdf` is indexed with the Part 1 on-IFC procedure, but its
direct URL was inaccessible during this audit. The cached similarly named DOCX is
only an inventory, not the Methods body. The author bioRxiv supplement
`829960_media-6.pdf.txt`, lines 1562–1665, contains the complete earlier recipe.
Exact values below are sourced to that preprint unless a published-version value
is stated. Full published off-IFC methods still need a direct check.

1. **Capture** one cell per chamber (C1 IFC, custom "ASTAR" scripts, Supplementary File 1).
2. **Lyse and tagment together**: 0.15 % NP-40, 1.5× Tagment DNA Buffer, 1.5× Nextera Tn5
   (TDE1), 37 °C 30 min. Open chromatin becomes Tn5 fragments: insert flanked by
   **ME, with s5 or s7** on either end and the usual 9-nt gap — see
   [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). mRNA is untouched.
3. **Kill Tn5**: EDTA (18.75 mM in the inactivation mix), 50 °C 30 min. The same mix
   brings **dNTPs and C1-P2-T31**; 72 °C 3 min / 4 °C 10 min / 25 °C 1 min anneals the
   oligo-dT.
4. **RT with template switching**: SuperScript IV + C1-P2-RNA-TSO, **MgCl₂ added to quench
   the EDTA** 🟢, 50 °C 60 min, 80 °C 10 min. First strand: `5'-biotin-handle-CG-T31-…cDNA…-CCC`,
   then the TSO's `GGG` pairs with the non-templated `CCC` and the RT copies the TSO,
   putting the handle (as its complement) on the 3' end. 🟡 (the standard
   [template-switching](../ref/concepts/template-switching.md) mechanism; the paper
   only says "double-stranded cDNA")
5. **On-chip cDNA PCR, 5 cycles** with biotinylated C1-P2-PCR-2 (Q5 Ultra II): 98 °C 3 min;
   5 × (98 °C 20 s, 58 °C 4 min, 68 °C 6 min); 72 °C 10 min. Every cDNA strand now has a
   5' biotin. The ATAC fragments are **not** amplified: no primer matches s5/s7. 🟡
6. **Harvest**: preprint Supplementary Experimental Procedures lines 1603–1605
   print a total of **5.5 µL** C1 Harvest Reagent and lines 1612–1615 mix 5.5 µL
   harvested sample with 5.5 µL washed streptavidin beads. 🟢 (preprint) An earlier
   audit recorded **approximately 3.5 µL** from a search-index snippet of the
   publisher's later `Supplemental_Material.pdf`; that PDF returned HTTP 404/403 on
   2026-10-17, so the later value is **unverified** and no model state depends on
   either volume. Whether any volume changed on publication remains open.

Off chip:

The preprint's original Supplementary Experimental Procedures lines 1606–1665
support the cycle counts and handling quantities in steps 7–10. They are verified
for the preprint version; the inaccessible later publisher PDF prevents a full
check of whether those exact off-IFC conditions changed on publication.

7. **Separate**: MyOne Streptavidin C1 beads retain biotinylated cDNA while the
   supernatant contains ATAC-DNA. The preprint Methods lines 1618–1621 specify
   elution of cDNA in 10 µL water at 90 °C for 10 min and transfer of the eluate. 🟢
   After on-chip PCR with the biotinylated primer both cDNA strands carry a 5'
   biotin, and the source names neither the release mechanism nor the strand form
   of the eluate; the model draws the eluted cDNA as inferred, without a strand
   assignment. (Corrected 2026-10-17: an earlier model released a
   "non-biotinylated strand", which on-chip PCR with one biotinylated primer would
   not leave.) 🟡
8. **cDNA, second PCR**: 15 µL Q5 Ultra II mix with 1.6 µM C1-P2-PCR2 (now given
   without "Bio-") on the 10 µL eluate; 9 cycles (64 °C 30 s, 68 °C 6 min) then 11
   cycles (68 °C 7 min) 🟢 (lines 1622–1628); then 20 µL AMPure XP, eluted in 10 µL
   🟢 (lines 1629–1632) — about 0.8× of the 25 µL reaction 🟡 (computed; the
   earlier "1×" in this note was not supported).
9. **RNA library**: standard C1 / Nextera XT on 0.15–0.2 ng/µL cDNA — tagment 55 °C 10 min,
   Nextera XT index PCR, 72 °C 5 min then 22 cycles; pooled (~1.1 mL) and purified
   twice with 0.75× AMPure XP 🟢 (lines 1636–1648). Only fragments with s5 at one end and s7 at the
   other amplify ([Tn5](../ref/concepts/tn5-tagmentation.md)); fragments that kept the
   C1 handle at one end and Tn5 at the other carry the handle **inside** the read, not
   as a primer site, because no handle-specific primer is used here. 🟡
10. **ATAC library**: the 10 µL supernatant + Q5 Ultra II with **one v2_Ad1.N and one
    v2_Ad2.N per cell**; 72 °C 5 min (gap fill-in) then 22 cycles. Pooled, ethanol
    precipitated, MinElute, AMPure 1.2×.

## 4. Final libraries — 🟡 (assembled from the oligos above)

**ATAC** (standard dual-indexed Nextera; the custom primers end inside the ME, and the
remaining ME bases come from the tagmented fragment itself):

```
5'- P5 · i5 · s5 · ME · <genomic insert> · ME' · s7' · i7' · P7' -3'
```

**RNA**: a Nextera XT library from full-length cDNA, indexed with the kit's i5/i7. Reads
start at random positions along the cDNA; a minority start in the C1 handle or poly(A).

## 5. Sequencing

🟢 HiSeq 4000. RNA: 192 libraries per lane, **2 × 101**. ATAC: 192 libraries per lane,
**2 × 50**. Standard Nextera read primers: Read 1 = `nextera.READ1_PRIMER`, Read 2 =
`nextera.READ2_PRIMER` 🟡 (not stated; implied by Nextera chemistry and no custom read
primers listed).

## 6. Open questions

- 🔴 Why the **`CG` spacer** between the handle and T31, and why the TSO is **all-RNA**.
- 🔴 Whether the second-round "C1-P2-PCR2" (no "Bio-") is a different, unbiotinylated
  oligo or the same one named loosely. Table 5 lists only the biotinylated one.
- 🟡 Identity of the index set with Buenrostro's — unchecked.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the preprint and all six supplements; `tools/doctext.py`
turned them into text; `tools/scrape_primers.py` on the directory put Supplementary
Table 5 first (its score ranks the biotinylated, RNA-modified oligos highest) and showed
all 16 i5 and the i7 primers with `illumina.P5` / `nextera.ADAPTOR_S7` already
recognised. The motif-enrichment tables (Homer, degenerate IUPAC) scored lowest and were
cut by `--max-hits`. The methods themselves had to be **read**, not scraped — the
scraper finds oligos, not reaction order.
