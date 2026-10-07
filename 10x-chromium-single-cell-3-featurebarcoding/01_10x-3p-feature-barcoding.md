# 10x Chromium Single Cell 3' Feature Barcoding — antibody tags and sgRNAs on the 3' v3 bead

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`).
> **Source situation:** this is a vendor kit with no defining paper. The 10x Feature
> Barcoding user guide (CG000185-family) and its oligo appendix were **not** fetched; the
> only complete account of the feature-barcode oligos is the upstream scg_lib_structs page,
> so most of §2 is 🟡 (secondary source). The 10x **gene-expression** user guides that
> were fetched confirm the bead oligos, TSO and index primer.

**10x Chromium Single Cell 3' Feature Barcoding** (10x Genomics, with the v3 3' chemistry,
2018; carried into GEM-X 3' v4, 2024). No defining publication. The scg_lib_structs page
cites as background:

- Stoeckius M *et al.* "Simultaneous epitope and transcriptome measurement in single
  cells." *Nat Methods* 2017. doi:[10.1038/nmeth.4380](https://doi.org/10.1038/nmeth.4380) — CITE-seq, the antibody-tag idea.
- Briner AE *et al.* "Guide RNA functional modules direct Cas9 activity and orthogonality."
  *Mol Cell* 2014. doi:[10.1016/j.molcel.2014.09.019](https://doi.org/10.1016/j.molcel.2014.09.019) — sgRNA scaffold module names.
- Sanjana NE, Shalem O, Zhang F. "Improved vectors and genome-wide libraries for CRISPR
  screening." *Nat Methods* 2014. doi:[10.1038/nmeth.3047](https://doi.org/10.1038/nmeth.3047) (GeCKOv2).
- Koike-Yusa H *et al.* "Genome-wide recessive genetic screening in mammalian cells with a
  lentiviral CRISPR-guide RNA library." *Nat Biotechnol* 2014. doi:[10.1038/nbt.2800](https://doi.org/10.1038/nbt.2800).
- Datlinger P *et al.* "Pooled CRISPR screening with single-cell transcriptome readout."
  *Nat Methods* 2017. doi:[10.1038/nmeth.4177](https://doi.org/10.1038/nmeth.4177) (CROP-seq).
- Mimitou EP *et al.* "Multiplexed detection of proteins, transcriptomes, clonotypes and
  CRISPR perturbations in single cells." *Nat Methods* 2019. doi:[10.1038/s41592-019-0392-0](https://doi.org/10.1038/s41592-019-0392-0) (ECCITE-seq; 10x 5').
- Frangieh CJ *et al.* "Multimodal pooled Perturb-CITE-seq screens in patient models define
  mechanisms of cancer immune evasion." *Nat Genet* 2021. doi:[10.1038/s41588-021-00779-1](https://doi.org/10.1038/s41588-021-00779-1) (preprint doi:10.1101/2020.09.01.267211).

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/10x-chromium-single-cell-3-featurebarcoding/`, never committed):

| File | What | Used for |
|---|---|---|
| `upstream_10xChromium3fb.html(.txt)` | scg_lib_structs page | **all feature-barcode oligos, steps, libraries, read layout** (secondary) |
| `upstreamdata_CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf` | 10x 3' v3 GEX user guide | bead oligos (dT, Capture Seq 1), TSO, cDNA primers, i7 index primer, read lengths |
| `upstreamdata_CG000731_ChromiumGEM-X_SingleCell3_ReagentKits_v4_UserGuide_RevA.pdf` | 10x GEM-X 3' v4 GEX user guide | v4 bead has **only Capture Sequence 1**; dual-index layout |
| `006726_v1.full.pdf`, `nmeth.3047_*MOESM647_ESM.pdf` | Sanjana 2014 (preprint, SI) | GeCKO sgRNA scaffold (DNA) |
| `nmeth.4380_PMC5669064.html` | CITE-seq | poly(A) ADT design, Nextera-handle ADTs |
| `s41588-021-00779-1_*`, `2020.09.01.267211_*` | Perturb-CITE-seq (paper, SI, Table 8) | a 10x 3' v3 *alternative*: poly(A) ADT + CROP-seq dial-out |
| `s41592-019-0392-0_*` | ECCITE-seq | 5' relative, only skimmed |
| `nmeth.4177_*`, `nbt.2800_*` | CROP-seq, Koike-Yusa | vector background, not read in detail |

Not fetched (listed `(manual)` in `MANIFEST.tsv`): Briner 2014 (Mol Cell, publisher page
only); the PMC supplements of nmeth.3047, nmeth.4380, s41588-021-00779-1 and
s41592-019-0392-0 (download gate). **Most important gap, not in the manifest at all:** the
10x *Feature Barcoding* user guides (e.g. CG000185 "Chromium Single Cell 3' Reagent Kits
v3 with Feature Barcoding technology for Cell Surface Protein", CG000184 "… for CRISPR
Screening", and the v4 counterparts CG000317-like documents) and the 10x "Guide RNA
specifications" technical note. Those are the primary source for every 🟡 oligo below. 🔴

---

## 1. What it is

The 3' v3 Gel Bead carries **three** primer species on the same bead, all with the same
16-nt cell barcode *pattern* and a 12-nt UMI 🟢 (CG000183, intro and appendix):

1. **TruSeq Read 1 handle + barcode + UMI + poly(dT)** — the ordinary 3' gene-expression primer.
2. **Nextera Read 1 (Read 1N) handle + barcode + UMI + Capture Sequence 1**.
3. **Read 1N handle + barcode + UMI + Capture Sequence 2**.

The gene-expression protocol uses only (1); Feature Barcoding uses (2)/(3) to prime from
any molecule that ends in the **reverse complement of a capture sequence**: an
oligo-conjugated antibody (TotalSeq-B type) or an sgRNA engineered to carry the capture
sequence. All three bead primers extend in the same RT, so mRNA and feature molecules from
one cell share a GEM but **not the same barcode sequence**: the upstream page states the
feature-primer barcodes differ from the dT-primer barcodes in the same bead and must be
translated with the Cell Ranger file `translation/3M-february-2018.txt.gz` (v3) or
`3M-3pgex-may-2023.txt.gz` (v4). 🟡 (upstream only)

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | A **second, non-poly(dT) capture primer on the bead** with a different Read 1 handle (Nextera Read 1N), so feature libraries are separable by PCR from gene-expression cDNA | CITE-seq/Perturb-CITE-seq instead put a poly(A) tail on the ADT and share the dT primer |
| 2 | **Direct sgRNA capture** at the 3' bead via a capture-sequence insert in the scaffold, no polyadenylated reporter needed | ECCITE-seq does direct sgRNA capture on the 5' bead with a scaffold-specific RT primer; CROP-seq/Perturb-seq read a Pol II reporter |
| 3 | Antibody tags and sgRNAs **share** Capture Sequence 1; the read-2 side differs (TruSeq Read 2 in the ADT vs. TSO for the sgRNA) | — |

Builds on [reverse transcription](../ref/concepts/reverse-transcription.md) and, for the
sgRNA arm and the mRNA arm, [template switching](../ref/concepts/template-switching.md).
The antibody arm does **not** need template switching (§3, step 4).

## 2. Oligos

### Confirmed in the 10x 3' v3 user guide (CG000183, appendix "Oligonucleotide Sequences") — 🟢

```
Gel Bead, dT       5'-CTACACGACGCTCTTCCGATCT-NNNNNNNNNNNNNNNN-NNNNNNNNNNNN-TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3'
Gel Bead, CS1      5'-GTCAGATGTGTATAAGAGACAG-NNNNNNNNNNNNNNNN-NNNNNNNNNNNN-TTGCTAGGACCGGCCTTAAAGC-3'
                   ("Primer not used for 3' Gene Expression")
TSO (PN-3000228)   5'-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3'
cDNA Primers (PN-2000089)  Fwd 5'-CTACACGACGCTCTTCCGATCT-3'   Rev 5'-AAGCAGTGGTATCAACGCAGAG-3'
Chromium i7 Sample Index (PN-220103)  5'-CAAGCAGAAGACGGCATACGAGAT-NNNNNNNN GTGACTGGAGTTCAGACGTGT-3'
```

The guide's intro figure calls the bead dT "Poly(dT)VN", the appendix writes T30 without VN;
the v4 guide (CG000731) writes `TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN`. 🟢 both; the v3
inconsistency is the guide's own.

### From the upstream page only — 🟡 (secondary source; 10x FB guide not fetched)

```
Gel Bead, CS2      |--5'- GTCAGATGTGTATAAGAGACAG[16-bp cell barcode][12-bp UMI]CCTTAGCCGCTAATAGGTGAGC -3'

Capture Sequence 1 rc:  5'- GCTTTAAGGCCGGTCCTAGCAA -3'
Capture Sequence 2 rc:  5'- GCTCACCTATTAGCGGCTAAGG -3'

Antibody DNA oligo:  5'- GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[random 9-mer][15-bp antibody barcodes][random 9-mer]GCTTTAAGGCCGGTCCTAGCAA -3'

Feature cDNA Primers 1 (PN-2000096, CRISPR): mRNA pair = PN-2000089 pair above, plus
   sgRNA Fwd   5'- GCAGCGTCAGATGTGTATAAGAGACAG -3'
   sgRNA Rev   5'- AAGCAGTGGTATCAACGCAGAG -3'
Feature cDNA Primers 2 (PN-2000097, surface protein): mRNA pair, plus
   ADT Fwd     5'- GCAGCGTCAGATGTGTATAAGAGACAG -3'
   ADT Rev     5'- GTGACTGGAGTTCAGACGT -3'
Feature SI Primers 1 (PN-2000098, CRISPR only):
   Fwd 5'- AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -3'
   Rev 5'- GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGCAGTGGTATCAACGCAGAG -3'
Feature SI Primers 2 (PN-2000099):
   5'- AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -3'
```

sgRNA scaffolds as the upstream page writes them (RNA, then the capture sequence in DNA
letters as on the page). GeCKOv2 / mouse v1 / CROP-seq backbone:

```
5'- GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU -3'
```

(🟢 as DNA in Sanjana 2014, `006726_v1.full.pdf` and `nmeth.3047_*MOESM647_ESM.pdf`.)
Perturb-seq backbone (Adamson 2016, Addgene #85967; paper not in sources, 🟡):

```
5'- GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU -3'
```

and the four capture-sequence placements (upstream, following the 10x guide-RNA note):

```
CS1 in hairpin  5'- [Protospacer]GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGGCCGCTTTAAGGCCGGTCCTAGCAAGGCCAAGUGGCACCGAGUCGGUGCUUUUUUU -3'
CS2 in hairpin  5'- [Protospacer]GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGGCCGCTCACCTATTAGCGGCTAAGGGGCCAAGUGGCACCGAGUCGGUGCUUUUUUU -3'
CS1 at 3' end   5'- [Protospacer]GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCGCTTTAAGGCCGGTCCTAGCAAUUUUUUU -3'
CS2 at 3' end   5'- [Protospacer]GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCGCTCACCTATTAGCGGCTAAGGUUUUUUU -3'
```

### Other 10x-3'-v3 readouts in the sources (not Feature Barcoding) — 🟢

Perturb-CITE-seq (Frangieh 2021, Extended Data / Supplementary Table 8, `MOESM3_ESM.xlsx`)
ran 10x 3' v3 *without* the capture sequences: CITE-seq ADTs (poly(A)-tailed, captured by
the dT primer) and a CROP-seq "dial-out" PCR from the whole-transcriptome amplification:

```
CROPDialOut_R1       CTACACGACGCTCTTCCGATCT
CropDialOut_U6_F     GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTGTGGAAAGGACGAAACACC
CropDialOut_P5_R1    AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC
Shari oligos 96 well CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

## 3. How the oligos interlock — 🟡 (computed with `lib/`)

Handles:

- Bead dT primer handle `CTACACGACGCTCTTCCGATCT` = **`illumina.TRUSEQ_READ1[-22:]`** (partial Read 1).
- Bead CS1/CS2 primer handle `GTCAGATGTGTATAAGAGACAG` = **`nextera.READ1_PRIMER[-22:]`**
  = the last 3 nt of `nextera.S5` + the full 19-nt **`nextera.ME`**. "Read 1N" is the
  Nextera Read 1 primer, so the feature reads use a *different* Read 1 primer from the
  gene-expression reads, on the same flow cell.
- Feature cDNA Fwd `GCAGCGTCAGATGTGTATAAGAGACAG` (27 nt) = `nextera.READ1_PRIMER[-27:]`:
  5 nt longer than the bead handle, so it extends the handle on the first PCR cycles.
- Feature SI Fwd (PN-2000098/99) = **`illumina.P5` + `nextera.READ1_PRIMER`** exactly (62 nt).
- Antibody oligo 5' end = **`illumina.TRUSEQ_READ2`** exactly (34 nt); ADT cDNA Rev
  `GTGACTGGAGTTCAGACGT` = `TRUSEQ_READ2[:19]`.
- TSO DNA part `AAGCAGTGGTATCAACGCAGAGTACAT` starts with **`rt.SMART_HANDLE`**; the
  cDNA Rev / "partial TSO" primer = `SMART_HANDLE[:22]`. Copied into cDNA, the TSO end
  reads `CCCATGTACTCTGCGTTGATACCACTGCTT` (30 nt) on the top strand.
- Feature SI Primers 1 Rev = **`illumina.TRUSEQ_READ2` + partial TSO** (56 nt): it
  converts the sgRNA amplicon's TSO end into a TruSeq Read 2 end, so that **both** feature
  library types end in Read 2 and take the same Chromium i7 primer.
- Chromium i7 primer = **`illumina.P7`** + 8-nt index + `TRUSEQ_READ2[:21]`.
- The 10x GEX sample-index forward primer (and Frangieh's `CropDialOut_P5_R1`) =
  `illumina.P5` + `TRUSEQ_READ1[4:24]` — P5 ends in `ACAC`, which is also the first 4 nt of
  TruSeq Read 1, so the two overlap.
- `CropDialOut_U6_F` = `TRUSEQ_READ2` + `TGTGGAAAGGACGAAACACC`, which ends in
  `crispr.U6_3PRIME` (the U6 promoter just before the spacer).

Capture sequences and sgRNA designs:

- CS1rc = revcomp(CS1) and CS2rc = revcomp(CS2); both 22 nt. The ADT's 3' end is CS1rc,
  so the **antibody oligo hybridises to the CS1 bead primer** and the bead primer extends
  across it (DNA template; MMLV RT has DNA-dependent polymerase activity).
- "3'-end" sgRNA designs = Perturb-seq scaffold minus its 7 terminal U + CSrc + 7 U
  (verified for CS1 and CS2). The capture sequence sits immediately before the Pol III terminator.
- "Hairpin" designs = scaffold up to `…CAACUUG`, then `GCC` + CSrc + `GGCC` replacing the
  `AAA` of the loop, then the rest of the scaffold (verified for CS1 and CS2). The CSrc is
  displayed in the loop of a hairpin closed by GGCC/GGCC.
- The upstream Perturb-seq scaffold differs from `crispr.SCAFFOLD_FE` (Chen 2013 F+E) at one
  position: `GCUAAGCUGG` vs `GCTATGCTGG` (index 12, A vs T). Which is the Addgene #85967
  sequence was not checked against a primary source. 🔴
- The scaffold part of the sgRNA cDNA (top strand, between CS1 and the protospacer) is
  exactly revcomp of the scaffold up to `…GGUGC` (86 nt) — i.e. in the 3'-end design the bead
  primer copies the whole scaffold before reaching the spacer.

Lengths (computed from the segments, 20-nt protospacer):

| Product | Computed | Upstream says |
|---|---|---|
| ADT cDNA after Feature cDNA PCR | **144 bp** | 129 bp — **disagrees** (exactly 15 = the antibody barcode short) |
| sgRNA cDNA after Feature cDNA PCR | 213 bp | ~210 bp — agrees |
| Final ADT library | 211 bp | 211 bp — agrees |
| Final CRISPR library | 314 bp (311 with a 17-nt spacer) | ~314 bp — agrees |

## 4. Step by step

Reaction conditions for the feature arms are in the 10x FB guide, which was not fetched;
order and primers are 🟡 from the upstream page, checked against the GEX guide where the
steps coincide.

1. **Stain cells** with TotalSeq-B-type antibodies (surface protein) and/or use cells
   expressing capture-sequence sgRNAs (CRISPR). 🟡
2. **GEM generation** on Chip B (v3) / GEM-X chip (v4); the bead dissolves, cell lyses. 🟢
3. **RT in the GEM** (MMLV-type RT, 53 °C 45 min, 85 °C 5 min in the GEX guide 🟢 for v3):
   - mRNA: primed by the dT primer → template switch onto the TSO (standard 10x 3').
   - Antibody oligo: its 3' CS1rc anneals to the bead CS1 primer, which extends over the
     DNA oligo to its 5' end (copying barcode, random 9-mers and the TruSeq Read 2
     handle). The antibody oligo's 3' end is assumed blocked. 🟡 (upstream's assumption)
   - sgRNA: the CSrc anneals to the CS1 or CS2 primer; the primer copies the scaffold and
     the protospacer to the 5' end of the RNA, adds non-templated C's, and template
     switches onto the TSO. 🟡
4. **cDNA amplification** with Feature cDNA Primers (mixture with the mRNA pair):
   - surface protein (PN-2000097): ADT Fwd (Read 1N side) + ADT Rev (TruSeq Read 2 side).
     The antibody arm is amplified **between the bead handle and the antibody oligo's own
     Read 2 handle** — the TSO end, if made, is outside the amplicon. 🟡 (computed from primer positions)
   - CRISPR (PN-2000096): sgRNA Fwd (Read 1N side) + partial TSO.
5. **Size split**: mRNA cDNA mainly > 1 kb, feature products ~144 / ~213 bp (sizes as
   computed in §3); upstream only says "size selection". That the small feature fraction is
   recovered from an SPRI supernatant is an inference, not in the sources. 🟡
6. **CRISPR only — Feature PCR** with Feature SI Primers 1: adds P5 + full Read 1N and
   swaps the TSO end for a TruSeq Read 2 end. 🟡
7. **Sample-index PCR** with Feature SI Primers 2 (P5 + Read 1N) and the Chromium i7
   Sample Index primer (P7 + i7 + Read 2[:21]). 🟡 The gene-expression library is made
   separately by fragmentation, A-tailing and Read 2 adaptor ligation (the 3' GEX
   chemistry; not repeated here).

## 5. Final libraries — 🟡 (assembled from the oligos above)

Gene expression (v3, single index; for reference):

```
5'- P5 · TruSeq Read 1 · <cbc 16> · <umi 12> · dT30(VN) · <cDNA> · TruSeq Read 2' · <i7 8> · P7' -3'
```

Surface protein (antibody-derived), 211 bp:

```
5'- P5 · Nextera Read 1 (33) · <cbc 16> · <umi 12> · CS1 (22) · <N9> · <antibody bc 15> · <N9> · TruSeq Read 2' (34) · <i7 8> · P7' -3'
```

CRISPR (sgRNA, 3'-end CS design), ~314 bp:

```
5'- P5 · Nextera Read 1 (33) · <cbc 16> · <umi 12> · CS1 or CS2 (22) · scaffold' (86) · <protospacer' 17–20> · TSO' (30) · TruSeq Read 2' (34) · <i7 8> · P7' -3'
```

Here `'` means "as its reverse complement on this strand". In the hairpin design the scaffold
segment would be only the part 3' of the insertion. 🔴 (not drawn upstream).

## 6. Sequencing

🟢 from the GEX guides: Read 1 28 cycles (16-nt cell barcode + 12-nt UMI); v3: i7 8,
i5 0, Read 2 91; v4 (dual index): i7 10, i5 10, Read 2 90.

🟡 (upstream): the Illumina Read 1 primer mix contains both TruSeq and Nextera Read 1
primers, so mixed GEX + feature pools need no custom primer. Index read with
`illumina.INDEX1_PRIMER`. Read 2 = TruSeq Read 2 for all three library types:

- ADT: Read 2 bases 10–24 are the 15-nt antibody barcode (after the 9-nt random block). 🟡
- sgRNA: Read 2 starts with 30 nt of TSO complement, then the protospacer at ~31–50,
  then scaffold. 🟡 (consistent with the computed 30-nt TSO' segment)

## 7. Agreements and disagreements with the upstream page

- **Agree** 🟢: bead dT and CS1 primers, TSO, cDNA primers, i7 sample-index primer, read 1
  28 cycles — identical in CG000183.
- **Disagree**: upstream says the ADT cDNA is **129 bp**; its own segments sum to **144 bp**
  (the final 211-bp library is consistent with 144, not 129).
- **Disagree / incomplete**: upstream says v4 library structure is "exactly the same" as
  v3. CG000731 says the v4 bead carries **only Capture Sequence 1** (no CS2), and v4
  libraries are **dual-indexed** (10 + 10). The v4 feature libraries presumably follow the
  same pattern with a dual-index i5 primer; not in the sources. 🔴
- **Typo upstream**: the bottom strand of the final CRISPR library (and the sequencing panels
  that reuse it) ends "`-3'`" where it must be "`-5'`".
- **Not checkable here**: CS2 sequence, the antibody oligo layout, all PN-200009x primer
  sequences — need the 10x FB user guide appendix. 🔴

## 8. Open questions

- 🔴 Primary sequences of Feature cDNA Primers 1/2, Feature SI Primers 1/2 and Capture
  Sequence 2 (10x FB user guides CG000185 / CG000184 appendices, and v4 equivalents).
- 🔴 Whether the antibody oligo's 3' end is blocked (upstream assumes so) and whether the
  RT also template-switches off the ADT (irrelevant to the amplicon, but it is the
  difference between "TotalSeq-B" and "TotalSeq-A" behaviour).
- 🔴 v4: feature-library index primer (dual-index i5/i7) and translation file details.
- 🔴 Which exact scaffold Addgene #85967 carries (A/T at scaffold index 12, §3).
- 🟡 Reaction conditions for the Feature cDNA / Feature PCR / SI PCR steps.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` had already fetched 90+ files, but for a kit protocol the useful
ones were the upstream page and the two 10x **gene-expression** guides; the Feature
Barcoding guides themselves were neither fetched nor listed as manual downloads. Most of
the bibliography rows are CRISPR-vector background papers, whose large xlsx tables (one
36 MB) made `tools/scrape_primers.py` on the whole directory exceed two minutes;
`--find` checks were fast (≈ 18 s). All §2/§3 relationships were computed with `lib/`.
