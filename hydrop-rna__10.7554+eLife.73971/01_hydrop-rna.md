# HyDrop-RNA — Drop-seq chemistry on dissolvable, split-pool barcoded hydrogel beads

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/`; the model
> is `tools/hydrop_rna.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length, kit or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence. Every HyDrop-specific base below comes from the
> authors' own eLife supplementary files — Supplementary file 4 ("Reagents and
> oligonucleotide list", the ordered oligos) and Supplementary file 3 ("Molecular
> sequence description of HyDrop-RNA") — and every generic Illumina base from
> `lib/illumina.py`, i.e. from *Illumina Adapter Sequences* #1000000002694. A filename
> prefix is not provenance: the `scg_elife-73971-supp{1,2,3}-v4.docx` and
> `upstream_data/elife-73971-supp{1,3}-v4.docx` files on disk are the publisher's own
> supplements, byte-identical (md5) to the copies fetched from `cdn.elifesciences.org`,
> and are cited here as such.

**HyDrop** — De Rop FV, Ismail JN, Bravo González-Blas C, Hulselmans GJ, Flerin CC,
Janssens J, Theunis K, Christiaens VM, Wouters J, Marcassa G, de Wit J, Poovathingal S,
Aerts S. "Hydrop enables droplet-based single-cell ATAC-seq and single-cell RNA-seq using
dissolvable hydrogel beads." *eLife* 2022;11:e73971.
doi:[10.7554/eLife.73971](https://doi.org/10.7554/eLife.73971) (PMC8993220, PMID 35195064).

Other papers / documents for this method (from `catalogue/scg_lib_structs.tsv`):

- Protocol: "HyDrop-RNA v1.0", protocols.io, doi:[10.17504/protocols.io.b4xwqxpe](https://doi.org/10.17504/protocols.io.b4xwqxpe) — not fetched (see below).
- Sister assay from the same paper: **HyDrop-ATAC** (`hydrop-atac__10.7554+eLife.73971`;
  protocols.io b4xvqxn6). Same beads up to the last barcoding round, same `HYi7` primers.
- `catalogue/scg_lib_structs.tsv` carries a `cited` row for Plessy et al.,
  nanoCAGE/CAGEscan, *Nat Methods* 2010, doi:10.1038/nmeth.1470 ("semi-suppressive PCR").
  That association was scraped from the drawn page's own citation list, and the eLife
  paper does not cite Plessy anywhere 🟡 — it is background reading for suppression PCR,
  not a HyDrop source, and nothing in this note rests on it.

Sources read (in `$CHEM_DATA/sources/hydrop-rna__10.7554+eLife.73971/`, never committed):

| File | What | Used for |
|---|---|---|
| `eLife.73971_PMC8993220.html.txt` | full text (PMC) | methods: bead making, RT, cDNA PCR, library prep, sequencing; optimisation trials |
| `elife_cdn/elife-73971-supp1-v3.docx.txt` | Supplementary file 1, "Molecular sequence description of HyDrop bead barcoding" | bead oligo after each extension, with one example barcode |
| `elife_cdn/elife-73971-supp3-v3.docx.txt` | Supplementary file 3, "Molecular sequence description of HyDrop-RNA" | every step of the RNA chemistry as sequences; read primers |
| `elife_cdn/elife-73971-supp4-v3.xlsx.txt` | Supplementary file 4, "Reagents and oligonucleotide list" | **every oligo**, incl. the 3 × 96 barcoding plates |
| `elife_cdn/elife-73971-supp2-v3.docx.txt` | Supplementary file 2 (HyDrop-ATAC sequences) | only to confirm the shared `HYi7` design |
| `elife_cdn/elife-73971-supp5-v3.xlsx.txt` | Supplementary file 5, cost sheets | enzyme list (Maxima H–, NEB Ultra II FS, KAPA) |
| `elife_cdn/elife-73971-fig1-figsupp{2,3,4}-data1-v3.docx.txt` | bead QC source data | not chemistry |
| `eLife.73971_elife-73971-supp4-v5.xlsx.txt` | Supplementary file 4 of the **version of record** (v5), direct from `cdn.elifesciences.org` | the ordered-oligo table this note and the model are built on; supersedes the v3 copy (see §6) |
| `scg_elife-73971-supp3-v4.docx.txt`, `upstream_data/elife-73971-supp{1,3}-v4.docx.txt` | the publisher's Supplementary files 1 and 3 (v4), mirrored in the scg_lib_structs `data/` tree — **genuine eLife originals**, md5-identical to the eLife CDN copies | cited as the authors' own documents; the `scg_`/`upstream_data` prefix is a download route, not provenance |
| `scg_20210712_supp_methods_table_hydrop_oligonucleotide_list.xlsx.txt` | the **authors'** protocols.io oligonucleotide list (dated 2021-07-12), mirrored in the same tree | oligo content identical to Supplementary file 4 v3 🟡 (diffed); carries the same v3 name swap |
| `protocolsio_hydrop-rna-v1-0-dm6gpwqjjlzp.pdf.txt` (v5, Dec 2023), `protocolsio_hydrop-bead-generation-amp-pcr-barcoding-v1-0-n92ld9np9g5b.pdf.txt` | the authors' protocols.io protocols (PDF export), now fetched | bench conditions only — DTT/GTP/PEG, Exo I, 13-cycle PCRs, 20 °C ligation; they print **no oligo sequence** 🟢 (checked) |
| `upstream_HyDrop_RNA.html.txt` | scg_lib_structs page | **inspiration only, not a source** — nothing in this note or the model is cited from, or verified against, it |

Not fetched:

| What | URL | Why it would matter |
|---|---|---|
| LNA-TSO and dN-SMRT primer sequences (optimisation trials) | — | named in the methods, absent from Supplementary file 4 and from both protocols.io protocols 🔴 — they are published nowhere we can reach, and the model does not draw them |

(The PMC supplement links were behind the download gate, and Europe PMC's
supplementary zip was cut off at 50 MB; the same files came from
`cdn.elifesciences.org` — see §8.)

---

## 1. What it is

A droplet scRNA-seq method built from **inDrop's bead making** and **Drop-seq's RNA
chemistry** 🟡. Cells and RT mix are co-encapsulated with polyacrylamide hydrogel beads that
carry the barcoded oligo-dT primer; DTT in the RT mix **dissolves** the beads (disulfide
cross-linker BAC; the oligo is also anchored through an internal disulfide,
`/iThioMC6-D/`), so the primer is free in the droplet rather than tethered 🟢 (bead monomer
mix and DTT in the RT mix are in the methods; "dissolvable" is the paper's word) / 🟡 (that
DTT is what dissolves them in the droplet is inferred from the chemistry).

| | New thing here | Builds on |
|---|---|---|
| 1 | **Dissolvable** hydrogel beads (BAC cross-linker + disulfide-linked acrydite primer), barcoded by **3 × 96 split-pool primer extension** → 884,736 barcodes 🟢 | inDrop bead barcoding (Klein 2015; Zilionis 2017), acrydite hydrogel beads (Ren 2021 / Wang 2020) |
| 2 | One unbarcoded bead → **RNA or ATAC** beads by choosing only the last barcoding plate 🟢 | — |
| 3 | RNA side copied from Drop-seq: same 3' handle on the bead, same TSO, same single-primer cDNA PCR 🟡 *(computed, see §2)* | Drop-seq ([note](../drop-seq__10.1016+j.cell.2015.05.002/01_drop-seq-seq-well.md)); [template switching](../ref/concepts/template-switching.md) |
| 4 | Library made by **NEB Ultra II FS** (enzymatic fragmentation + dA-tail + NEB hairpin adapter) instead of Nextera, and the **cell barcode read as Read 2** with a custom primer; the sample index sits on **both** P5 and P7 🟢 | — |
| 5 | **Exo I** after droplet breaking, to remove unused bead primers before the bulk PCR; **GTP + PEG** in the RT 🟢 | GTP/PEG from Smart-seq3 🟢 (paper); Exo I clean-up as in Drop-seq 🟡 |

## 2. Oligos

🟢 Verbatim from Supplementary file 4 (`supp4-v3.xlsx`, sheet "Protocol_primers" and the
three barcoding-plate sheets). IDT notation as given: `/5Acryd/` = 5' acrydite,
`/iThioMC6-D/` = internal thiol-modifier C6 S-S (disulfide), `rG` = ribo-G, `*` =
phosphorothioate, `/56-FAM/` = 5' FAM.

```
Acrydite_primer        /5Acryd//iThioMC6-D/TTTTTTTTAATACGACTCACTATAGGGAAGCAGTGGTATCAACGCAGAGTAC     HPLC
20200130_plate-1-96    GCAGTAGCTG <BC1> GTACTCTGCG           96 oligos, 30 nt; e.g. GCAGTAGCTGTGTAGCAAGTGTACTCTGCG
20200130_plate-2-96    AGGGTACTCG <BC2> GCAGTAGCTG           96 oligos, 30 nt; e.g. AGGGTACTCGTTAGTTGGACGCAGTAGCTG
20200130_plate-3-96-RNAseq
                       AAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN <BC3> AGGGTACTCG
                                                             96 oligos, 58 nt; e.g. AAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNTGACCGTACTAGGGTACTCG

TSO                    AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG                          HPLC
TSO-P                  AAGCAGTGGTATCAACGCAGAGT                                    "primer for post-RT ISPCR"

HYi7_N_<i7>            CAAGCAGAAGACGGCATACGAGAT <i7> CTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C
                       N = 1..16, e.g. HYi7_1_CGCTCAGTTC = CAAGCAGAAGACGGCATACGAGATCGCTCAGTTCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C
HYi5_TruSeq_N_<i5>     AATGATACGGCGACCACCGAGATCTACAC <i5> ACACTCTTTCCCTACACGACGCT
                       N = 1..5, e.g. HYi5_TruSeq_1_TCGTGGAGCG = AATGATACGGCGACCACCGAGATCTACACTCGTGGAGCGACACTCTTTCCCTACACGACGCT

HyDrop_CustSeq_R2      CTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC     HPLC   "custom read 2 read primer for hydrop-rna"
HyDrop_CustSeq_Short   GTACTCTGCGTTGATACCACTGCTTCCGCGGACAG     HPLC   "custom i7 read primer for hydrop-rna"
```

(These are the names and sequences of the **version of record**, `supp4-v5`. The earlier
`supp4-v3` copy, and the authors' protocols.io oligo list, have the two sequences under
each other's names — see the box below.) Bead QC probes (not used in the library) 🟢:

```
Anti-Acrydite_FAM   /56-FAM/TTTTTGTACTCTGCGTTGATACCAC
Anti-RNA_FAM        /56-FAM/AAAAAAAAAAAAAAAAAAAA
Anti-BC1_FAM        /56-FAM/TTTTTCTATCCGTCAGTAC
Anti-BC2_FAM        /56-FAM/TTTTTACACGTTGTGGCAG
Anti-BC3_FAM        /56-FAM/TTTTTCTCCTATCATAGGG
```

Kit oligos used but not listed in the paper: the **NEBNext Adapter for Illumina** (hairpin,
dU in the loop). The paper names the kit ("NEBNext Ultra II FS", "NEBNext Adapter for
Illumina") 🟢, and the bases come from NEB's own adapter specification by way of
`illumina.NEBNEXT_HAIRPIN` 🟢 (bases) / 🟡 (that this is the adapter version the authors
received) — paper → named vendor kit → vendor document, with nothing taken from the drawn
page. The sequencing run's **TruSeq Read 1** and **Index 2** primers
are drawn in Supplementary file 3: `ACACTCTTTCCCTACACGACGCTCTTCCGATCT` and
`AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` 🟢.

### How the oligos interlock — 🟡 (computed)

**The bead primer.** `Acrydite_primer` (52 nt) = `TTTTTTT` + **T7 promoter**
`TAATACGACTCACTATAGGG` + **`rt.SMART_HANDLE`** (`AAGCAGTGGTATCAACGCAGAGT`) + `AC`. The
3'-terminal 25 nt `AAGCAGTGGTATCAACGCAGAGTAC` are base for base the 5' part of the
**Drop-seq bead oligo** (and of Drop-seq's custom Read 1 primer, see below). The T7 promoter
is never used for in-vitro transcription here — vestigial from the inDrop bead design. 🟡
(the inDrop lineage is stated by the paper; "vestigial" is inferred: no IVT step exists)

**Split-pool extension.** Each barcoding oligo is a *template* that anneals by its 3'-terminal
10 nt to the current 3' end of the bead primer, which is then extended across it
(KAPA HiFi; 5 cycles of 98/38/72 °C + one final cycle) 🟢 conditions:

- plate 1 3' end `GTACTCTGCG` = revcomp of the bead's 3' `CGCAGAGTAC` → bead gains
  `<BC1>` + `CAGCTACTGC` (= revcomp of plate-1's 5' `GCAGTAGCTG`).
- plate 2 3' end `GCAGTAGCTG` = revcomp of `CAGCTACTGC` → bead gains `<BC2>` + `CGAGTACCCT`.
- plate 3 3' end `AGGGTACTCG` = revcomp of `CGAGTACCCT` → bead gains `<BC3>` + `N8` + `T30`.
- All 96 oligos in each plate share these 10-nt ends; BC1, BC2, BC3 are 10 nt each.
- Between rounds: EDTA stop, then NaOH strips the template strand 🟢.

**Final RNA bead oligo**, 5'→3' (supp4 design): `acrydite · S–S · T7 · SMART handle · AC ·
BC1(10) · CAGCTACTGC · BC2(10) · CGAGTACCCT · BC3(10) · UMI(8) · T30` = **140 nt** of DNA. 🟡

**TSO / TSO-P / HYi7 — one handle, three different 3' ends.**

- `TSO-P` **is** `rt.SMART_HANDLE` (23 nt). It anneals to the complement of the handle at
  **both** cDNA ends (bead end and TSO end), so the cDNA PCR is single-primer
  ("IS-PCR" / suppression PCR), exactly as in Drop-seq.
- `TSO` = handle + **`GAAT`** + `rGrGrG` (30 nt) — the Drop-seq TSO.
- `HYi7` = **`illumina.P7`** + 10-nt i7 + `CTGTCCGCGG` + handle + **`AC`** (69 nt), with
  phosphorothioates on the last two linkages. Its 3' end matches the **bead** end
  (`…GCAGAGTAC`) but not the **TSO** end (`…GCAGAGTGA`): the two 3'-terminal bases mismatch.
  So in the library PCR only fragments carrying the bead end (barcode, UMI) amplify; 5'-end
  fragments that kept the TSO handle do not. The `*A*C` presumably stops KAPA HiFi's
  3'→5' proofreading from trimming the mismatch away and rescuing those fragments. 🟡
  (the mismatch is computed from the two ordered oligos; the role of the phosphorothioates
  is inferred, not stated anywhere)
- The 3' 35 nt of `HYi7` (`CTGTCCGCGG` + handle + `AC`) equal the **Drop-seq custom Read 1
  primer** minus its first two bases `GC` (compared with the Drop-seq note). HyDrop moved
  that spacer from the P5 side (Drop-seq's P5-SMART hybrid primer) to the P7 side and put an
  index in front of it. 🟡
- The i7 is written **as named**: `HYi7_1_CGCTCAGTTC` contains `CGCTCAGTTC` — true for all 16.
- `HYi7` is shared with HyDrop-ATAC (same bead handle) 🟢 (table description; Supp. file 2).

**HYi5_TruSeq.** = **`illumina.P5`** + 10-nt i5 + `illumina.TRUSEQ_READ1[:23]`
(`ACACTCTTTCCCTACACGACGCT`), 62 nt — the Read 1 arm stops **10 nt short** of
`CTTCCGATCT`, which the NEB hairpin has already installed in the fragment. The model
(`hydrop_rna.HYI5_READ1_ARM`) now draws that truncation; an earlier revision gave the
primer the full 33-nt arm, which `tools/selftest_source_oligos.py` flagged as absent from
the source text. Its 3' end pairs inside the read-1 arm of the NEB
adapter. The i5 is written **as named** for all five, and the five i5 indices are the first
five **`HYi5` (ATAC) indices** — the ATAC primers have the Nextera s5 tail
(`TCGTCGGCAGCGTCAGATGTG` = `nextera.ADAPTOR_S5[:21]`) instead.

**The custom read primers.** `HyDrop_CustSeq_R2` = the 3' 35 nt of `HYi7`;
`HyDrop_CustSeq_Short` = its reverse complement. 🟡

> **Resolved disagreement — read-primer names across supplement versions.** In the final
> library the barcode/UMI region follows `CTGTCCGCGG·handle·AC` on the P7-side strand. A
> primer with **that** sequence (`CTGTCCGCGG…GCAGAGTAC`) extends into BC1 → it is the
> **Read 2** (barcode) primer; its reverse complement (`GTACTCTGCG…CGCGGACAG`) extends the
> other way into the i7 → it is the **Index 1** primer. Two independent grounds settle
> this, neither of them the drawn page:
>
> 1. **The version of record.** `supp4-v5` names them this way 🟢
>    (`HyDrop_CustSeq_R2` = `CTGTCC…AGTAC`, "custom read 2 read primer";
>    `HyDrop_CustSeq_Short` = `GTACTC…GACAG`, "custom i7 read primer"), and
>    Supplementary file 3's drawing labels the same two arrows the same way 🟢. The
>    *earlier* `supp4-v3` copy and the authors' protocols.io oligo list carry the two
>    sequences under each other's names — a defect corrected between eLife versions, not
>    a live conflict.
> 2. **Computation.** `tools/hydrop_rna.py` declares both primers by reference and
>    `seqprimers.verify` locates them on the modelled library; `final_library()` raises if
>    either fails to land. The assignment is therefore enforced by the model rather than
>    asserted in prose. 🟡 *(computed)*

### Bead QC probes — 🟡 (computed)

`Anti-Acrydite_FAM` = `TTTTT` + revcomp of the bead primer's 3'-terminal 20 nt. Each
`Anti-BCn_FAM` = `TTTTT` + 14 nt taken from **one** plate oligo (BC1: plate-1 oligo #3;
BC2: plate-2 #28; BC3: plate-3 #19), i.e. complementary to one of 96 barcodes — the "1 in 96"
sub-barcode test of the methods (counts in the Fig. 1 figure-supplement source data). `Anti-RNA_FAM` (A20) checks the poly(T).

## 3. Step by step

Volumes and conditions 🟢 from "HyDrop-RNA single-cell library preparation" unless marked.

**Beads** (once, in bulk):

1. **Polymerise**: 6 % acrylamide, 0.55 % BAC, 12 µM acrydite primer, 0.6 % APS, in 10 % TBSET,
   emulsified into ~50 µm droplets, 65 °C 14 h; break, wash (PFO, hexane/SPAN-80).
2. **Barcode, × 3**: split into 96 wells, 2.5 µL 100 µM plate oligo + KAPA HiFi HotStart;
   95 °C 3 min, 5 × (98 °C 20 s, 38 °C 4 min, 72 °C 2 min), 98 °C 1 min, 38 °C 10 min,
   72 °C 4 min. Stop with EDTA, pool, NaOH (150 mM) to remove the template strand,
   neutralise. Plate 3 = `plate-3-96-RNAseq` for RNA beads.
3. **Store**: glycerol Bead Freezing Buffer, –80 °C.

**In the droplet:**

4. **Co-encapsulate** ~3,800 cells (for 2,000 recovered) in 85 µL RT mix — 1× Maxima RT
   buffer, 0.9 mM dNTPs, **25 mM DTT**, **1.3 mM GTP**, 15 % OptiPrep, 1.3 U/µL RNase
   inhibitor, 15 U/µL **Maxima H– RT**, **12.5 µM TSO**, **4.4 % PEG-8000** — with 35 µL
   beads, on the Onyx platform (Droplet Genomics). DTT dissolves the bead and frees the
   primer 🟡.
5. **RT + template switching**: 42 °C 90 min; 11 × (50 °C 2 min, 42 °C 2 min); 85 °C 5 min.
   The bead primer's T30 primes on poly(A); at the 5' end of the mRNA MMLV adds `CCC`, the
   TSO's `rGrGrG` pairs and the RT copies the TSO, giving first strand
   `bead oligo · cDNA' · CCC · (TSO handle)'` 🟢 (Supp. file 3 draws this) — see
   [template switching](../ref/concepts/template-switching.md),
   [reverse transcription](../ref/concepts/reverse-transcription.md).

**Bulk:**

6. **Break** with 20 % PFO + GITC buffer + DTT; AMPure (99 µL), elute in EB-DTT-Tween.
7. **Exo I**, 37 °C 5 min / 80 °C 1 min — removes unused barcode primers so they cannot prime in the
   bulk IS-PCR (the paper's reason: barcode purity, Fig. 8—fig. suppl. 1); that Exo I acts
   because they are single-stranded 3' ends is the enzyme's known specificity 🟡.
   0.8× AMPure.
8. **cDNA PCR (IS-PCR)** with **TSO-P only**: KAPA HiFi, 10 µL of 10 µM TSO-P per 100 µL (1 µM final 🟡); 95 °C 3 min,
   13 × (98 °C 20 s, 63 °C 20 s, 72 °C 3 min), 72 °C 5 min. 0.6× AMPure. Product:
   `handle · AC · BC1 · L1 · BC2 · L2 · BC3 · UMI · T30 · cDNA · CCC · ATTCACTCTGCGTTGATACCACTGCTT`
   — the T7/T7-prefix of the bead oligo is lost because TSO-P primes inside it (Supp. file 3
   draws the shortened product 🟢).
9. **Fragment + dA-tail**: 80 ng cDNA, NEBNext Ultra II FS, 37 °C 10 min, 65 °C 30 min.
   0.8× AMPure.
10. **Ligate** NEBNext hairpin adapter (2.5 µL), 20 °C 15 min. 0.8× AMPure. The
    supplement's text says "NEB USER enzyme treatment" follows 🟢; the methods text lists
    no separate USER step, and in the Ultra II FS kit USER is added after ligation 🟡.
    Each fragment strand becomes `NEB read-1 arm · insert · A · NEB read-2 arm` 🟡.
11. **Index PCR**: KAPA HiFi, `HYi7` + `HYi5_TruSeq` (10 µL of 10 µM each per 100 µL, 1 µM final 🟡); 95 °C 3 min,
    13 × (98 °C 20 s, 64 °C 30 s, 72 °C 30 s), 72 °C 5 min; 0.8× AMPure. `HYi7` primes only
    on the bead-end fragments (step 8 reasoning); `HYi5_TruSeq` on the NEB read-1 arm.
    Internal fragments (no handle) and TSO-end fragments are not exponentially amplified;
    the NEB read-2 arm on the barcode side is **not** in the final molecule — `HYi7`
    primes inside the handle, 5' of it. 🟡 *(computed from the two ordered primers and the
    adapter)*

Optimisation trials (methods, 🟢): no-Exo I, an **LNA TSO**, GTP/PEG, and two
second-strand-synthesis alternatives (RNase H, then a random-primed "dN-SMRT" primer with
Klenow exo– or Bst 2.0, followed by the same ISPCR). TSO + ISPCR + GTP/PEG performed best
and is the protocol above. Sequences of the LNA TSO and dN-SMRT primer: 🔴.

## 4. Final library — 🟡 (assembled from the oligos above)

5'→3', top strand written from P5, the orientation these notes use throughout;
Supplementary file 3 writes the same molecule from P7:

```
5'- P5 · i5 · TruSeq Read 1 · <cDNA, sense> · A30 · UMI(8) · BC3(10) · AGGGTACTCG · BC2(10) · GCAGTAGCTG · BC1(10) · GTACTCTGCG TTGATACCACTGCTT CCGCGGACAG · i7' · P7' -3'
```

Segments: `illumina.P5` (29) · i5 (10) · `illumina.TRUSEQ_READ1` (33) · insert ·
A30 · N8 · BC3 · L2' · BC2 · L1' · BC1 · revcomp(`CTGTCCGCGG` + handle + `AC`) (35) · i7'
(10) · `illumina.P7_RC` (24). With `HYi5_TruSeq_1` and `HYi7_1`, i5 = `TCGTGGAGCG` and
i7' = `GAACTGAGCG`. The 3' part is computed from the two ordered library-PCR primers —
`GTACTCTGCGTTGATACCACTGCTTCCGCGGACAG <i7'> ATCTCGTATGCCGTCTTCTGCTTG` — and is the
reverse complement of `HYi7` read back across the amplicon, so the model derives it rather
than transcribing it. 🟡 *(computed)*

The cDNA insert is the 3' end of the transcript (the bead end is the only one carried
through), read in **sense** from Read 1 towards the poly(A).

## 5. Sequencing

🟢 NextSeq 2000: **Read 1 50 cycles** (3' cDNA, standard TruSeq Read 1), **Index 1
10 cycles** (sample index, custom i7 read primer), **Index 2 10 cycles** (sample index;
Supp. file 3: standard TruSeq Index 2 primer `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`),
**Read 2 58 cycles** (barcode + UMI, custom read 2 primer).

Read 2, 🟡 (computed): `BC1(10) · CAGCTACTGC · BC2(10) · CGAGTACCCT · BC3(10) · UMI(8)`
= 58 nt — the cycle count exactly covers an **8-nt** UMI. The linkers are trimmed
(mawk) before STARsolo `CB_UMI_Complex` with a 884,736-barcode whitelist, ≤1 mismatch 🟢.

Index reads 🟡: Index 1 primes on `GTACTCTGCG…GGACAG` (the table's "CustSeq_R2", the
drawings' "CustSeq_Short") and reads revcomp of the i7 as named (`HYi7_1` → `GAACTGAGCG`).
Index 2 (`illumina.INDEX2_PRIMER_RC` = revcomp of TruSeq Read 1) on the NextSeq 2000
reverse-complement workflow reads revcomp of the i5 as named (`HYi5_TruSeq_1` →
`CGCTCCACGA`).

## 6. Where the authors' own files disagree, and how each case is settled

The eLife package contains both **ordered oligos** (Supplementary file 4, an IDT order
sheet) and **schematic drawings** (Supplementary files 1–3), and the two do not always
match. The rule used here is: an ordered oligo is the molecule that existed; a drawing is
an illustration. Where that is not enough, the model's own computation decides. Nothing in
this table is settled by the drawn page.

| Point | Disagreement inside the authors' files | Settled by | Status |
|---|---|---|---|
| UMI / poly(T) width | Supp. 4 ordered plate-3-RNAseq oligos: `A30 + N8 + BC3 + AGGGTACTCG` 🟢 · Supp. 1 drawing: N9 + T25 ("136 bp") · Supp. 3 drawing: N10 + T25 | the ordered oligos, independently confirmed by the **58-cycle Read 2** (10+10+10+10+10+8 = 58) 🟡 *(computed)* | model draws N8 + T30 🟡 |
| Custom read-primer names | `supp4-v5`: R2 = `CTGTCC…`, Short = `GTACTC…` 🟢 · `supp4-v3` and the protocols.io list: swapped | the version of record, plus `seqprimers.verify` on the modelled library (§2 box) | resolved 🟡 |
| `HYi5_TruSeq` Read 1 arm | Supp. 4 prints a 23-nt arm; the full TruSeq arm is 33 nt | the ordered oligo — the remaining 10 nt come from the NEB hairpin, not from the primer | model corrected to 23 nt 🟢 |
| `HYi7` 3' end | Supp. 4: `…CAGAGT*A*C` 🟢 | modifications are recorded in §2 but not drawn as bases | 🟡 (phosphorothioates not modelled) |
| NEB adapter bases | paper names the kit only | NEB's own adapter specification via `illumina.NEBNEXT_HAIRPIN` | 🟡 (which adapter version) |
| TSO-end fragments not amplified | not stated anywhere | the 2-nt 3' mismatch between `HYi7`'s `…GCAGAGTAC` and the TSO end's `…GCAGAGTGA` 🟡 *(computed)* | inference, flagged as such |
| Supp. 1 example barcodes | `TTCCTGTGAG`, `TCGGACTTAT`, `GGCTGAATTA` are **not** in the plate-1/2/3 lists | the plate lists are the ordered material; the drawing's barcodes are illustrative | 🟡 |
| cDNA PCR | Supp. 3 draws one primer at both ends; Supp. 4 names one oligo, `TSO-P` | agree; the model uses `PCRGeometry.symmetric` so the single-species claim is checked | 🟢/🟡 |

## 7. Open questions

- 🔴 Why the **`AC`** after the handle on the bead (inherited from Drop-seq's
  `…CAGAGTAC`) and the **`GAAT`** in the TSO — the paper does not say; the computed effect is
  that they let `HYi7` tell the two cDNA ends apart.
- 🔴 Whether the USER step is a separate addition (Supp. 3) or omitted (methods text).
- 🔴 LNA-TSO and dN-SMRT primer sequences.
- 🟡 UMI length: 8 nt by the ordered oligos and the Read 2 length; the supplementary drawings
  show 9 and 10.
- 🔴 protocols.io v1.0 not read — conditions there may differ from the paper.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py "hydrop-rna__10.7554+eLife.73971"` matched **no protocol** (the full
slug did not match); `get_sources.py HyDrop-RNA` worked. It saved the PMC text and XML,
but: the Europe PMC supplementary zip was **truncated at exactly 50,000,000 bytes** (no
central directory; a few entries were recoverable from local headers); the
"transrepform1.docx" it saved was a **reCAPTCHA HTML page** with a .docx name; all other
supplements were `(manual)` behind the PMC gate. The v3 supplements came instead from
`cdn.elifesciences.org/articles/73971/elife-73971-<name>-v3.<ext>` (fetched by hand) and
the authors' protocols.io oligo list from the scg_lib_structs `data/` mirror; a later
fetch added the **v5** (version-of-record) supplements from the same eLife CDN, and §6
uses those. `scrape_primers.py` found every oligo
table (Supp. 4 first), but read `TSO-P` as a 5'-phosphorylated `TSO` (the `-P` of the name
taken as `5' P`), and joined a duplex's two strands into one hit in Supp. 1. Its `--find`
reported **0 locations** for sequences containing `N` and missed the Supp. 4 rows written
with `*` (phosphorothioate), finding only the unmodified copies elsewhere; those were
checked by hand with `grep -F`.
