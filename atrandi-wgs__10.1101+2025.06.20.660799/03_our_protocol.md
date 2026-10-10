# Our protocol (Bascet/Zorn preprint) — quick reference

**Gourlé H\*, Yakovenko I\*, Verma J\*, Dicken J\*, Albrecht F, Stenlund A, Sharaf C, Dwibedi C,
Mažutis L, Normark J, Sheng N, Strömberg N, Löfstedt T, Carroll LM<sup>#</sup>, Henriksson J<sup>#</sup>.
"Scalable single-cell metagenomic analysis with Bascet and Zorn." bioRxiv (2025/2026).**
doi:[10.1101/2025.06.20.660799](https://doi.org/10.1101/2025.06.20.660799) · CC-BY 4.0
v1 posted 21 Jun 2025; **v2 posted 28 Sep 2026** (the version documented here).
Full text: <https://www.biorxiv.org/content/10.1101/2025.06.20.660799v2.full>


> **Evidence marking.** 🟢 = verbatim from the source document · 🟡 = derived or inferred
> (reasoning given; arithmetic verified where possible) · 🔴 = not published anywhere.
> Any claim not marked 🟢 must reach the final page carrying its uncertainty.

> **Scope warning.** This is primarily a *bioinformatics* paper. The wet-lab chemistry is ~6
> paragraphs and delegates heavily to the Atrandi manuals. Much of what a base-by-base page needs
> is **not stated** — see §7.

> v2 change log: "Improved in-SPC lysis protocol and developed new **PTA-based scMetaG** method";
> compared scMetaG to bulk metaG; full software re-write. The PTA arm is new in v2.

---

## 1. Workflow

Six steps (Fig 2a), "adapted from Atrandi":

(a) encapsulate cells in SPCs at Poisson singlet dilution → (b) lysis → (c) WGA (**MDA or PTA**) →
(d) debranching → (e) combinatorial split-pool ligation barcoding (4 rounds) → (f) release,
fragmentation, index PCR.

## 2. Samples

- **Mock:** ATCC **MSA-2003**, ten strains, equal abundance.
- **Saliva:** healthy donor, ethical permit #08-047M. 5 mL into a 50 mL falcon, vortexed 3 min with
  glass beads (Fisherbrand 41401001), filtered 35 µm (Corning 352235), 5,000 g × 10 min, pellet
  washed 3× in 1 mL PBS, resuspended in 10% glycerol, −80 °C.

## 3. Encapsulation

- Dilute in PBS to **λ ≤ 0.1**, counted by hemocytometer after 10× SYTO9 (Invitrogen S34854).
- **ONYX platform (Atrandi, #MHN-ONYX1)** — *not* FLUX — with the **SPC Innovator Kit (#CKN-G-12)**.
  C2-chip.
- Saliva pretreated with **DNase I (0.1 U/µL, NEB, 37 °C, 30 min)** to destroy host/free DNA.
- **0.05% xanthan gum (Sigma G1253)** added to the core solution to stop bacteria aggregating or
  sticking to channel walls. ⚠ This introduces *Xanthomonas campestris* DNA — filtered out
  bioinformatically against GCF_001372255.1.
- Encapsulation run **3 h**; excess oil removed; shell crosslinked with **405 nm UV, 30–60 s**
  (Atrandi Light Exposure Device); emulsion broken; washed in Capsule Wash Buffer.
- Occupancy checked by fluorescence imaging with 10× SYTO9.

## 4. Lysis — custom protocol "R4"

Replaces the Atrandi lysis entirely; inspired by Li et al. (bioRxiv 2023.08.08.551713).
The Atrandi-style 0.8 M KOH / 20 mM EDTA / 200 mM DTT alkaline buffer is **dropped**.

1. Fix SPCs in **chilled 100% methanol**.
2. Enzyme cocktail in **TE (pH 8)**: 5 U lysostaphin (Sigma L9043), 50 U mutanolysin (Sigma SRE0007),
   50 U lysozyme (Biosearch E-0057-D2), 0.5 mg/mL achromopeptidase (Sigma A3547),
   2.5 mM EDTA, 10 mM NaCl. **37 °C, 800 rpm, overnight.**
3. Wash 3× WB; then proteinase K solution in TE: 4 U proteinase K (NEB P8107AA), 1% Triton X-100
   (Sigma X100), 100 mM NaCl. **55 °C, 1 h, 800 rpm.**
4. Pellet, resuspend in **1 mL 4 M guanidine thiocyanate, RT 30 min.**
5. Wash 3× WB, 2× nuclease-free water.

WB = 10 mM Tris-HCl pH 7.5 + 0.1% Triton X-100.

Dataset ↔ lysis protocol map (Supplemental Note Table 3):

| Protocol | Datasets |
|---|---|
| R1 | mock1_MDA |
| R2 | mock2_MDA |
| R3 | mock3_MDA, saliva1_MDA |
| **R4** | mock4_MDA, **mock5_PTA**, saliva2_MDA |

⚠ **The only PTA dataset is `mock5_PTA`.** The headline 15k-cell saliva result (Fig 4) is
`saliva2_MDA` — generated with **MDA**, even though the paper recommends against it.

## 5. Whole-genome amplification

### Option A — MDA (Atrandi kit)
Single-Microbe DNA Barcoding Kit (**#CKP-BARK1**). **45 °C 15 min → 65 °C 10 min.**
Wash 3× WB; verify by SYTO9 imaging. Then "DNA **debranching** and DNA end preparation as per the
instructions in the manual."

> 🟡 **INFERRED.** 45 °C is far above classic phi29 MDA (30 °C, hours), suggesting a thermostable
> phi29 variant — the preprint's own cost table lists "Equiphi29 250 U" under a DIY-MDA line item.
> Not stated as such anywhere.

### Option B — PTA ← **the key deviation**
**BioSkryb Genomics ResolveDNA Whole Genome Single-Cell Core Kit, #100954**
(cost table names it "ResolveDNA Whole Genome Amplification Kit v2.0"), per manufacturer protocol:

1. Lysis mix, **RT, 1400 rpm, 15 min.**
2. Add reaction mix, mix **1400 rpm, 1 min.**
3. Aliquot into PCR tubes: **30 °C, 2.5 h → 65 °C, 5 min → 4 °C hold.**
4. Then "SPCs were processed for **DNA end preparation** as per the instructions in the Atrandi manual."

> ⚠ **Debranching is absent from the PTA arm.** The MDA section names "debranching *and* end
> preparation"; the PTA section names end preparation only. Chemically coherent — PTA's terminators
> prevent the hyperbranching that debranching exists to resolve. Quotable deviation.
>
> 🟢 **Independently confirmed.** Negreira et al. (*bioRxiv* 2025.09.10.675331) substituted PTA into
> the same Atrandi SPC + split-pool workflow and state it outright: *"The samples submitted to PTA
> … **did not require debranching** and were directly submitted to barcoding instead."* They did
> still run end-prep (20 °C 30 min → 65 °C 30 min). See `06_pta.md` §6.
>
> Rationale given: "MDA library complexity saturates at 1x coverage… PTA achieved more even genomic
> coverage… We thus **strongly recommend against MDA for scMetaG**, and instead urge users to adopt
> PTA." Fig 2d: "PTA improves upon MDA by adding a blocking dNTP, preventing cyclic amplification,
> resulting in more amplicons of input template origin."
>
> Volumes and SPC input per PTA reaction: **NOT STATED.**

## 6. Barcoding, library prep, sequencing

### Barcoding — the entire Methods section, verbatim
> "Four-step combinatorial split-and-pool barcoding was performed to label each SPC with a unique
> variant. The barcodes were added by ligation following the Atrandi barcoding user guide, yielding
> up to 24⁴ barcode combinations."

4 rounds × 24 = **331,776** combinations. Ligase, temperatures, times, barcode and linker
lengths: **NOT STATED** — see `01_barcoding_kit.md`.

### Library prep
Release reagent (Atrandi) → **0.8× AMPure XP (Beckman A63881)** →
**NEBNext Ultra II FS DNA Library Prep Kit for Illumina (NEB #E7805S)** + **custom PCR indexing
primers** → double-sided AMPure → Qubit dsDNA HS (Thermo Q32854) → Agilent Bioanalyzer.

> 🟢 Confirms library prep is run **separately** with the NEBNext FS kit. But 🔴 **which adapter was
> ligated is NOT STATED**, nor the custom index primer sequences, PCR cycle count, fragmentation
> time, or target insert size. 🟡 **INFERRED** from the use of E7805S: the NEBNext hairpin adapter +
> USER workflow rather than the Atrandi single-tailed adapter — but this is only an inference from
> the kit choice, and the two imply *different* constructs. Worth confirming before drawing.

### Our index PCR primers 🟢 (as ordered — `ref/our_index_primers.tsv`)

Shared by both our protocols: this MDA/PTA scWGS and the work-in-progress
[florian-PTA-rnaseq design](https://henriksson-lab.github.io/chem_florian/01_primers.html).

```
i7 + P7:  CAAGCAGAAGACGGCATACGAGAT  [10-nt i7]  GTGACTGGAGTTCAGACGTGTGCTCT*T     61 nt
          +------- P7, 24 nt ------+             +-- Read-2 primer, 5' 27 nt --+
P5:       AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA*C                       45 nt
```

| Primer | i7 in oligo | ✅ i7 read / sample sheet |
|---|---|---|
| UDP0005_p7 | `TAATCTCGTC` | `GACGAGATTA` |
| UDP0006_p7 | `GCGCGATGTT` | `AACATCGCGC` |
| UDP0007_p7 | `AGAGCACTAG` | `CTAGTGCTCT` |
| UDP0008_p7 | `TGCCTTGATC` | `GATCAAGGCA` |
| generic P5 | — (no i5) | i5 read gives `ACCGAGATCT` (!) |

- ✅ **Atrandi's own primers** (`04_nebnext_illumina.md` §5) with one change: the 6-nt i7
  `AACCTG` becomes a **10-nt IDT UDP i7**. P5 is byte-identical to Atrandi's. `*` = 3'-PS.
- 🟡 The names match IDT for Illumina UD Indexes UDP0005–0008; the oligo carries the revcomp of
  the sample-sheet i7.
- ✅ **No i5 index.** A 10-cycle Index 2 read runs into P5 itself and reports
  `ACCGAGATCT` (sample-sheet orientation — the 10 nt of P5 before the shared `ACAC`), identical
  for every sample. Enter that as i5, or run i5 at 0 cycles. Consistent with §8's
  bcl2fastq `NNNNNN` + in-house demultiplexing.
- 🔴 PCR cycle count still not stated.

### Sequencing
🟢 **Illumina NovaSeq X**, 1.5B flow cell (one lane). 🔴 Read lengths for the single-cell libraries
are **NOT STATED**; bulk metagenomics was 2×150. 🟡 **INFERRED** 2×150 (R2 must exceed 45 nt).

## 7. Final construct (Fig 2c)

```
5'- P7 - i7 - [unlabelled block] - D - C - B - A - [genomic DNA] - Adapter - i5 - P5 -3'
```

🟢 Barcodes run **D, C, B, A** from the P7 side — reverse of ligation order. 🟡 The unlabelled block
between i7 and D is **inferred** to be the Read-2 primer site, and "Adapter" on the P5 side the
Read-1 site — the figure labels neither. 🟢 Cell barcode is read in **R2**; 🟢 **no UMI**
(`umi_from == umi_to == 0` in Bascet).

## 8. Read structure (authoritative — from Bascet source)

`crates/bascet-cli/src/barcode/atrandi_wgs_barcode_illumina.rs`:

```rust
pools[3].pos_anchor = (8 + 4) * 0;   // barcode D at R2 offset 0
pools[2].pos_anchor = (8 + 4) * 1;   // barcode C at R2 offset 12
pools[1].pos_anchor = (8 + 4) * 2;   // barcode B at R2 offset 24
pools[0].pos_anchor = (8 + 4) * 3;   // barcode A at R2 offset 36
trim_bcread_len = 8+4+8+4+8+4+8+1;   // = 45
//8 barcodes, 3 spacers, and 1 to account for ligation
```

| R2 offset (0-based) | Len | Content |
|---|---|---|
| 0–7 | 8 | Barcode **D** (round 4) |
| 8–11 | 4 | linker |
| 12–19 | 8 | Barcode **C** (round 3) |
| 20–23 | 4 | linker |
| 24–31 | 8 | Barcode **B** (round 2) |
| 32–35 | 4 | linker |
| 36–43 | 8 | Barcode **A** (round 1) |
| 44 | 1 | 🟡 ligation junction (inferred: the dA/dT T); trimmers add 1-2 nt of safety margin |
| 45– | — | genomic insert |

> The paper says R2 is trimmed by "the expected barcode size, **plus 2bp**"; the code trims **+1**.
> **Not a discrepancy — this is safety margin**, and the manual is not exact here. The construct
> itself has a single ligation-junction base. No action needed.

Matching: ≤1 mismatch per 8 nt barcode (`part_distance_cutoff = 1`), ≤4 total
(`total_distance_cutoff = 4`); each barcode searched at anchor and anchor+1
(`pos_rel_anchor = [0,1]`). Pairs discarded if post-trim length ≤ 60.
Barcode A is looked up first — "the most likely one to fail in case of over-fragmentation".

Barcode positions are not hard-coded blindly: they are **detected from the first 10k reads** by
scanning for whitelist entries, with ±1 bp slack per barcode and 4 bp total.

### Trimming
R1 is the genomic read and must be trimmed when the insert is short enough to read into the
barcode. "This analysis is complicated by the fact that **there is no fixed adapter sequence**."
fastp (v0.23.4) was insufficient; they wrote their own trimmer that reverse-complements a **16 bp
window** and scans for it in the opposite read, deliberately scanning into the barcode for
confidence on short reads.

### Sample index
bcl2fastq run with **`NNNNNN` as the library index** plus `--create-fastq-for-index-reads`, so
R1/R2/I1 come out and the sample index is demultiplexed in-house.

## 9. Barcode whitelist

96 sequences, 24 per round, 8 nt, in
<https://github.com/henriksson-lab/bascet/blob/main/crates/bascet-cli/src/barcode/atrandi_barcodes.tsv>
(archived as `ref/bascet_atrandi_barcodes.tsv`). Plate: round 1 = wells A1–H3, round 2 = A4–H6,
round 3 = A7–H9, round 4 = A10–H12.

⚠ **This list disagrees with the Atrandi manual for round 1** (12/24 sequences differ) and uses a
transposed well mapping throughout. See `02_library_prep.md` §10 for the full comparison.

These are the sequences **as read in R2**, not the ordered oligos — no 5'-phosphate, no linker arms,
no partner strands.

## 10. Software, data, references

| | |
|---|---|
| Bascet (Rust, MIT) | <https://github.com/henriksson-lab/bascet> — version #c378038, v0.0.2 |
| Zorn (R, MIT) | <https://github.com/henriksson-lab/zorn> · docs <http://zorn.henlab.org/> |
| Study code | <https://github.com/henriksson-lab/scwgs2026> |
| Saliva scMetaG + bulk | BioStudies **S-BSST3416** |
| Mock scMetaG | BioStudies **S-BSST3417** |
| Simulated mock | Zenodo [10.5281/zenodo.15813069](https://doi.org/10.5281/zenodo.15813069) |

Chemistry selectors in the CLI: `AtrandiWGS` (paired-end Illumina) and `AtrandiWGSLR`
(long-read PacBio/Nanopore, single-read — **present in code, no long-read experiment in the paper**).

Contaminant filtering: human T2T-CHM13v2.0; *X. campestris* B-1459 GCF_001372255.1.

## 11. Deviations from the stock Atrandi protocol

1. **WGA replaced** — Atrandi MDA → BioSkryb ResolveDNA **PTA** (#100954), 30 °C/2.5 h.
2. **Debranching omitted** in the PTA arm.
3. **Custom lysis R4** — methanol fixation + 4-enzyme cocktail + proteinase K + 4 M GTC;
   Atrandi alkaline lysis dropped.
4. **Library prep run separately** with NEBNext Ultra II FS (#E7805S) + custom index primers.
5. **0.05% xanthan gum** in the core solution (Pranauskaite et al. 2024).
6. **DNase I pre-treatment** of saliva before encapsulation.
7. **Custom debarcoder/trimmer (Bascet)**; fastp explicitly rejected.
8. **Sample index demultiplexed in-house.**

Explicitly *unchanged*: SPC generation (ONYX + CKN-G-12 per manual), the 4-round split-pool
ligation barcoding, DNA end-prep, and the Release reagent step.

## 12. Still missing for a base-by-base page

**Resolved by decision — not blockers:**

- **Barcode sequences.** We do *not* need the real A/B/C/D 8-mers. The page uses placeholders,
  exactly as scg_lib_structs does (`<cbc>[8-bp Round1 barcode]</cbc>`). Draw them as
  `AAAAAAAA` / `BBBBBBBB` / `CCCCCCCC` / `DDDDDDDD` so the four rounds stay distinguishable, with
  bracketed prose in the oligo list. The whitelist conflict in `02_library_prep.md` §10 is therefore
  irrelevant to the page — it matters only for demultiplexing.
- **The 4 nt linkers** are likewise drawn as placeholders unless real sequences turn up.

**Genuinely still open — these change what the diagram shows:**

1. **The barcode cassette architecture.** Not the barcode bases, but the *oligo design*: how many
   strands, where the overhangs and 5'-phosphates sit, and — critically — **the Read-2 / P7 arm that
   the cassette installs**, since `04_nebnext_illumina.md` §5 shows the ligation adapter carries no
   P7 arm at all. Without this the left half of the construct cannot be drawn.
2. ~~Our custom i5/i7 index primer sequences~~ — resolved, §6. PCR cycle count still open.
3. NEBNext FS fragmentation time / target insert size.
4. Single-cell run read lengths (R1/R2/i5/i7).

Items 2 (cycle count) – 4 are ours to state. Item 1 is the one real external unknown.

**Best candidate sources for item 1:** Baronas et al., *Science* (2026),
doi:[10.1126/science.ady7227](https://doi.org/10.1126/science.ady7227) — the Atrandi SPC platform
paper; Mazelis et al., *Science* (2026),
doi:[10.1126/science.ady7209](https://doi.org/10.1126/science.ady7209); Ling et al.,
*Front. Microbiol.* **15**, 1516656 (2025); or Atrandi directly.

## 13. qPCR primers (Supplemental Note Table 2)

Lysis-verification primers only — not part of the library chemistry. Master mix:
PCRBIO HS Taq Mix Red (PB10.23-02).

| Species | Fwd (5'→3') | Rev (5'→3') | Target | bp |
|---|---|---|---|---|
| *C. beijerinckii* | TGACACGATTTTTCATTCTCCA | TCCATTGCCTTAATGACAGGT | nifH | 448 |
| *E. coli* | CGTGGTGATTGATGAAACTG | TGATACATATCCAGCCATGC | uidA | 564 |
| *L. gasseri* | ATCACATTCAACTCTCGCTG | TCATTCATCTTCATCGTCCT | LGAS_0517 | 400 |
| *S. epidermidis* | GATATTCGCGATGAACTTGC | ATCAGGTGTTGCAAATAGGG | fibrinogen-binding | 414 |
| *S. mutans* | TCGCGAAAAAGATAAACAAACA | GCCCCTTCACAGTTGGTTAG | species-specific | 479 |
| *Rhodobacter*/*Ceribacter* sp. | ATGTTGGACCTCCGCAAAG | TGCCATCTGATGCCGTATTG | Rsph17029_0004 | 341 |
| *E. faecalis* | AGTACCATTCGTGCCAGTTT | GCGTATTCTTGCGCTTGATG | ddl | 399 |
| *D. radiodurans* | GAAGTCGAGGTGGCGTTTAT | TCTTGCCGATCTTGGGATTT | gyrB | 361 |
| *B. adolescentis* | GGTGATTACGCAGCATCCTT | CTTCCTCACAAACGTCAGCA | — | — |
| *B. pacificus* (*cereus*) | GTGGTTCTGCTGTATCTA | CAGCACCAGTAACGTTTA | entFM | 183 |
