# FIPRESCI — Tn5 pre-indexing of RNA/cDNA hybrids, then overloaded 10x 5' droplets

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`). Claims that rest only
> on the upstream scg_lib_structs page are 🟡 (secondary source).

**FIPRESCI** — Li Y, Huang Z, Zhang Z, Wang Q, Li F, Wang S, Ji X, Shu S, Fang X, Jiang L.
"FIPRESCI: droplet microfluidics based combinatorial indexing for massive-scale 5′-end
single-cell RNA sequencing." *Genome Biology* 24:70 (2023).
doi:[10.1186/s13059-023-02893-1](https://doi.org/10.1186/s13059-023-02893-1) ·
PMID 37024957 · PMC10078054 (CC-BY 4.0).

Other papers on the method:

- Publisher Correction, *Genome Biol.* 24:88 (2023),
  doi:[10.1186/s13059-023-02944-7](https://doi.org/10.1186/s13059-023-02944-7) (PMC10122391).
  It fixes a typesetting error in Fig. 1 (figure layers lost in production; the legend is
  reprinted unchanged) and adds a Human Cell Atlas acknowledgement — nothing about the
  chemistry. 🟢
- Preprint: "Droplet microfluidics based combinatorial indexing for massive-scale 5′-end
  single-cell RNA sequencing", Research Square,
  doi:[10.21203/rs.3.rs-1510072/v1](https://doi.org/10.21203/rs.3.rs-1510072/v1) — not fetched;
  the DOI is not in any source read here and was not checked. 🔴
- The same group later generalised the idea to several modalities as UDA-seq
  (*Nat Methods* 2024, doi:10.1038/s41592-024-02586-y) — not read; mentioned only as
  the follow-up; citation not checked against any fetched source. 🔴

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/fipresci__10.1186+s13059-023-02893-1/`, never committed):

| File | What | Used for |
|---|---|---|
| `s13059-023-02893-1_PMC10078054.xml` → `.xml.txt` | full text, JATS XML from Europe PMC (text twin made by hand, see §9) | **methods**, results, conditions |
| `…_supp_13059_2023_2893_MOESM2_ESM.xlsx.txt` | Supplementary Table 1, sheets "Fipresci-Seq oligos" and "96 i7-Tn5" | **every oligo** |
| `…_MOESM1_ESM.pdf.txt` | Additional file 1 (Supplementary Figs.), esp. **Fig. S2** "Detailed assay design" | reaction-by-reaction sequence diagram, bead TSO |
| `…_MOESM3_ESM.xlsx.txt` | Supplementary Table 2: round-1 barcode → sample, per experiment | which wells were what |
| `…_MOESM4_ESM.xlsx.txt` | Supplementary Table 3: 16 tagmentation buffers; per-experiment conditions | buffer, RT primer and loading per experiment |
| `…_MOESM5_ESM.xlsx.txt` | Supplementary Table 4: donor clinical data | not chemistry |
| `…_MOESM6_ESM.docx.txt` | legends of processed data files on figshare | not chemistry |
| `…_MOESM7_ESM.pdf.txt` | Review history (reviewer comments, authors' replies) | TCR barcode loss, which strand Tn5 tags |
| `upstream_FIPRESCI.html.txt` | scg_lib_structs page | second account, checked in §7 |
| `s13059-023-02944-7_PMC10122391_correction.xml(.txt)` | the Publisher Correction, fetched by hand from Europe PMC | what it corrects (Fig. 1 artwork, acknowledgement) |

Could not be read:

| What | URL | Why |
|---|---|---|
| PMC article HTML | https://pmc.ncbi.nlm.nih.gov/articles/PMC10078054/ | saved file is a reCAPTCHA page; the Europe PMC XML replaces it |
| Figure images (`Fig1–5_HTML.jpg`) | in the supplementary zip | images; not read |

---

## 1. What it is

**Two rounds of combinatorial indexing, the second done by an unmodified 10x Genomics 5'
gel bead** 🟢 (main text, "FIPRESCI overview"; Fig. S2):

1. Fixed/permeabilised cells or nuclei are reverse transcribed **in situ**, in bulk, with
   oligo-dT and/or random hexamers. The RNA stays hybridised to its first-strand cDNA.
2. They are split over a 96-well plate. Each well holds a Tn5 loaded with **one** adapter
   type (an "i7-only" homodimer) carrying a **6-nt well barcode (round 1)**. Tn5 cuts the
   **RNA/cDNA hybrid** inside the cell.
3. All wells are pooled and **overloaded** into one 10x Chromium 5' channel (15,300 up
   to ~200,000 cells/nuclei). In the droplet the cDNA's 3' C-tail template-switches onto
   the gel-bead TSO, adding the **16-nt droplet barcode (round 2) and the 10-nt UMI**.
4. Only the **5'-most fragment** of each transcript carries both handles (round 1 at the
   cDNA 5' end, round 2 at its 3' end) and is amplified. 🟡

Cell = (round-1 barcode, round-2 barcode). Droplets with several cells are resolved, not
discarded; species mixing at 15,300 loaded gave a 0.2 % collision rate. 🟢

| | New thing here | Builds on / where else |
|---|---|---|
| 1 | **Tn5 tagmentation of RNA/cDNA hybrids inside intact cells** as the first-round barcode, and the observation that Tn5 stays bound (fragments only after SDS), so the tagged molecule stays in the cell 🟢 | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md); Tn5 on RNA/DNA hybrids from refs 19–20 of the paper |
| 2 | **Round 2 from the droplet TSO**, so the 5' end of the transcript is what is read, and TSS / eRNA positions come out of the data | [template switching](../ref/concepts/template-switching.md); the 10x 5' GE kit (`../10x-chromium-single-cell-5-ge/01_10x-5p-ge.md`) |
| 3 | Overloading a droplet system with pre-indexed cells | scifi-RNA-seq (`../scifi-rna-seq__10.1101+2019.12.17.879304/01_scifi-rna-seq.md`), which pre-indexes at the RT primer and is therefore 3'-end |
| 4 | RT **before** Tn5, and the RT is done once in bulk, not in the droplet | [reverse transcription](../ref/concepts/reverse-transcription.md) |

## 2. Oligos

🟢 Verbatim from Supplementary Table 1 (`MOESM2_ESM.xlsx`, sheet "Fipresci-Seq oligos"),
modifications as written there: `[5Phos]` = 5' phosphate, `/3ddC/` = 3' dideoxy-C
("to block 3' extension"), `[5Bio]` = 5' biotin, `-s-` = phosphorothioate bond. All HPLC
purified.

```
TN5_A_ME              [5Phos]CTGTCTCTTATACACATCT/3ddC/
TN5_R2_index          GACGTGTGCTCTTCCGATCT[NNNNNN]AGATGTGTATAAGAGACAG
RT PolyT Primer       d(T)23VN                      (NEB #S1327S)
RT random Primer      d(NNNNNN)                     (Thermo Scientific # SO142)
S5R-P5                AATGATACGGCGACCACCGA-GATCTACACTCTTTCCCTACACGACGCTC
S5R-P5-bio            [5Bio]AATGATACGGCGACCACCGA-GATCTACACTCTTTCCCTACACGACGCTC
S-P7                  CAAGCAGAAGACGGCATACGAGAT[NNNNNN]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
Human TCR outer-1     TGAAGGCGTTTGCACATGCA
Human TCR outer-2     TCAGGCAGTATCTGGAGTCATTGAG
Humane TCR inner-1    CTGGTTGCTCCAGGCAATGG
Human TCR inner-2     TGTAGGCCTGAGGGTCCGT
```

(The table's hyphen inside S5R-P5 is a typographic break, not a modification; the
methods write the same primer without it: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`.)

**Twelve sample-index primers** 🟢 (same sheet), 64 nt each:

```
S-P7-index1   CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index2   CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index3   CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index4   CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index5   CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index6   CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index7   CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index8   CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index9   CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index10  CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index11  CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
S-P7-index12  CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T
```

The in-vitro hybrid test in the methods uses an **un-indexed** S-P7 🟢:
`CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC-s-T`.

**96 round-1 Tn5 oligos** 🟢 (sheet "96 i7-Tn5", named `n-Tn5-bo_i7-index1..96`), each
`GACGTGTGCTCTTCCGATCT` + barcode + `AGATGTGTATAAGAGACAG`, e.g.
index1 = `GACGTGTGCTCTTCCGATCTAAAGAAAGATGTGTATAAGAGACAG`. The barcodes, in order:

```
 1-12  AAAGAA AACAGC AACGTG AAGCCA AAGTAT AATTGG ACAAGG ACCCAA ACCTTC ACGGAC ACTGCA AGACCC
13-24  AGATGT AGCACG AGGTTA AGTAAA AGTCTG ATACTT ATAGCG ATATAC ATCCGG ATGAAG ATTAGT CAACCG
25-36  CAAGTC CACCAC CACTGT CAGACT CAGGAG CATAGA CCACGC CCGATG CCGTAA CCTCTA CGAAAG CGAGCA
37-48  CGCATA CGGCGT CGGTCC CGTTAT CTAGGT CTATTA CTCAAT CTGTGG CTTACG CTTGAA GAAATA GAAGGG
49-60  GACTCG GAGCTT GAGGCC GAGTGA GATCAA GCCAGA GCCGTT GCGAAT GCGCGG GCTCCC GCTGAG GCTTGT
61-72  GGACGA GGATTG GGCCAT GGGATC GGTAGG GGTGCT GTACAG GTCCTA GTCGGC GTGGTG GTTAAC GTTTCA
73-84  TAAGCT TAATAG TACCGA TAGAGG TATTTC TCAGTG TCATCA TCCAAG TCGCCT TCGGGA TCTAGC TGAATT
85-96  TGAGAC TGCGGT TGCTAA TGGCAG TGTGTA TGTTCG TTAAGA TTCGCA TTCTTG TTGCTC TTGGAT TTTGGG
```

**10x 5' gel-bead TSO** (not a FIPRESCI oligo; drawn in Fig. S2 step 4) 🟢:
`5’-CTACACGACGCTCTTCCGATCT-(10X-round2)-NNNNNNNNNN-TTTCTTATAT-rGrGrG` — round 2 is the
16-nt 10x barcode (length from the upstream page and the 10x guide, 🟡).

### How they interlock — 🟡 (computed)

- **TN5_A_ME is exactly `nextera.ME_RC`** (19 nt) and pairs with the 3'-terminal ME of
  every TN5_R2_index oligo: the usual Tn5 adapter duplex, with the 5'-phosphorylated,
  3'-blocked non-transferred strand.
- **TN5_R2_index = `illumina.TRUSEQ_READ2[14:]` (its last 20 nt) + 6-nt barcode +
  `nextera.ME`** — 45 nt, for all 96. The adapter is **TruSeq Read 2**, not Nextera s7:
  the Tn5 is loaded with a single adapter type, so every cut end gets Read 2. The paper
  calls it "i7-only" and Fig. S2 "Tn5 s7 homodimer"; "s7" there means "the i7 side".
- The 96 barcodes are all distinct, minimum pairwise Hamming distance **3**.
- **S-P7-index = `illumina.P7` + 6-nt index + the full 34-nt `illumina.TRUSEQ_READ2`**
  (last T phosphorothioate-protected). Its 3' 20 nt are the TN5_R2_index head, so it
  primes directly on the Tn5 end of the fragment.
- The 12 indices as written in the oligos (`CGTGAT`, `ACATCG`, …) are the **reverse
  complements** of what the i7 read reports (`ATCACG`, `CGATGT`, `TTAGGC`, `TGACCA`,
  `ACAGTG`, `GCCAAT`, `CAGATC`, `ACTTGA`, `GATCAG`, `TAGCTT`, `GGCTAC`, `CTTGTA`) — these
  read-orientation sequences look like Illumina TruSeq LT indices 1–12 (recognised by
  eye; no TruSeq index list in `lib/` to check against).
- **S5R-P5 = `illumina.P5` + `TCTTTCCCTACACGACGCTC`** — 49 nt; equals
  `illumina.NEBNEXT_UNIVERSAL_PRIMER` minus its last 9 nt (`TTCCGATCT`). (P5's last four
  bases `ACAC` are TruSeq Read 1's first four; the 20 nt after them are
  `TRUSEQ_READ1[4:24]`.)
- The bead TSO's 22-nt head `CTACACGACGCTCTTCCGATCT` is the 3' end of
  `illumina.TRUSEQ_READ1`. **S5R-P5's 3' end overlaps it by only 13 nt**
  (`CTACACGACGCTC`): the biotinylated primer anneals by 13 nt and adds the rest of P5 /
  Read 1 as a 5' tail.
- The two "Human TCR outer" primers are identical to 10x **Human T Cell Mix 1 v2**
  reverse outer primers; the two "Human/Humane TCR inner" primers are identical to the
  two primers listed under 10x **Mouse T Cell Mix 1 v2** (reverse outer) in the 10x 5' v2
  user guide CG000331 Rev E, not to the human inner (Mix 2) primers — string match against
  `_data/sources/10x-chromium-single-cell-5-vdj/CG000331_5v2_UserGuide_RevE.pdf.txt` and
  `../10x-chromium-single-cell-5-vdj/01_chromium-5-vdj.md`. The methods say the VDJ
  enrichment used the 10x kit; whether these table entries are what was used, or a
  labelling error in the table, is not stated. 🔴

## 3. Step by step

Volumes and conditions 🟢 from the main-text Methods unless noted.

1. **Transposome assembly.** TN5_A_ME + one TN5_R2_index, 10 µM each 1:1 in Vazyme
   Annealing Buffer: 75 °C 15 min, 60 °C 10 min, 50 °C 10 min, 40 °C 10 min, 25 °C 30
   min. Then 7 µl cassette + 4 µl TruePrep Tagment Enzyme (Vazyme S601-01) + 39 µl
   coupling buffer, 30 °C 1 h. One transposome per barcode (96).
2. **Cells or nuclei.** Either methanol-fixed permeabilised cells (90 % ice-cold methanol,
   −20 °C 10 min) or nuclei (NP-40/Tween-20/digitonin lysis; for the E10.5 embryo the
   nuclei were additionally fixed in 2 % formaldehyde 10 min on ice).
3. **In-situ RT.** 100,000 cells/nuclei in 7 µl + 3 µl 25 µM RT primer (oligo-dT, random
   hexamer, or both); 55 °C 5 min, ice. Add 40 µl: 10 µl 5× RT buffer, 2.5 µl 100 mM DTT,
   2.5 µl 10 mM dNTPs, 2.5 µl RNaseOUT, 3.5 µl **Maxima H Minus** RT, 21.5 µl water.
   50 °C 10 min; 3 × [8 °C 12 s, 15 °C 45 s, 20 °C 45 s, 30 °C 30 s, 42 °C 2 min,
   50 °C 3 min]; 50 °C 5 min. (The low-temperature ramp presumably lets random hexamers
   anneal and extend 🟡.)
   Product: RNA/cDNA hybrid, and the M-MLV-type RT leaves 2–5 untemplated C on the cDNA
   3' end 🟢 (overview) — the C's that later pair with the bead TSO.
   No TSO in this step: template switching happens only in the droplet. 🟡
4. **Tagmentation in the plate.** 2,000–4,000 cells/nuclei (1–2 µl) into each well of
   18 µl: 1× tagmentation buffer + **1 µM** indexed transposome; 37 °C 30 min, 1000 rpm.
   Stop with 5 µl 0.5 M EDTA, 10 min on ice. Buffer: early experiments used Vazyme
   "Tagment buffer L"; after the 16-buffer screen (Supplementary Table 3) the default was
   **Tris-DMF**, 10 mM Tris-HCl pH 7.5, 5 mM MgCl₂, 10 % DMF. Crowding agents (PEG 8000 /
   PEG 200) did *worse* in cells, unlike on free hybrids. 🟢
   On the 5'-terminal fragment the transferred strand (TN5_R2_index) is joined to the
   **5' end of the cDNA** fragment; TN5_A_ME stays annealed but not ligated (9-nt gap).
   🟡 (standard [Tn5](../ref/concepts/tn5-tagmentation.md) strand transfer; the authors'
   reply to Reviewer 3, `MOESM7`, argues the same)
5. **Pool**, add 240 µl 10 % BSA (to 1 %), spin 600 g 5 min, wash twice with PBS-BSA-RNase
   inhibitor, resuspend in 20 µl.
6. **GEMs** on one Chromium lane with the 10x 5' kit (V(D)J v1.1, CG000207; or 5' v2,
   CG000331). Loaded: 15,300 (species mix), 100,000 (cell lines), 80,000 (T cells),
   ~200,000 nuclei (embryo). Custom incubation: **25 °C 30 min, 42 °C 90 min, 53 °C 10
   min**. In the droplet the C-tailed cDNA anneals to the bead TSO's `rGrGrG` and is
   extended across the TSO, so the cDNA 3' end gains barcode + UMI + Read 1 head 🟢
   (Fig. S2 step 4–5). Fig. S2 also draws the TSO strand extended over the whole insert
   and Tn5 end. 🟢 (drawn) / mechanism not discussed 🔴
7. **GEM cleanup**: Recovery Agent, Dynabeads MyOne Silane, elute 35.5 µl (kit procedure).
8. **Linear enrichment with S5R-P5-bio** (single primer): 35 µl + 50 µl NEBNext HiFi 2× +
   0.5 µl 100 µM S5R-P5-Bio + 14.5 µl water; 72 °C 3 min, 98 °C 45 s, 13–16 × [98 °C
   20 s, 67 °C 30 s, 72 °C 1 min], 72 °C 1 min. 0.8× AMPure. Only molecules with the TSO
   handle are copied, and every copy is 5'-biotinylated. 🟡
9. **Streptavidin capture**: 10 µl MyOne Streptavidin C1 in B&W, 60 min RT rotating,
   wash, resuspend in 20 µl water. The beads go straight into the next PCR.
10. **Library PCR on beads**: 20 µl beads + 50 µl KAPA HiFi HotStart ReadyMix + 0.5 µl
    100 µM S5R-P5 + 5 µl 10 µM S-P7-index + 24.5 µl water; 98 °C 45 s, 16 × [98 °C 20 s,
    54 °C 30 s, 72 °C 20 s], 72 °C 1 min. 0.75× AMPure.
11. **Optional TCR**: 10x nested V(D)J PCR (two rounds) on the FIPRESCI cDNA, then the 10x
    V(D)J library prep; sequenced at ~5,000 reads/cell. 🟢

Which fragments survive 🟡 (upstream page and Fig. S2 agree): internal fragments of the
hybrid have Tn5 adapters but no C-tail and are not template-switched; the 3'-terminal
fragment ends in oligo-dT, not C's. Only the **5'-terminal fragment** (cDNA 3' end with
the untemplated C's) gets both the Tn5 Read-2 handle and the TSO Read-1 handle.

## 4. Final library — 🟡 (assembled from the oligos above; matches Fig. S2 "Final Library" and upstream)

```
5'- P5 (illumina.P5, 29)
  · TCTTTCCCTACACGACGCTC (rest of S5R-P5, 20)
  · TTCCGATCT (from the bead TSO; P5's last 4 nt `ACAC` + these 20 + 9 nt = TruSeq Read 1, 33)
  · <round-2 droplet barcode, 16>
  · <UMI, 10>
  · TTTCTTATATGGG (TSO linker + GGG, 13)
  · <5' end of transcript>
  · CTGTCTCTTATACACATCT (ME_RC, 19)
  · <revcomp of round-1 barcode, 6>
  · AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC (revcomp TruSeq Read 2, 34)
  · <i7 sample index, 6>
  · ATCTCGTATGCCGTCTTCTGCTTG (illumina.P7_RC, 24) -3'
```

The cDNA insert is in **sense** orientation on this strand (Read 1 starts at the TSS side).
The round-1 barcode on this strand is the reverse complement of the barcode written in the
oligo table; Read 2 reads it as written. 🟡 (computed)

TCR library: the nested PCR primes in the TCR constant region, downstream of the 5'
variable region, so the **Tn5 end — and with it the round-1 barcode — is cut off** 🟢
(main text and authors' reply). Only droplets that the gene-expression data show to hold a
single cell are kept; for those the 10x barcode alone identifies the cell. 🟢

## 5. Sequencing

🟢 NovaSeq 6000 or MGISeq-T7, **2 × 150**, 10,000–50,000 reads per cell. Read lengths
for the index read are not given. 🔴

Read layout 🟡 (computed from §4; agrees with upstream):

| Read | Primer | Reads |
|---|---|---|
| Read 1 (150) | `illumina.TRUSEQ_READ1` | 16 barcode · 10 UMI · `TTTCTTATATGGG` · transcript from its 5' end |
| i7 (≥6) | `illumina.INDEX1_PRIMER` | sample index, as the reverse complement of the oligo's 6 nt |
| Read 2 (150) | `illumina.TRUSEQ_READ2` | **round-1 barcode (6, as written in the table)** · ME `AGATGTGTATAAGAGACAG` (19) · cDNA (antisense) |

Processing 🟢: reads are split by round-1 barcode (no mismatch allowed), each subset run
through Cell Ranger 6.0.2 `--chemistry=fiveprime --include-introns` as if it were an
ordinary 10x 5' run, and the matrices merged in Seurat with the round-1 barcode added to
the 10x barcode.

## 6. Experiments, as conditions

From Supplementary Table 3 "Reaction conditions" 🟢:

| Experiment | Sample | RT primer | Tagmentation buffer | Loaded |
|---|---|---|---|---|
| Species mixing | Jurkat + NIH-3T3, permeabilised cells | oligo-dT | Tagment buffer L | 15,300 |
| 3 cell lines, nuclei / cells | HEK293T, Jurkat, K562 | oligo-dT | Tagment buffer L | 100,000 each |
| Buffer tests | HEK293T / HeLa cells | oligo-dT | 16 buffers | 15,300 |
| RT primer test | HeLa cells and nuclei | dT / random / mix | Tris-DMF | 20,000 |
| E10.5 embryo | fixed nuclei | oligo-dT | Tris-DMF | 200,000 |
| PBMC T cells | permeabilised cells | oligo-dT | Tris-DMF | 80,000 |

Results quoted for orientation 🟢: 8,049 cells from 15,300 loaded (0.2 % collisions);
58,771 cells and 59,622 nuclei from 100,000 loaded each; 117,804 nuclei from one channel
for the embryo. Note: the main text says the three cell lines were HEK293T, Jurkat and
K562 (as does Supplementary Table 3), but the Fig. 1 legend, the correlation values and
the data-file legends (`MOESM6`) say HEK293(T), HeLa and K562 — an internal
inconsistency of the paper. 🟡

## 7. Upstream (scg_lib_structs) checked against the paper

| Point | Upstream | Paper / supplement | Verdict |
|---|---|---|---|
| Oligo source | "Supplementary Table 1" | MOESM2 is Supplementary Table 1 | agree |
| TN5_A_ME | `/Phos/ CTGTCTCTTATACACATCT /3ddC/` | `[5Phos]…/3ddC/` | agree (notation only) |
| TN5_R2_index | Read-2 head + 6-nt barcode + ME | same; the 96 barcodes are listed only in the paper | agree; upstream omits the list |
| RT primers | `TTTTTTTTTTTTTTTTTTTTTTTVN` (T23VN), `NNNNNN` | `d(T)23VN`, `d(NNNNNN)` | agree |
| Bead TSO | `CTACACGACGCTCTTCCGATCT[16][10]TTTCTTATATrGrGrG`, "PN-220112" | Fig. S2: same sequence; no part number given | agree; part number is upstream's |
| S5R-P5 / -bio | as paper | as written | agree |
| S-P7-index | ends `…CTCTTCCGATCT`, no modification | ends `…CTCTTCCGATC-s-T` (phosphorothioate) | **upstream drops the phosphorothioate** |
| Sample indices | "[6-bp i7]" | 12 explicit S-P7-index sequences | upstream omits the list |
| TCR primers | not mentioned | 4 primers in Table 1 | upstream omits |
| Streptavidin step | step (6) "cDNA purification", no beads mentioned | MyOne C1 capture, library PCR on beads | upstream omits the bead step |
| GEM incubation, PCR cycles, buffers | not given | given (§3) | paper only |
| Read layout | R1 150: barcode, UMI, spacer, cDNA; i7 6; R2 150: Tn5 barcode + cDNA | PE150 only | consistent; index length is upstream's |
| Fragment fates | three products, only the 5' one captured | Fig. S2 shows the 5' product only | agree |

## 8. Open questions

- 🔴 The "Human TCR inner" primers in Table 1 are 10x's mouse TCR outer primers (§2); which
  primers were actually used for the human TCR enrichment.
- 🔴 How reads were handled before Cell Ranger: Read 2 starts with 6 nt barcode + 19 nt ME
  before the cDNA; whether these 25 nt were trimmed is not said.
- 🔴 How MGISeq-T7 sequencing was done with Illumina P5/P7 libraries (no conversion step is
  described).
- 🔴 Purpose of the 72 °C 3 min at the start of the S5R-P5-bio linear PCR (no gap needs
  filling on the strand that is copied).
- 🔴 Tn5 also cuts the RNA strand of the hybrid; the paper does not quantify the 5'
  fragments lost because the RNA, not the cDNA, carries the transferred adapter at the 5'
  end (Reviewer 3 asked; the authors argue it cannot happen with this adapter design).
- 🟡 Index identity with Illumina TruSeq LT 1–12 — recognised by eye, not checked in `lib/`.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` got the Europe PMC XML, the supplementary zip (7 files + figure
images) and the upstream page; the PMC HTML it saved is a reCAPTCHA page. `doctext.py`
made no text twin for the JATS XML (`--stdout` passes the XML through raw), so the
`.xml.txt` was made with a small ElementTree script. The Publisher Correction was found
by a Europe PMC search and fetched by hand. `scrape_primers.py --max-hits 80` put
Supplementary Table 1 first and recognised `illumina.P5`, `P7`, `TRUSEQ_READ2` and
`nextera.ME` in the oligos, but it did not parse the `-s-T` phosphorothioate (the S-P7
sequences come out one base short, ending `…CCGATC`) and did not report the `[5Phos]` of
TN5_A_ME; the 96 Tn5 oligos used most of the hit budget. `--find` likewise reports 0 locations for
the full S-P7 / S-P7-index sequences ending in `…CCGATCT`, because the source writes the
last base as `-s-T`; the sequences were confirmed with `--find` on `…CCGATC` plus a plain
grep for `CCGATC-s-T` (13 occurrences: S-P7 and the 12 indexed primers). The methods were read, not
scraped.
