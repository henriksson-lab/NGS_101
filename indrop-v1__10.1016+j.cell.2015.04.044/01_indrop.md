# inDrop: droplet barcoding with photo-released hydrogel-bead primers and CEL-Seq-style IVT

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or not in anything we could fetch. Relationships marked 🟡 *(computed)* were
> worked out with `lib/` while writing this note. They are not yet asserted in a
> self-test, because this protocol has no `tools/` module yet (status `notes`).
> Claims taken only from the upstream scg_lib_structs page are 🟡 *(upstream)*: a
> secondary source, not checked against the authors' own files.

**inDrop V1.** Klein AM, Mazutis L, Akartuna I, Tallapragada N, Veres A, Li V, Peshkin L,
Weitz DA, Kirschner MW. "Droplet barcoding for single-cell transcriptomics applied to
embryonic stem cells." *Cell* 161, 1187–1201 (2015).
doi:[10.1016/j.cell.2015.04.044](https://doi.org/10.1016/j.cell.2015.04.044) ·
PMID 26000487 · PMC4441768 (author manuscript). Data: GEO GSE65525 🟢 (authors and
accession from the fetched PMC page).

**inDrop V2 (the "inDrops" protocol paper).** Zilionis R, Nainys J, Veres A, Savova V,
Zemmour D, Klein AM, Mazutis L. "Single-cell barcoding and sequencing using droplet
microfluidics." *Nature Protocols* 12, 44–73 (2017).
doi:[10.1038/nprot.2016.154](https://doi.org/10.1038/nprot.2016.154). The title and
journal come from `catalogue/scg_lib_structs.tsv`. The author list is from memory 🟡,
because the main text was not fetched.

Later papers in the catalogue that build on it, not covered here: **inDrops-2**
(doi:10.1101/2023.09.26.559493) and **spinDrop** (doi:10.1101/2023.01.12.523500).

## Sources read

Fetched by `tools/get_sources.py "indrop__10.1016+j.cell.2015.04.044"` into
`_data/sources/indrop__10.1016+j.cell.2015.04.044/` (gitignored):

| File | What | Used for |
|---|---|---|
| `upstream_inDrop.html.txt` | scg_lib_structs page "inDrop V1 / inDrop V2&3" | the only source for the **V1 library oligos** (RLO, 2nd RT primer, V1 PCR and read primers), and the step drawings |
| `j.cell.2015.04.044_PMC4441768.html.txt` | Cell 2015 main text (PMC author manuscript) | bead design (Fig. 2 legend), 147,456 barcodes, encapsulation, UV release, RT incubation, demulsification, "processed as per CEL-Seq" |
| `nprot…MOESM455_ESM.xlsx.txt` | Nat Protoc Supplementary Table: sheet `oligos_except_barcodes` | **every V2 oligo** except the barcode plates |
| `nprot…MOESM456_ESM.xlsx.txt` | sheet `P1=W1rc-bc1.n-PE1rc`, 384 oligos | barcode-1 plate |
| `nprot…MOESM457_ESM.xlsx.txt` | sheet `P2=BA19-N6-bc2.n-W1rc`, 384 oligos | barcode-2 plate |
| `nprot…MOESM454_ESM.pdf.txt` | Supplementary Methods 1–2, Table 1, Fig. 1 | bead QC gel (expected band sizes), cost |
| `nprot…MOESM460_ESM.zip` → `MOESM460_unzipped/Supplementary_Script_v3/` | analysis pipeline `indrops.py`, `template.yaml`, `barcode_lists/` | **read layout** (which FASTQ holds what, the barcode-read anatomy) |
| `nprot…MOESM459_ESM.zip`, `…458_ESM.zip` (.dwg), `…461–465` (movies) | Hamilton robot methods, chip drawing, videos | not chemistry; not read |

Could not be fetched (these need a manual download):

| Missing | URL | What it would settle |
|---|---|---|
| Cell 2015 **Extended Experimental Procedures** (`NIHMS685995-supplement-1.pdf`) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4441768/bin/NIHMS685995-supplement-1.pdf | the **V1 oligo sequences** from the authors, bead synthesis, RT/lysis buffer, the CEL-Seq library steps and conditions, sequencing setup. Without it every V1-only oligo below is 🟡 (upstream). |
| Cell 2015 Tables S1/S2 (`supplement-11/12.xlsx`), Movies, Fig. images (`supplement-2..10`) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4441768/bin/NIHMS685995-supplement-11.xlsx (and siblings) | run statistics; not chemistry |
| Nature Protocols 2017 **main text** (paywalled) | https://www.nature.com/articles/nprot.2016.154 | all **V2 reaction conditions** (RT in drops, IVT, fragmentation, PE2-N6 RT, PCR cycles), the sequencing recipe and cycle numbers |
| Not attempted: the indrops GitHub repo (source of upstream's "V3") | https://github.com/indrops/indrops | V3 read layout |

---

## 1. What it is, and what it builds on

inDrop puts single cells into ~4 nL droplets 🟢 together with RT/lysis mix and one
**barcoded hydrogel microsphere (BHM)**. Each bead carries many copies of one barcoded
oligo-dT primer, attached through an acrydite and a **photo-cleavable spacer**. After
encapsulation, UV light **releases the primers into the droplet** 🟢. The authors call the
release step critical for efficient RT, citing Figs. 1 and 3F 🟢.
RT happens in the droplets. The droplets are then broken, and the pooled cDNA is processed
"as per CEL-SEQ protocol" 🟢 (Cell, Library preparation), i.e. linear amplification by
**T7 in-vitro transcription**. The acknowledgements thank a colleague for guidance on
"the CEL-SEQ/MARS-SEQ protocol" 🟢.

| | New thing here | Builds on / compare |
|---|---|---|
| 1 | **Combinatorial split-and-pool barcoding of hydrogel beads**: two rounds of primer extension on 384-well plates, 384 × 384 = **147,456** barcodes 🟢 (Cell; the 384 × 384 is 🟡 computed from the plate sheets) | Drop-seq (same issue of *Cell*) builds its barcodes by split-pool phosphoramidite synthesis on hard beads instead: [Drop-seq](../drop-seq-seq-well__10.1016+j.cell.2015.05.002/01_drop-seq-seq-well.md) |
| 2 | **Deformable beads** that pack and flow in step, so almost every droplet gets exactly one bead 🟢 (Cell: "nearly 100% hydrogel droplet occupancy"), and **photo-release** of the primer into solution | Drop-seq / 10x Chromium: Drop-seq primes on the bead surface; 10x dissolves its gel beads |
| 3 | A barcode-1 of **variable length (8–11 nt)**, which staggers the bases that follow so the constant W1 region is not read by every cluster in the same cycle 🟡 (the length set is 🟢 from the plate sheet and the script; the purpose is not stated in what we read) | |
| 4 | Library chemistry: none. V1 is **CEL-Seq / MARS-seq**: T7 promoter on the RT primer, second-strand synthesis, IVT, aRNA fragmentation, 3'-adapter ligation (RLO), second RT. V2 swaps the ligation for **random-hexamer RT** (PE2-N6) | [CEL-Seq](../cel-seq-family__10.1016+j.celrep.2012.08.003/01_cel-seq.md), [MARS-seq](../mars-seq-mars-seq2-0__10.1126+science.1247651/01_mars-seq.md); concepts: [reverse transcription](../ref/concepts/reverse-transcription.md), [RNA adapter ligation](../ref/concepts/small-rna-ligation.md) |

**Variants.** Upstream distinguishes three variants 🟡 (upstream):

- **V1** is the Cell 2015 paper.
- **V2** is the Nature Protocols 2017 paper. The bead oligo is identical to V1, but the
  library primers are different and the library is flipped on the flow cell.
- **V3** uses the same oligos as V2 with a four-read sequencing layout. It has no paper in
  our sources, and upstream itself says its V3 drawing is a guess.

This note covers all three; §5 and §6 give each variant's differences.

## 2. Oligos

### 2.1 Bead synthesis and QC: identical in V1 and V2

🟢 From Nat Protoc Supplementary Table (`MOESM455`, sheet `oligos_except_barcodes`),
written exactly as there. The acrydite primer is printed in triplets with its
modifications spelled out:

```
Acrydate-modified primer  (Acrydite) (PC Spacer) CGA TGA CGT AAT ACG ACT CAC TAT AGG GAT ACC ACC ATG GCT CTT TCC CTA CAC GAC GCT CTT C    HPLC, Trilink
BA19 oligo                BAAAAAAAAAAAAAAAAAAA
PE1* probe                /56-FAM/AGATCGGAAGAGCGTCGTGTAGGGAAAGAG
W1* probe                 /56-FAM/AAGGCGTCACAAGCAATCACTC
BA19 probe                /56-FAM/BAAAAAAAAAAAAAAAAAAA
```

Upstream writes the same primer as `/acrydite/iSpPC/ CGATGACG…` and the bead product as
`/5Acryd/iSpPC/…`. These are the same modifications in IDT notation 🟡 (upstream).

The barcode plates (🟢, `MOESM456` / `MOESM457`, 384 wells each, spaces as in the sheet;
first well shown):

```
P1  W1*-bc1.1-PE1*      AAGGCGTCACAAGCAATCACTC AAACAAAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAG
    general form        AAGGCGTCACAAGCAATCACTC <bc1: 8-11 nt> AGATCGGAAGAGCGTCGTGTAGGGAAAGAG
P2  BA19-N6-bc2.1-W1*   BAAAAAAAAAAAAAAAAAAA NNNNNN AAACAAAC AAGGCGTCACAAGCAATCACTC
    general form        BAAAAAAAAAAAAAAAAAAA NNNNNN <bc2: 8 nt> AAGGCGTCACAAGCAATCACTC
```

### 2.2 V2 library and sequencing primers

🟢 `MOESM455`, written as there. Purification is standard desalting for the library
primers and HPLC for the three sequencing primers.

```
PE2-N6 primer          TCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNNN
PE2 primer             AATGATACGGCGACCACCGAGATCTACACGGTCTCGGCATTCCTGCTGAAC
PE1 primer index 1     CAAGCAGAAGACGGCATACGAGATCGTGATCTCTTTCCCTACACGA     index read: ATCACG
PE1 primer index 2     CAAGCAGAAGACGGCATACGAGATACATCGCTCTTTCCCTACACGA     index read: CGATGT
  ...                  (24 in all, the last:)
PE1 primer index 24    CAAGCAGAAGACGGCATACGAGATGCTACCCTCTTTCCCTACACGA     index read: GGTAGC
  general form         CAAGCAGAAGACGGCATACGAGAT <revcomp of index> CTCTTTCCCTACACGA     46 nt
Custom Read 1 primer   GGCATTCCTGCTGAACCGCTCTTCCGATCT
Custom Index Read primer  AGATCGGAAGAGCGTCGTGTAGGGAAAGAG
Custom Read 2 primer   CTCTTTCCCTACACGACGCTCTTCCGATCT
```

Upstream lists the same nine V2 sequences, character for character. It writes the PE1
PCR primer as `…GAGAT[6-bp sample index]CTCTTTCCCTACACGA` without saying which
orientation the index is in; §3 settles that.

### 2.3 V1 library primers: 🟡 (upstream only)

The authors' list is in the Cell Extended Experimental Procedures, which we could not
fetch. These sequences come from the scg_lib_structs page only:

```
RNA Ligation Oligo (RLO)   /5Phos/ AGATCGGAAGAGCGGTTCAGCAGGAATGCC /3SpC3/
2nd RT primer              GTCTCGGCATTCCTGCTGAAC
PCR enrichment primer 1    AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA
PCR enrichment primer 2    CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAAC
Read 1 sequencing primer   TCTTTCCCTACACGACGCTCTTCCGATCT
Read 2 sequencing primer   CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
```

## 3. How the oligos interlock: 🟡 (computed with `lib/`: `chemdraw.revcomp`, `illumina`)

**Acrydite primer (64 nt)** = `CGATGACG` (an 8-nt leader) + **T7 promoter**
`TAATACGACTCACTATAGGG` (nt 9–28) + `ATACCACCATGG` + `CTCTTTCCCTACACGACGCTCTTC`.

- The last 24 nt equal **PE1[:24]**, where PE1 = `CTCTTTCCCTACACGACGCTCTTCCGATCT` is
  the V2 Custom Read 2 primer.
- **PE1 = `illumina.TRUSEQ_READ1[3:]`**: the TruSeq Read 1 site without its first `ACA`.
- Upstream and the Cell legend name the elements (T7 promoter and sequencing primer). The
  role of the 8-nt leader and of `ATACCACCATGG` is not stated in what we read 🔴. An
  `NcoI` site, `CCATGG`, sits inside `ATACCACCATGG`; any use of it is unstated 🔴.

**Round 1 (P1 plate).** The P1 oligo is `W1*` + bc1 + `PE1*`, and **`PE1*` =
revcomp(PE1)**, which also equals the first 30 nt of `illumina.INDEX2_PRIMER_RC`.

- The bead primer's 3' 24 nt pair with the 3' 24 nt of `PE1*`. Extension copies the rest,
  so the bead strand gains `CGATCT` (completing PE1), then **revcomp(bc1)**, then
  **W1** = `GAGTGATTGCTTGTGACGCCTT` (`W1*` = revcomp(W1), 22 nt).

**Round 2 (P2 plate).** The P2 oligo is `BA19` + `N6` + bc2 + `W1*`. Its 3' `W1*` pairs
with the W1 now on the bead, and extension adds **revcomp(bc2)** + **revcomp(N6)** (the
UMI) + **`T19V`**. That last part is the complement of `BA19`, with B = not-A giving
V = not-T as the anchor.

**The finished bead oligo** (the bead-side strand; the plate strand is removed with NaOH
per upstream 🟡):

```
5'- /acrydite/ /PC spacer/ · CGATGACG · T7 promoter (20) · ATACCACCATGG · PE1 (30) · bc1' (8–11) · W1 (22) · bc2' (8) · UMI (6) · T19 · V -3'
```

- Length 64 + 6 + bc1 + 22 + 8 + 6 + 20 = **134–137 nt**.
- The bead QC in Nat Protoc Supplementary Method 2 expects bands of **135, 102 (or 101)
  and 65 nt** before Exonuclease I clean-up and only 135 nt after 🟢. By our arithmetic
  those are the full oligo, the oligo after round 1 only (64 + 6 + bc1 + 22 = 100–103) and
  the unextended primer (64). This identification is 🟡 (inferred). Note that the
  Supplementary Method itself gives 102 nt in one place and 101 nt in the other 🟢.
- Barcodes are written in **plate orientation** in the sheets. The bead, and every read
  of it, carries their **reverse complement** (bc1 `AAACACGGT` becomes `ACCGTGTTT` on the
  bead).

**Barcode sets.**

- P1: 384 oligos, all with the same `W1*`/`PE1*` arms, bc1 lengths **8, 9, 10 and 11 nt
  × 96 each**.
- P2: 384 oligos, bc2 always 8 nt.
- The pipeline's `gel_barcode1_list.txt` / `gel_barcode2_list.txt` are **identical, in
  order, to the plate barcodes**.
- `indrops.py` reverse-complements each list entry before matching (`barcode =
  rev_comp(line…)`) 🟢. That agrees with the orientation above.

**V2 library primers.**

- **PE2-N6** = `TCGGCATTCCTGCTGAACCGCTCTTCCGATCT` (32 nt, here called the PE2 site) + N6.
  The PE2 site ends in `illumina.STEM_COMPLEMENT` + `T` (`…GCTCTTCCGATCT`), like every
  Illumina read-primer site.
- **PE2 primer (51 nt)** = `illumina.P5` + `GGTC` + PE2 site[:18]. Its 3' 18 nt match the
  5' end of PE2-N6, so it primes on the complement of the RT primer and puts **P5 on the
  cDNA side**.
- **Custom Read 1 primer** = PE2 site[2:] (30 nt). It reads **into the cDNA insert**.
- **PE1 primer index N (46 nt)** = `illumina.P7` + 6-nt index + PE1[:16]. Its 3' 16 nt
  anneal to the `PE1*` that the PE2-N6 RT copies from the aRNA. So **P7 is on the bead /
  barcode side**.
- **For all 24 PE1 primers, the 6 nt inside the primer are the reverse complement of the
  index listed beside it** (index 1: `CGTGAT` in the primer, `ATCACG` listed). The listed
  index is what the sequencer reads (see §6). All 24 primers are 46 nt.
- **Custom Index Read primer** = `PE1*`. **Custom Read 2 primer** = PE1.

**V1 primers (upstream sequences).**

- **RLO** = revcomp(V2 Custom Read 1 primer), all 30 nt. That is, the RLO carries exactly
  the V2 PE2 read site, as its complement.
  - It is a 5'-phosphate donor with a 3' C3 block, so it ligates onto aRNA 3' ends and
    cannot concatenate 🟡 (inferred; see [small-RNA ligation](../ref/concepts/small-rna-ligation.md),
    trick 1, here without pre-adenylation).
- **2nd RT primer (21 nt)**: only its **3' 16 nt pair with the RLO** (with the RLO's
  3'-terminal `GTTCAGCAGGAATGCC`). The 5' `GTCTC` overhangs and adds `GAGAC` to the
  library.
- **PCR primer 2 (47 nt)** = `illumina.P7` + `CGGTCTCGGCATTCCTGCTGAAC`. It ends in the
  2nd RT primer, so **P7 goes on the RLO / cDNA side in V1**.
- **PCR primer 1 (44 nt)** = `illumina.P5` + `TCTTTCCCTACACGA` =
  `illumina.NEBNEXT_UNIVERSAL_PRIMER[:44]`. It primes on the PE1 site, so **P5 goes on
  the bead / barcode side in V1**.
- **Read 1 primer** = `illumina.TRUSEQ_READ1[4:]` (29 nt).
- **Read 2 primer (37 nt)** = revcomp(RLO + `GAGACCG`). It contains the whole V2 Custom
  Read 1 primer as its 3' 30 nt.

So **V1 and V2 use the same two adapter sites, PE1 on the bead and PE2 on the cDNA, but
the flow-cell ends are swapped.** V1 is P5–PE1–barcodes…cDNA–PE2–P7. V2 is
P5–PE2–cDNA…barcodes–PE1–index–P7. V1 has **no sample index**. V2 has a 6-nt i7 index
read with a custom primer.

`lib/` has no constants for the T7 promoter, W1, or this PE2 site (the old Illumina
paired-end "PE 2.0" read-2 sequence, by its look 🟡). The comparisons above used them as
literals.

## 4. Step by step

### 4.1 Bead (BHM) synthesis: same for V1 and V2

Order from upstream's drawing 🟡 (upstream), consistent with the Cell Fig. 2 legend 🟢
("microfluidic preparation of hydrogel microspheres containing a common DNA primer…
combinatorial barcoding") and with the Nat Protoc plate sheets and QC 🟢. Conditions
(enzyme, temperatures) are 🔴: they are in the unfetched Extended Procedures and Nat Protoc
main text.

1. **Make acrylamide hydrogel beads** in a microfluidic device, copolymerising the
   acrydite primer into the gel 🟢 (Fig. 2A–B legend).
2. **Split into the 384 wells of plate P1** and extend the bead primer on the
   `W1*`-bc1-`PE1*` template. The bead strand gains `CGATCT`-bc1'-W1.
3. **Pool, denature with NaOH, wash** away the plate strand 🟡 (upstream).
4. **Split into the 384 wells of plate P2** and extend on `BA19`-N6-bc2-`W1*`. The bead
   strand gains bc2'-UMI-T19V.
5. **Pool, NaOH strip, neutralise and wash.** Nat Protoc adds an **Exonuclease I**
   clean-up, after which only the 135-nt band remains 🟢 (Supplementary Method 2). That it
   works by digesting the incompletely extended single-stranded primers is 🟡 (inferred
   from ExoI's activity on free 3' ends).
6. **QC**: UV-release primers from 1 µL of packed beads (7 min on ice), then 6 % TBE-urea
   gel 🟢 (Supplementary Method 2). FAM probes against PE1, W1 and BA19 🟢 (oligo sheet,
   "Barcoded hydrogel bead QC") check each synthesis stage 🟡 (their exact use is in the
   main text, not fetched).

### 4.2 Library, V1 (Cell 2015)

Steps 1–3 are 🟢 from Cell "Experimental Procedures". Steps 4–8 follow upstream's drawing
🟡 (upstream); their conditions are 🔴 (Extended Methods not fetched).

1. **Encapsulate** in 4-nL drops. Flow rates are 100 µL/h cells, 100 µL/h RT/lysis mix,
   10–20 µL/h BHMs and 90 µL/h carrier oil. The oil is HFE-7500 with 0.75 % EA
   surfactant. Cells are at 50–100 k/mL in 16 % Optiprep. Cells, RT/lysis mix and
   collection tubes are kept on ice.
2. **Release the primers**: 8 min UV (365 nm, ~10 mW/cm²) on ice.
3. **RT in the drops**: 50 °C 2 h, then 70 °C 15 min, then on ice. The bead oligo's
   `T19V` primes on the poly(A) tail, so the first strand is bead oligo + antisense cDNA 🟡
   (inferred from the oligo design; the conditions are 🟢). Then
   demulsify: split into aliquots of 100–3,500 cells, add 0.2× volume of 20 %
   perfluorooctanol / 80 % HFE-7500, spin, store at −20 °C.
4. **Second-strand synthesis** (RNase H / DNA Pol I) gives a double-stranded T7 promoter
   🟡 (upstream; standard CEL-Seq).
5. **IVT** from T7. The aRNA starts `GGGAUACCACCAUGGCUCUUUCC…` (transcription from the
   first G of `GGG`) and runs PE1 · bc1' · W1 · bc2' · UMI · U19 · V · antisense mRNA
   🟡 (computed). Upstream starts the aRNA at `GAUACC…`, dropping the first two Gs. It is
   a drawing simplification, but strictly the +1 is the first G of the promoter's `GGG`.
6. **Fragment the aRNA, then ligate the RLO** to the 3' ends of the fragments. The 3' C3
   block allows one ligation product only 🟡 (upstream).
7. **Second RT** with the 2nd RT primer on the RLO. The cDNA runs from the RLO back to
   the aRNA 5' end, i.e. through UMI, barcodes and PE1 🟡.
8. **PCR** with PCR primers 1 (P5, on PE1) and 2 (P7, on the RT2 / RLO end) 🟡.

### 4.3 Library, V2 (Nature Protocols 2017)

The reaction order is from upstream's V2 drawing 🟡 (upstream), and it agrees with the
oligo sheet's grouping 🟢 ("Library preparation": PE2-N6, PE2, PE1 index 1–24). Every
condition is 🔴 because the main text was not fetched.

1–5 as V1: same beads, encapsulation, UV release, in-drop RT, second strand, IVT.

6. **Random-primed RT on the aRNA with PE2-N6.** The N6 primes anywhere on the antisense
   aRNA. The first-strand cDNA reads PE2 site · sense mRNA · poly(A) · N6 · bc2 · `W1*` ·
   bc1 · `PE1*` · `CCATGGTGGTATC(CC)` 🟡 (computed). This replaces V1's fragmentation +
   RLO ligation + RT2. Whether V2 still fragments the aRNA before this RT is not stated in
   anything we read 🔴. Upstream's V2 drawing goes straight from IVT to PE2-N6 RT.
7. **PCR** with the PE2 primer (P5) and one PE1 index primer (P7 + i7). The PE1 primer
   anneals inside `PE1*`, so the `CCATGG…` tail left over from the T7 leader is trimmed
   off the final product 🟡 (computed).

The protocol prices reagents at about **$600 per 10,000 cells** (about $0.06 per cell),
$400 of it library preparation 🟢 (Supplementary Table 1).

## 5. Final libraries: 🟡 (assembled from the oligos above with `lib/`)

**V1** (top strand = the P5 strand, as upstream; built by simulating steps 5–8):

```
5'- P5 · TCTTTCCCTACACGACGCTCTTCCGATCT (PE1[1:], = Read 1 primer) · bc1' (8–11) · W1 (22) · bc2' (8) · UMI (6) · T19 · V · <antisense mRNA fragment> · RLO (30) · GAGACCG · P7' -3'
```

The upstream final structure has these segments in this order. Its spelled-out 3' end
`…AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG` equals RLO + `GAGACCG` +
`illumina.P7_RC` in our simulation. **Agreement** (both rest on upstream's oligos).

**V2**:

```
5'- P5 · GGTC · TCGGCATTCCTGCTGAACCGCTCTTCCGATCT (PE2 site) · <sense mRNA insert> · A19 · B · N6 (UMI, plate orientation) · bc2 (8) · W1* (22) · bc1 (8–11) · PE1* (30) · <i7, as listed> · P7' -3'
```

- Simulated with a test insert. The upstream V2 final structure matches segment for
  segment, including the 6-bp index between `PE1*` and `ATCTCGTATGCCGTCTTCTGCTTG`
  (= `illumina.P7_RC`). **Agreement.**
- Fixed parts: 29 + 4 + 32 + 20 + 6 + 8 + 22 + bc1 + 30 + 6 + 24 = **181 + bc1 = 189–192
  nt plus the insert** 🟡 (computed).

## 6. Read layout and sequencing

**V2**, from the oligo sheet and the analysis script:

| Read | Primer (🟢 sheet) | Reads (🟡 computed on the V2 construct) |
|---|---|---|
| Read 1 | Custom Read 1 = PE2 site[2:] | **sense mRNA** from the random-primed end, the "bio" read |
| Index 1 (i7) | Custom Index = `PE1*` | the **6-nt index exactly as listed** in the sheet (simulated: index 1 reads `ATCACG`) |
| Read 2 | Custom Read 2 = PE1 | **bc1' · W1 · bc2' · UMI · T…**: the bead orientation, the reverse complement of the plate barcodes |

- `template.yaml` 🟢: `raw_read_bio_fastq: R1.fastq` is the "sequencing read of
  transcript"; `raw_read_meta_fastq: R2.fastq` is the "sequencing read of cell barcode
  and UMI". **Agreement** with the table above and with upstream.
- `indrops.py` 🟢 gives the barcode-read anatomy as bc1 (8, 9, 10 or 11) · W1 (22) · bc2
  (8) · UMI (6) · polyT.
  - It finds W1 at position 8–11, allowing Hamming ≤ 3.
  - It requires ≥ 7 T's after the UMI (Hamming ≤ 3).
  - It rejects pairs whose bio read contains `W1*` ("almost certainly empty library").
  - It slices the UMI as **7 nt** (`name[umi_pos:umi_pos+7]`), although the docstring
    says 6. This is a quirk of the published script 🟢.
  - Internally the barcode read is called "R1" and is fed from the `R2` file. This
    naming is a leftover of V1, where the barcodes were Read 1 🟡 (inferred).
- Cycle numbers and instrument: 🔴, not in what we read. Upstream says the barcode read
  needs "at least 51 cycles" 🟡 (upstream). Arithmetic: the longest bc1 (11) + 22 + 8 + 6
  = 47 nt, plus ≥ 7 T for the script's poly(T) check = **≥ 54 cycles** to pass the filter
  on every bead 🟡 (computed).

**V1**, all 🟡 (upstream):

- Read 1 (TruSeq-type primer `TCTTTCCC…CGATCT`) = bc1' · W1 · bc2' · UMI · T. Upstream
  says 51 cycles.
- Read 2 (`CGGTCTCGG…CGATCT`) = the cDNA, reading the sense strand from the fragment end
  toward the poly(A) 🟡 (computed).
- **No index read.**

**V3**: upstream's guess 🟡 (upstream, explicitly "not confirmed"), on the V2 library.

- Read 1 is the cDNA.
- A read primed with `W1*` reads bc1 (8 cycles).
- The index read is primed with `PE1*` (6 cycles).
- Read 2 is primed with W1 and reads bc2 + UMI (≈21 cycles).
- Our computation on the V2 construct puts both W1-primed reads where upstream says, but
  only for the V2 constructs as assembled above. No author source for V3 was read 🔴.

## 7. Agreements and disagreements with upstream (scg_lib_structs)

- **Agree** 🟢 vs 🟡: all V2 oligos (acrydite primer, BA19, the three FAM probes, PE2-N6,
  PE2, PE1 index form, three custom sequencing primers) and both barcode-plate general
  forms are identical in the Nat Protoc sheets and upstream.
- **Agree**: bc1 has variable length (upstream: "variable length"; sheets: 8–11). There are
  384 bc1 and 384 bc2.
- **Agree**: the read assignment in V2 (cDNA = Read 1, barcodes = Read 2) matches the
  pipeline's `template.yaml`.
- **Adds**: the PE1 index is written in the primer as the **reverse complement** of the
  listed (read) index. Upstream leaves the orientation unspecified.
- **Minor**: upstream's aRNA starts `GAUACC…`. Transcription from this T7 promoter starts
  at the first G of `GGG`, so the aRNA starts `GGGAUACC…` 🟡.
- **Unverifiable here**: every V1-only oligo (RLO, 2nd RT, V1 PCR 1/2, V1 read primers),
  the V1 steps, and the V3 layout. The author source (Cell Extended Methods) is behind
  PMC's download gate.

## 8. Open questions

- 🔴 The V1 oligo sequences from the authors (Cell Extended Experimental Procedures). Are
  they exactly upstream's? In particular, are the RLO's `/5Phos/` and `/3SpC3/`
  modifications as upstream writes them?
- 🔴 All V2 conditions: RT enzyme and buffer in drops, IVT kit and time, whether the aRNA
  is fragmented before PE2-N6 RT, RT enzyme for PE2-N6, PCR cycles, size selection.
- 🔴 Why V2 flipped the library relative to V1 (barcodes moved from Read 1 to Read 2, cDNA
  to Read 1). A plausible reason is cluster calling on the low-diversity W1 / poly(T)
  stretch, which the variable-length bc1 also addresses, but this is unstated.
- 🔴 The role of the 8-nt leader `CGATGACG` and of `ATACCACCATGG` between the promoter and
  PE1.
- 🟡 The read-cycle settings for V2 and the V3 four-read layout. Check against the
  indrops GitHub repository and the Nat Protoc main text.

## 9. How this note was made (tool evaluation)

- `get_sources.py` fetched the upstream page, the PMC full text, and every Nat Protoc
  supplement from Springer. It could not get the Cell supplements.
  - It reports supplement-1, 2, 10, 11 and 12 correctly as `(manual)`.
  - But it saved **supplement-3 (.avi) and 4–9 (.png) as "saved"**, and each of those
    files is a ~21 KB Google reCAPTCHA HTML challenge page, not the media. They are in the
    manifest with md5s as if real.
  - It also left the zip supplements unpacked. `MOESM460` (the analysis scripts and
    barcode lists) had to be unzipped by hand into `MOESM460_unzipped/`.
- `scrape_primers.py --max-hits 80` ranked the oligo sheet first and recognised P5/P7 in
  every PCR primer.
  - It **missed the acrydite primer in the Nat Protoc sheet**, which is written in
    space-separated triplets. It found it only in the upstream page.
  - Its `--find` did find the 64-mer in the sheet.
  - The plate sheets (384 lines each) were cut by `--max-hits`, as intended.
- `scrape_primers.py --find` returned **0 locations for the full P2 oligo**
  `BAAAAAAAAAAAAAAAAAAANNNNNNAAACAAACAAGGCGTCACAAGCAATCACTC`. The oligo is verbatim (with
  spaces) on `MOESM457` line 3, and each half of it is found. It also reported
  `BAAAAAAAAAAAAAAAAAAANNNNNN` at lines that contain no N6, which suggests the N's are
  treated as wildcards. The 🟢 for that oligo rests on a direct `grep`.
