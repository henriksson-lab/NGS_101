# Drop-seq / Seq-Well — 3' mRNA counting on barcoded beads, in droplets or in sealed microwells

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published but not obtained. Relationships marked 🟡 *(computed)* were
> worked out with `lib/`; the model is `tools/drop_seq.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not verification. Nothing it alone asserts is treated as published.
>
> **Where the bases come from.** The available Drop-seq main article gives bead geometry
> but sends exact primers to Extended Experimental Procedures/Table S6, which remains
> behind the PMC download gate 🔴. The inventor-lab [Drop-seq protocol v3.1](https://mccarrolllab.org/wp-content/uploads/2015/05/Online-Dropseq-Protocol-v.-3.1-Dec-2015.pdf)
> by Macosko and Goldman (28 Dec 2015, p. 3 required-primer table) prints the bead,
> TSO, SMART PCR, P5 hybrid and custom Read 1 oligos and maps their names to Macosko et al. 🟢
> The later **Seq-Well Supplementary Table 1** corroborates the same mRNA capture bead
> family 🟢. Assigning the later printed strands to the original Cell experiment remains
> an **inference** until its Table S6 is checked
> 🟡, and every segment it supplies is created `inferred=True` and drawn dotted.

**Drop-seq** — Macosko EZ, Basu A, Satija R, Nemesh J, Shekhar K, Goldman M, Tirosh I,
Bialas AR, Kamitaki N, Martersteck EM, Trombetta JJ, Weitz DA, Sanes JR, Shalek AK,
Regev A, McCarroll SA. "Highly parallel genome-wide expression profiling of individual
cells using nanoliter droplets." *Cell* 161:1202–1214 (2015).
doi:[10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002) · PMID 26000488 ·
PMC4481139 (author manuscript). Data: GEO GSE63473.

**Seq-Well** — Gierahn TM, Wadsworth MH II, Hughes TK, Bryson BD, Butler A, Satija R,
Fortune S, Love JC, Shalek AK. "Seq-Well: portable, low-cost RNA sequencing of single
cells at high throughput." *Nat Methods* 14:395–398 (2017).
doi:[10.1038/nmeth.4179](https://doi.org/10.1038/nmeth.4179) · PMID 28192419 · PMC5376227.

Also listed for this method in `catalogue/scg_lib_structs.tsv`: nanoCAGE / CAGEscan
(Plessy et al., *Nat Methods* 2010, doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470)).
Not fetched for this protocol; nothing here depends on it.

Sources read (cached under `$CHEM_DATA/sources/`, never committed):

| File | What | Used for |
|---|---|---|
| `j.cell.2015.05.002_PMC4481139.html.txt` | Drop-seq full text (PMC author manuscript) | concept, bead synthesis, Experimental Procedures (short form), read layout in Fig. 2C legend |
| `crop-seq__10.1038+nmeth.4177/Online-Dropseq-Protocol-v.-3.1-Dec-2015.pdf.txt` | **Macosko and Goldman inventor-lab Drop-seq protocol v3.1**, 28 Dec 2015 (the scoped Drop-seq route currently fails DNS; identical original cached in the CROP-seq source bundle) | p. 3 maps exact bead, TSO, PCR and custom Read 1 primer names to Macosko et al.; pp. 18–20 give Nextera N70X PCR and read settings. Its opening notice says the optimized protocol does not precisely match the experiments in the defining Cell paper. |
| `nmeth.4179_..._MOESM237_ESM.pdf.txt` | Seq-Well Supplementary Figures + **Supplementary Table 1** | **every oligo** (Seq-Well) |
| `nmeth.4179_..._MOESM238_ESM.pdf.txt` | **Seq-Well Master Protocol v1.0** (Shalek & Love labs, 13 Feb 2017) | **all Seq-Well conditions**: arrays, sealing, lysis, RT, ExoI, WTA, Nextera, sequencing, buffers |
| `upstream_Drop-seq.html.txt` | scg_lib_structs method page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `nmeth.4179_..._MOESM234/235/236/240/242–246_ESM.xlsx.txt` | Seq-Well source data (barcode counts, tSNE, gene lists, GSEA) | not chemistry |

Could not be fetched (🔴 for everything only they contain):

| What | URL | Why it matters |
|---|---|---|
| Drop-seq **Supplemental Information** (Extended Experimental Procedures; contains Table S6, the primers) | https://pmc.ncbi.nlm.nih.gov/articles/instance/4481139/bin/NIHMS687993-supplement-supp_data_1.pdf | the exact defining-paper oligos and reaction conditions — the main text defers details to it; the later inventor-lab protocol is corroboration, not this original table |
| Drop-seq Fig. S1, S6, Tables S2/S3 | `.../NIHMS687993-supplement-Figure_S1.pdf`, `..._Figure_S6.pdf`, `..._supp_data_2.xlsx`, `..._supp_data_3.xlsx` | bead QC; not chemistry |
| Drop-seq Fig. S2–S5, Table S4, Movie S1, DataFile 1 (device .dwg) | `.../NIHMS687993-supplement-Figure_S2.pdf` … `..._supp_data_6.dwg` | saved by the tool, but each file is a ~21 kB reCAPTCHA page, not the file |
| Seq-Well old HTML route | https://pmc.ncbi.nlm.nih.gov/articles/PMC5376227/ (saved file was a reCAPTCHA page) | full-text JATS is now cached as `nmeth.4179_PMC5376227.xml.txt`; it explicitly compares with Drop-seq using the same mRNA capture beads (line 76) |
| Seq-Well MOESM239 / MOESM241 (zips), MOESM247/248 (videos) | Springer ESM | zips truncated at 50 MB; videos not chemistry |

---

## 1. What it is

Both methods put **one cell and one barcoded bead** in a small compartment, lyse the cell
there, and let its poly(A) mRNA hybridise to the bead's oligo-dT. Everything after that is
done **in bulk on the pooled beads** ("STAMPs", single-cell transcriptomes attached to
microparticles): reverse transcription from the bead-bound primer, template switching,
exonuclease I, whole-transcriptome PCR, then a Nextera XT library that keeps only the
**3' end** of each cDNA. 🟢 (Drop-seq main text, Fig. 2A legend and Experimental
Procedures; Seq-Well Master Protocol)

The two differ only in **the compartment**:

| | Drop-seq | Seq-Well |
|---|---|---|
| Compartment | ~1 nL aqueous droplet in oil, made on a co-flow microfluidic device 🟢 | PDMS array of microwells; beads and cells pipetted dropwise onto the array and allowed to settle (array rocked) 🟢, then **sealed with a plasma-treated polycarbonate membrane** (37 °C 30 min in a clamp) 🟢 |
| Lysis | lysis buffer flows in with the beads, so lysis starts at droplet formation 🟢 (buffer composition 🔴, in the unfetched supplement) | after sealing: 5 M guanidine thiocyanate, 1 mM EDTA, 0.5 % Sarkosyl, 1 % β-ME, 20 min, with the membrane still on 🟢 (that the buffer reaches the cells through the 0.01 µm membrane is our inference 🟡) |
| Hybridisation | inside the droplet; droplets then broken with perfluorooctanol in 30 mL 6× SSC — the large volume slows re-hybridisation to the wrong bead 🟢 | 2 M NaCl, 3 mM MgCl₂, 0.5 % Tween-20, 40 min, still sealed 🟢 |
| Bead recovery | from the broken emulsion 🟢 | array inverted onto lifter slips and spun 1000 × g 5 min 🟢 |
| Beads | made by ChemGenes (reverse-direction phosphoramidite synthesis, split-and-pool) 🟢 | the same beads, ChemGenes cat. no. MACOSKO-2011-10 🟢 |

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Split-and-pool synthesis of the cell barcode directly on the bead**: 12 rounds of split-into-four / add one base / pool → 4¹² = 16,777,216 barcodes, then 8 degenerate rounds for the UMI, then T30. Each bead carries > 10⁸ primers with one cell barcode 🟢 | inDrop (hydrogel beads, barcodes ligated), 10x Gel Beads, CEL-seq-style UMIs |
| 2 | **Bulk RT after cell isolation**: the mRNA only hybridises in the compartment; the covalent capture (RT) happens after all beads are pooled 🟢 | Seq-Well, Microwell-seq, BD Rhapsody |
| 3 | **3'-end selection at the Nextera step** by a P5 primer that only matches the bead end ("custom primers that enabled the specific amplification of only the 3' ends", Drop-seq main text 🟢) | 10x 3' v1–v3, Seq-Well, most bead-based 3' methods |

Builds on: [reverse transcription](../ref/concepts/reverse-transcription.md) from a
solid-phase oligo-dT, [template switching](../ref/concepts/template-switching.md) with the
SMART handle (cited as Zhu et al. 2001 in Fig. 2A), and
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) (Nextera XT) of the amplified cDNA.

## 2. Oligos

### Seq-Well — 🟢 verbatim from Supplementary Table 1 (`MOESM237`)

IDT notation as given: `rG` = ribo-G, `*` = phosphorothioate bond. `J` = split-and-pool
base (same on every primer of a bead), `N` = degenerate base (different on each primer).

```
1. Barcoded Bead Sequence          5'–Bead–Linker-TTTTTTTAAGCAGTGGTATCAACGCAGAGTACJJJJJJJJJJJJNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3'
2. Template Switching Oligo        AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
3. SMART PCR Primer                AAGCAGTGGTATCAACGCAGAGT
4. New P5-SMART PCR Hybrid Oligo   AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C
5. Custom Read 1 Primer            GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
```

(The PDF wraps 1, 4 and 5 across lines and has a stray space in 4 and 5; joined here. 12
`J`, 8 `N` and 30 `T` counted in the source. 🟢) Not in the table, from the Master
Protocol's shopping list: **Nextera N70X Oligo** (Illumina, i.e. the Nextera XT i7 index
primers). 🟢 Their sequence is not given; the standard form is `nextera.n7xx_primer()` =
`CAAGCAGAAGACGGCATACGAGAT <i7> GTCTCGTGGGCTCGG`. 🟡

### Drop-seq's own oligo table — 🔴 not obtained

The paper's Table S6 ("custom primers that enabled the specific amplification of only the
3' ends" 🟢, main text) sits in the Supplemental Information behind the PMC download gate.
The inventor-lab v3.1 required-primer table **does** print the SeqB bead constant `AC`,
TSO, P5 hybrid with two terminal phosphorothioates, and custom Read 1 primer, using
Macosko et al. names. It is a later original from the inventors, not the defining
paper's Table S6. Therefore, for the **original experiment**, still unresolved:

- 🔴 whether the defining 2015 Table S6 used the v3.1 SeqB bead constant `AC` or another revision;
- 🔴 whether the defining Table S6's P5-SMART hybrid had the v3.1 two 3′ phosphorothioates;
- 🔴 the exact Drop-seq Read 1, index and Read 2 primers.

An earlier revision of this note reproduced a second bead batch ("seqA", handle + `ACGT`)
and a matching Read 1 primer, plus pre-amp and library primers, from the scg_lib_structs
page. None of that has a primary source, so none of it is asserted any more, and
`tools/drop_seq.py` no longer models a second bead batch. 🔴

### How the oligos interlock — 🟡 (computed)

One 23-nt handle, **`rt.SMART_HANDLE` = `AAGCAGTGGTATCAACGCAGAGT`**, is in every
non-Illumina oligo:

- **Bead** = `TTTTTTT` (linker) + **handle** + `AC` + J12 + N8 + T30 🟢 (Seq-Well
  Supplementary Table 1).
- **TSO** = **handle** + `GAAT` + `rGrGrG` (27 DNA + 3 RNA = 30 nt).
- **SMART PCR primer / ISPCR** *is* the handle (23 nt).
- **P5-SMART hybrid** = **`illumina.P5`** (29) + `GCCTGTCCGCGG` (12) + **handle** (23) +
  `AC` = 66 nt.
- **Custom Read 1 primer** = the P5 hybrid minus P5, exactly (37 nt). It ends on the last
  constant bead base, so read 1 starts on the first cell-barcode base. 🟡 *(computed)*
- Read 2 and the i7 index read use the stock Nextera primers `nextera.READ2_PRIMER` and
  `nextera.INDEX1_PRIMER`, because the only Nextera-side oligo the Master Protocol names
  is the "Nextera N70X Oligo" and no custom Read 2 or index primer is listed 🟡
  *(derived; the lib constants come from Illumina's adapter document, not from any
  protocol-specific claim)*.

The consequences:

1. After RT and template switching, **both ends of every cDNA carry the handle** — the
   bead end as `handle-AC-…`, the TSO end as `handle-GAATGGG-…` — so the single SMART
   PCR primer amplifies the whole cDNA (whole-transcriptome amplification, WTA).
2. The two ends differ **right after the handle**: `AC` at the bead end, `GA` at the TSO
   end. The P5 hybrid ends in `…GCAGAGT*A*C`, so its 3'-terminal two bases match only the
   **bead end**. On the TSO end they are mismatched (`AC` against a template that calls for
   `GA`) and the primer does not extend. 🟡 (computed; the paper's wording is only "custom
   primers that enabled the specific amplification of only the 3' ends" 🟢)
3. In the Nextera PCR the other primer is a plain **N7xx** (s7) 🟢 (Master Protocol
   shopping list). No s5 primer is added, so of the bead-end fragments only those whose
   Tn5 end is s7 amplify exponentially — the 3'-most fragment of each cDNA. 🟡 *(computed
   from the two named primers)*. Internal Tn5–Tn5 fragments with s7 at both ends could in
   principle be primed by N7xx alone; no source we hold discusses them (presumably
   suppressed by the inverted-repeat panhandle, unverified 🟡).
   `nextera.amplifiable` models only the standard N5xx + N7xx case and does not apply
   directly.
4. The 12-nt `GCCTGTCCGCGG` spacer is 83 % GC and raises the custom Read 1 primer to a
   computed Tm of ~71 °C, against ~59 °C for the bare handle (`chemdraw.tm`, 0.5 µM primer,
   50 mM Na⁺). 🟡 (computed). Why the spacer was added is not stated 🔴.

## 3. Step by step

### Seq-Well — 🟢 from the Master Protocol (MOESM238), per array

1. **Membranes**: polycarbonate, 0.01 µm, 22 × 66 mm, oxygen-plasma treated 5–7 min, then
   hydrated in PBS (use within 48 h).
2. **Beads into wells**: ~110,000 beads in 200 µL bead loading buffer (100 mM sodium
   carbonate, 10 % BSA, pH 10) dropwise; wash off surface beads. (The buffer guide's own
   recipe, 2.5 mL of 100 mg/mL BSA in 25 mL, gives 10 mg/mL = 1 % BSA, not the 10 % it
   states 🟡 (computed); which is meant is unresolved.)
3. **Cells into wells**: 10,000 cells in 200 µL RPMI + 10 % FBS, rock 5 min; wash 4× PBS
   (FBS prevents sealing), leave in RPMI without FBS. (Optional variant: fix with
   CellCover, stain, image, *then* load beads.)
4. **Seal**: lay the plasma-treated membrane on the array, clamp, 37 °C 30 min.
5. **Lyse**: pre-lysis buffer until the top slide lifts, then complete lysis buffer (5 M
   GuSCN, 1 mM EDTA, 0.5 % Sarkosyl, 1 % β-ME), 20 min rotating. That the sealed membrane confines the mRNA to its well so
   it hybridises to that well's bead is the design rationale, not stated in these steps 🟡.
6. **Hybridise**: 2 M NaCl, 3 mM MgCl₂, 0.5 % Tween-20 in PBS, 40 min.
7. **Recover beads**: wash buffer (2 M NaCl, 3 mM MgCl₂, 20 mM Tris pH 8), peel membrane,
   invert array onto lifter slips, spin 1000 × g 5 min.
8. **RT with template switching** (200 µL): Maxima H⁻ RT in 1× Maxima buffer, **4 % Ficoll
   PM-400** (40 µL of 20 %), 1 mM dNTPs (20 µL of 10 mM), RNase inhibitor, **2.5 µM TSO** (5 µL of 100 µM); final
   concentrations 🟡 (computed from the 200 µL mix).
   30 min room temperature, then 52 °C 90 min, rotating. First strand runs from the bead's
   T30 to the 5' end of the mRNA; MMLV adds `CCC`, the TSO's `rGrGrG` pairs with it and
   the RT copies the TSO, putting the handle (as complement) on the 3' end of the
   bead-bound cDNA. 🟡 (the standard [template-switching](../ref/concepts/template-switching.md)
   mechanism; the protocol gives only reagents). Washes: TE-TW, TE-SDS, TE-TW ×2.
9. **Exonuclease I**, 200 µL (1× ExoI buffer, 10 µL ExoI), 37 °C 5 min on a heat block + 45 min rotating — removes bead
   primers that captured nothing (single-stranded 3' ends). 🟢 (Drop-seq Experimental
   Procedures: "to remove unextended primers")
10. **WTA PCR** on beads, 1,500–2,000 beads per 50 µL, KAPA HiFi HotStart ReadyMix,
    **0.8 µM SMART PCR primer** (0.4 µL of 100 µM; final concentration 🟡 computed): 95 °C 3 min; 4 × (98 °C 20 s, 65 °C
    45 s, 72 °C 3 min); 9–12 × (98 °C 20 s, 67 °C 20 s, 72 °C 3 min); 72 °C 5 min. Total
    13 cycles for cell lines, 16 for primary cells. Pool ~14–16,000 beads' worth, AMPure
    0.6×; expect 1–2 kb.
11. **Tagment** 750 pg cDNA with Nextera XT (TD buffer + ATM), 55 °C 5 min; neutralise
    (NT buffer) 5 min RT.
12. **Nextera PCR**: 15 µL Nextera PCR mix + **0.2 µM New-P5-SMART PCR hybrid** + 0.2 µM
    one Nextera **N70X** (1 µL of 10 µM each in 50 µL; final concentrations 🟡 computed): 95 °C 30 s; 12 × (95 °C 10 s,
    55 °C 30 s, 72 °C 30 s); 72 °C 5 min. AMPure 0.6×; expect 650–750 bp (420–800 bp
    sequenced fine). No separate 72 °C gap fill-in before denaturation is listed;
    tagmentation's 9-nt gap is presumably filled in the first extension. 🟡
13. **Sequence** (§5).

### Drop-seq — differences

From the main text only 🟢: droplets ~1 nL, cells and beads co-flowed at equal rates;
droplets broken with perfluorooctanol in 30 mL 6× SSC; beads washed; RT in bulk;
exonuclease I; beads counted, aliquoted into PCR tubes, PCR-amplified; cDNA quantified on
a BioAnalyzer; Nextera XT with "custom primers" (Table S6); sequenced on a NextSeq 500.
Lysis buffer, RT mix, cycle numbers, beads per PCR and Nextera input are in the unfetched
Extended Experimental Procedures 🔴. Seq-Well uses the same bead and oligos (§2), so the
molecular steps 8–12 above are the same chemistry; the published Drop-seq conditions may
differ in detail. 🟡

## 4. Final library — 🟡 (assembled from the Seq-Well oligos above)

Top strand, 5'→3' (read 1 reads this strand's sequence; read 2 reads its reverse complement — see §5):

```
5'- P5 · GCCTGTCCGCGG · SMART handle · AC · <cbc 12> · <umi 8> · T30 · <cDNA insert> · ME' · s7' · <i7'> · P7' -3'
```

Segment by segment:

| Segment | Sequence | nt |
|---|---|---|
| `illumina.P5` | `AATGATACGGCGACCACCGAGATCTACAC` | 29 |
| spacer | `GCCTGTCCGCGG` | 12 |
| `rt.SMART_HANDLE` | `AAGCAGTGGTATCAACGCAGAGT` | 23 |
| bead constant end | `AC` | 2 |
| cell barcode | `J` × 12 | 12 |
| UMI | `N` × 8 | 8 |
| oligo-dT | T × 30 (reads as poly(T)) | 30 |
| cDNA insert | from the poly(A) site toward the 5' end of the mRNA, antisense on this strand | variable |
| `nextera.ME_RC` | `CTGTCTCTTATACACATCT` | 19 |
| `nextera.S7_RC` | `CCGAGCCCACGAGAC` | 15 |
| i7 (as reverse complement of the N7xx index) | `<i7'>` | 8 |
| `illumina.P7_RC` | `ATCTCGTATGCCGTCTTCTGCTTG` | 24 |

The 5' end is exactly `illumina.P5` + the Custom Read 1 primer, and the 3' end exactly
`nextera.ME_RC` + `nextera.S7_RC` + 8 index bases + `illumina.P7_RC`. 🟡 *(computed from
the Seq-Well oligo table and the `lib/` constants; the whole construct is an inference for
Drop-seq, which publishes no library structure.)*

The library has **no i5 index**: the P5 end is the custom hybrid primer, identical for all
samples. Samples are multiplexed on i7 only.

## 5. Read layout / sequencing

🟢 Seq-Well Master Protocol: **Read 1 20 bp** with the **Custom Read 1 primer**, **Read 2
50 bp**, **Index 1 8 bp** (only when multiplexing), MiSeq or NextSeq 500; NextSeq loading
at 2.2 pM. 🟢 Drop-seq Fig. 2C legend: the first read gives the cell barcode and UMI; the
second, paired read gives cDNA, typically 50 bp.

| Read | Primer | Reads | Length |
|---|---|---|---|
| Read 1 | custom: spacer + handle + `AC` 🟢 (Seq-Well Table 1), 🟡 for Drop-seq | 12 cell barcode + 8 UMI 🟡 (computed from the primer ending on the last constant base) | 20 🟢 |
| Index 1 (i7) | `nextera.INDEX1_PRIMER` 🟡 *(derived: no custom index primer is named, so the stock Nextera one is implied)* | i7 | 8 🟢 |
| Read 2 | `nextera.READ2_PRIMER` 🟡 *(derived, as above)* | cDNA, sense strand of the mRNA, starting at the Tn5 cut | 50 🟢 |

A read 1 longer than 20 would run into T30. 🟡 (computed)

## 6. Open questions

- 🔴 The original Drop-seq Table S6/Extended Procedures: whether its bead, TSO and P5
  hybrid are base-for-base the later inventor-lab v3.1 primer table and Seq-Well oligos,
  and which original RT, PCR, Nextera and sample-index settings were used. Obtain
  `NIHMS687993-supplement-supp_data_1.pdf` for the direct comparison. Until then the
  model assigns later source-positive oligos to the defining experiment by inference.
- 🔴 The purpose of the 7-T linker and of the `GCCTGTCCGCGG` spacer; the purpose of `GAAT`
  between the handle and `rGrGrG` in the TSO (computed effect: it makes the TSO end
  differ from the bead end right after the handle, §2).
- 🔴 Whether Drop-seq used more than one bead design. The paper describes one bead; no
  primary source we hold mentions a second batch, so none is modelled.
- 🟡 Read 2 and index primers: inventor-lab v3.1 gives 20-base Read 1, 50-base Read 2
  and optional eight-base Index 1, and names the Nextera N70X oligo without assigning a
  particular index member. Exact original Cell run settings still require Table S6.
- 🔴 Seq-Well Online Methods (main paper) were not obtained; array geometry (number and
  size of wells) and the array functionalisation chemistry beyond Appendix A were not
  checked.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` got the Drop-seq PMC full text and all Springer supplements of
Seq-Well, but the PMC supplements of Drop-seq and the Seq-Well PMC page came back as
reCAPTCHA pages — several were saved and listed as fetched although they are ~21 kB HTML.
`tools/scrape_primers.py` found Seq-Well Table 1 in `MOESM237`, but split the bead oligo
at its line wraps (it reported the 23-nt "linker" and T30 separately, missing `J`/`N`),
and cut the P5 hybrid at `GCAGAGT*`, losing the phosphorothioate `*A*C` tail; `--find`
does not match across the `*`. It found nothing in the Drop-seq full text, correctly: the
sequences are only in the unfetched supplement. The Master Protocol has no sequences;
its conditions had to be read.
