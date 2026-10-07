# SCRB-seq and mcSCRB-seq — plate-based 3' tag scRNA-seq with an in-RT well barcode

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`). Claims that rest only on
> the upstream scg_lib_structs page are 🟡 (secondary source).

**SCRB-seq** (single-cell RNA barcoding and sequencing) — Soumillon M, Cacchiarelli D,
Semrau S, van Oudenaarden A, Mikkelsen TS. "Characterization of directed differentiation
by high-throughput single-cell RNA-Seq." *bioRxiv* 2014 (posted 5 Mar 2014).
doi:[10.1101/003236](https://doi.org/10.1101/003236). Data: GEO GSE53638. Never
published in a journal as far as the sources show.

**mcSCRB-seq** (molecular-crowding SCRB-seq) — Bagnoli JW, Ziegenhain C, Janjic A,
Wange LE, Vieth B, Parekh S, Geuder J, Hellmann I, Enard W. preprint "mcSCRB-seq:
sensitive and powerful single-cell RNA sequencing", *bioRxiv* 2017,
doi:[10.1101/188367](https://doi.org/10.1101/188367); published as "Sensitive and
powerful single-cell RNA sequencing using mcSCRB-seq", *Nat Commun* 9:2937 (2018),
doi:[10.1038/s41467-018-05347-6](https://doi.org/10.1038/s41467-018-05347-6). Step-by-step
protocol: protocols.io doi:[10.17504/protocols.io.nrkdd4w](https://doi.org/10.17504/protocols.io.nrkdd4w)
(V.1). Data: GEO GSE103568.

Also in the catalogue for this protocol: nanoCAGE / CAGEscan (Plessy et al., *Nat Methods*
2010, doi:10.1038/nmeth.1470), cited by upstream for "semi-suppressive PCR". An
intermediate SCRB-seq version (Ziegenhain et al. 2017, *Mol Cell*, KAPA HiFi instead of
Advantage 2) is the "SCRB-seq" baseline of the mcSCRB-seq paper; it was not fetched.

Sources read (`tools/get_sources.py` fetched the first three; the rest were fetched by
hand into the same directory, `$CHEM_DATA/sources/scrb-seq-mcscrb-seq__10.1101+003236/`,
and added to its `MANIFEST.tsv`; never committed):

| File | What | Used for |
|---|---|---|
| `003236_v1.full.pdf` | SCRB-seq preprint | overview only; no methods in the main text |
| `003236-1.pdf` | SCRB-seq Supplementary Information (bioRxiv adjunct file) | **all SCRB-seq methods, every oligo, Fig. S1 (strand diagram), Table S2 (384 well barcodes × 2 sets), read lengths** |
| `003236-2..5.xlsx` | Supplementary Data 1–4 (PCA weights, GSEA) | not chemistry |
| `188367_v1.full.pdf` | mcSCRB-seq preprint | methods (identical in substance to the journal version) |
| `s41467-018-05347-6.pdf` | mcSCRB-seq, *Nat Commun* 2018 | methods, TSO comparison |
| `41467_2018_5347_MOESM1_ESM.pdf` | mcSCRB-seq Supplementary Information | Supp. Table 2 (SCRB-seq vs mcSCRB-seq differences), Supp. Fig. 2 (blocked vs unblocked TSO) |
| `41467_2018_5347_MOESM2_ESM.pdf` | peer-review file | not chemistry |
| `protocolsio_mcscrb-seq-nrkdd4w.pdf` | protocols.io mcSCRB-seq V.1 (PDF export) | **mcSCRB-seq oligo table (incl. unblocked TSO, N7xx), volumes, cycling, read set-up** |
| `upstream_SCRB-seq.html` | scg_lib_structs method page | second account; checked in §6 |

Not obtained 🔴:

| What | URL |
|---|---|
| mcSCRB-seq cell-barcode list (`mcSCRBseq_oligodT.txt`, 44 KB, a protocols.io attachment) | https://www.protocols.io/view/mcscrb-seq-protocol-nrkdd4w (attachment, needs the web page / login) |
| Ziegenhain et al. 2017 *Mol Cell* SCRB-seq methods (the intermediate version) | https://doi.org/10.1016/j.molcel.2017.01.023 |
| nanoCAGE paper (only cited for semi-suppressive PCR) | https://doi.org/10.1038/nmeth.1470 |

---

## 1. What it is

A **plate-based, early-barcoding 3' tag** method. One FACS-sorted cell per well is reverse
transcribed with an oligo-dT primer that carries a **6-nt well barcode and a 10-nt UMI**
behind a TruSeq Read 1 handle; a template-switching oligo adds a handle at the other end.
Because every cDNA is barcoded during RT, a whole plate (384 wells in SCRB-seq, 96 or 384 in
mcSCRB-seq) is **pooled before amplification**, amplified by single-primer PCR, and turned
into **one Nextera XT library per plate** in which a custom P5 primer selects only the 3'
end of each cDNA. 🟢 (SCRB-seq Supplementary Information; mcSCRB-seq methods)

It is SMART-seq-style chemistry with the handle swapped for the TruSeq Read 1 sequence —
see [template switching](../ref/concepts/template-switching.md) and
[reverse transcription](../ref/concepts/reverse-transcription.md) — followed by a
**half-Nextera** library: Tn5 supplies only the far end
([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Barcode + UMI on the oligo-dT**, pool before PCR | CEL-seq, MARS-seq, Drop-seq, 10x 3' |
| 2 | **Single-primer cDNA PCR** off a 22-nt handle shared by oligo-dT and TSO (biotinylated) | SMART-seq2 single ISPCR primer (`rt.SMART_HANDLE`) |
| 3 | **3' selection inside Nextera XT**: i5 primer replaced by a P5 + TruSeq Read 1 primer with a 3'-phosphorothioate tail | Drop-seq (custom P5 "P5-TSO hybrid"), Seq-Well, mcSCRB-seq, many 3' plate methods |
| 4 (mcSCRB) | **7.5 % PEG 8000 in RT** (molecular crowding), oligo-dT in the lysis buffer, Terra polymerase, unblocked TSO, bead pooling | crowding agents are standard in ligation; new for template-switching RT here 🟡 |

## 2. Oligos

### SCRB-seq — 🟢 verbatim from the Supplementary Information (`003236-1.pdf`)

Notation as given there: `/5Biosg/` = 5' biotin, `iC` = iso-dC, `iG` = iso-dG,
`rG` = RNA G, `*` = phosphorothioate bond, `[BC6]` = 6-nt well barcode (Table S2),
`N 10` = UMI.

```
E5V6NEXT     5'-iCiGiCACACTCTTTCCCTACACGACGCrGrGrG-3'                         (Eurogentec, 1 pmol)
E3V6NEXT     5'-/5Biosg/ACACTCTTTCCCTACACGACGCTCTTCCGATCT[BC6]N10T30VN-3'      (IDT, 1 pmol)
SINGV6       5'-/5Biosg/ACACTCTTTCCCTACACGACGC-3'                              (IDT, 10 pmol)
P5NEXTPT5    5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG*A*T*C*T*-3'   (IDT, 5 µM)

E3V6NEXT, bulk version (no well barcode):
             5'-/5Biosg/ACACTCTTTCCCTACACGACGCTCTTCCGATCTN16T30VN-3'
```

(The PDF typesets the counts as `N 10 T 30 VN` and `N 16 T 30 VN`; the spaces are
subscript formatting, not bases.) The i7 side is the stock Nextera XT N7xx primer; its
sequence is not printed in the SCRB-seq supplement but is drawn in Fig. S1. Table S2 lists
**two sets of 384 6-nt well barcodes** (set 2 is headed "D3"; set 1 is headed "D2 and D2", read here as D1/D2 🟡 — evidently a typo in the table header).

### mcSCRB-seq — 🟢 verbatim from protocols.io V.1 (`protocolsio_mcscrb-seq-nrkdd4w.pdf`)

As printed in its oligo table (vendor, purification, stock):

```
barcoded oligo-dT (E3V6NEXT)   BiotinACACTCTTTCCCTACACGACGCTCTTCCGATCT[BC6][UMI10][T30]VN   IDT TruGrade, 2 µM
TSO unblocked (E5V6NEXT)       ACACTCTTTCCCTACACGACGCrGrGrG                                  IDT HPLC, 100 µM
PreAmp (SINGV6)                BiotinACACTCTTTCCCTACACGACGC                                  IDT desalted, 10 µM
3' enrichment primer (P5NEXTPT5)  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG*A*T*C*T   IDT HPLC, 5 µM
i7 Index Primer (N7XX)         CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG                    IDT TruGrade, 5 µM
```

The mcSCRB-seq cell barcodes are in an attachment that was not obtained 🔴; whether they
are one of the SCRB-seq Table S2 sets is unknown. The Nat Commun methods name the TSO only
as "template-switching oligo (IDT)"; the protocols.io table is the only place the
**unblocked** sequence appears.

### How the oligos interlock — 🟡 (computed with `lib/`)

- **E3V6NEXT handle = `illumina.TRUSEQ_READ1`** exactly (33 nt). Full E3V6NEXT:
  33 + 6 + 10 + 30 + 2 = **81 nt** (the bulk N16 version is also 81 nt).
- **SINGV6 = `TRUSEQ_READ1[:22]`** (`ACACTCTTTCCCTACACGACGC`, 22 nt, 5'-biotin).
- **Both TSOs carry the same 22 nt.** Unblocked = SINGV6 sequence + `rGrGrG` (25 nt).
  Blocked = `iCiGiC` + the same 22 nt + `rGrGrG`; read as plain bases that is
  `CGC` + 22 + `GGG` (28 nt), which is how upstream draws it. So after template switching
  both cDNA ends carry the 22-nt handle and **one primer (SINGV6) amplifies the cDNA**.
- The oligo-dT end carries 11 more handle bases, **`TCTTCCGATCT`** (`TRUSEQ_READ1[22:]`),
  that the TSO end lacks. That 11-nt difference is what later separates 3' ends from 5'
  ends (step 7 below).
- **P5NEXTPT5 = `illumina.P5` + `TRUSEQ_READ1[4:]`** (58 nt): P5 ends in `ACAC` and Read 1
  starts with `ACAC`, so the two overlap by 4 nt. It is base-for-base identical to
  `illumina.NEBNEXT_UNIVERSAL_PRIMER` / `TRUSEQ_P5_FULL`; what is special is the four
  phosphorothioates on the last 5 bases.
- **N7XX = `illumina.P7` + i7 (8 nt) + `nextera.S7`** (47 nt) — the standard Nextera i7
  primer.
- Fig. S1's full-length cDNA ends, on the top strand, with
  `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` = revcomp(TRUSEQ_READ1) =
  `illumina.INDEX2_PRIMER_RC` — i.e. the oligo-dT handle, copied.
- Well barcodes (Table S2, parsed from the PDF text): set 1 and set 2 are each **384
  distinct 6-mers with minimum pairwise Hamming distance 2**, and the two sets share no
  barcode. Distance 2 detects but cannot correct a single error; the analysis accordingly
  required exact barcode matches (Supplementary Information).

## 3. Step by step — SCRB-seq (Soumillon 2014)

🟢 from the Supplementary Information unless marked.

1. **Sort** single cells (Hoechst; optionally LipidTOX lipid gates) with a FACSAria II into
   384-well plates holding 5 µL 1:500 Phusion HF buffer; RNA stabilised beforehand in
   RNAprotect + RNaseOUT. Seal, spin, dry ice, −80 °C.
2. *(D3 only)* **Proteinase K** 200 µg/mL, 50 °C 15 min sealed, then 95 °C 10 min unsealed —
   kills the protease and **desiccates** the well to shrink the RT volume (2 µL per the
   mcSCRB-seq paper's description of SCRB-seq 🟢).
3. Add **ERCC** spike-ins (1 µL of 1:10⁷ for D1/D2, 1:10⁶ for D3).
4. **RT with template switching**: SmartScribe (D1/D2) or Maxima H Minus (D3), E3V6NEXT and
   E5V6NEXT 1 pmol each. Temperature/time not given in the supplement 🔴 (42 °C 90 min in
   the later SCRB-seq description, Bagnoli 2018 🟢 for that version). First strand 🟡
   (as drawn in Fig. S1, written 5'→3'):
   `5'-biotin · TruSeq R1 (33) · BC6 · UMI10 · T30VN · cDNA · CCC-3'`; the TSO's `rGrGrG`
   pairs with the CCC and the RT copies the TSO, adding the complement of the 22-nt handle
   at the 3' end. The iso-dC/iso-dG at the TSO 5' end cannot be copied by the RT, so
   copying stops there (purpose: block TSO concatemers, Kapteyn et al. 2010, ref. 5 of the
   preprint) 🟡.
5. **Pool** all 384 wells, one Zymo DNA Clean & Concentrator-5 column.
6. **Exonuclease I** (NEB) to remove unused oligo-dT, then **single-primer PCR** with
   SINGV6 (10 pmol) and Advantage 2. Cycle numbers and cycling not given 🔴. Because the
   amplicon carries the same 22 nt inverted at both ends, short products can fold into
   panhandles and amplify poorly — the "semi-suppressive PCR" upstream points to (nanoCAGE)
   🟡. AMPure XP 0.6×, Qubit.
7. **Nextera XT** on the full-length cDNA, per the kit, **except the i5 primer is replaced
   by P5NEXTPT5** (5 µM); i7 = kit N7xx. Only fragments with the oligo-dT end
   (full TruSeq Read 1) at one side and an s7 Tn5 end at the other are exponentially
   amplified 🟡 *(computed)*:
   - the TSO end matches P5NEXTPT5 only over its first 22 handle bases; the primer's last
     11 nt (`TCTTCCGATCT`) have nothing to pair with, so it does not prime there;
   - the 3' phosphorothioates presumably stop a proofreading polymerase from chewing back
     that mismatched 3' tail, which would otherwise let it prime on TSO ends 🟡 (reason not
     stated by the authors 🔴);
   - with no i5 primer, s5 ends are not amplified.
8. AMPure XP 0.6×, **size select 300–800 bp** on a 2 % E-Gel EX, QIAquick gel extraction, Qubit.

## 4. Step by step — mcSCRB-seq (Bagnoli 2018)

🟢 from the Nat Commun methods and protocols.io V.1; final concentrations in brackets are
🟡 *(computed)* from the protocols.io recipes and agree with the paper's Supplementary
Table 2 (0.2 µM oligo-dT, 2 µM TSO, 7.5 % PEG, 20 U RT).

1. **Lysis plates**: 4 µL lysis buffer (1:500 Phusion HF buffer, Proteinase K
   1.25 µg/µL in the 4 µL buffer per protocols.io [1.0 µg/µL in 5 µL 🟡; the paper quotes
   1.25 µg/µL for the whole 5 µL lysis buffer]) + 1 µL 2 µM barcoded
   E3V6NEXT per well [0.4 µM in 5 µL]. Plates with primer
   can be stored at −20 °C.
2. **Sort** one cell per well (Sony SH800, "3 drops" purity) from PBS (no RNAprotect).
   Seal (aluminium), spin, dry ice, −80 °C (≤ 6 months). ERCCs optional (0.1 µL of 1:80,000 Mix 1).
3. **Proteinase K** 50 °C 10 min, **heat-kill** 80 °C 10 min — no desiccation.
4. **RT mix**, 5 µL per well: Maxima H− [20 U], Maxima buffer [1×], dNTPs [1 mM each],
   **unblocked TSO** [2 µM], **PEG 8000** [7.5 %]. 42 °C 90 min. The RT mix is viscous —
   mix carefully, calibrate dispensers.
5. **Pool** the plate and bind to **SPRI-type beads** in 30 % PEG / 2 M NaCl (1:1); 2 × 80 %
   ethanol, elute 17 µL.
6. **Exonuclease I** (2 µL 10× buffer, 1 µL 20 U/µL): 37 °C 20 min, 80 °C 10 min. (The
   protocols.io step lists 80 °C before 37 °C; the paper gives 37 °C then 80 °C, which is the
   only order that works 🟡.)
7. **cDNA PCR**: + 30 µL mix (Terra direct buffer, SINGV6 [0.2 µM in 50 µL; paper "0.33 µM"
   refers to the 30 µL mix 🟡], Terra 1.25 U). 98 °C 3 min; **13–21 cycles** of 98 °C 15 s,
   65 °C 30 s, 68 °C 4 min (13–15 for ES cells; paper: 15); 72 °C 10 min.
8. **Clean-up** 0.8× beads, elute 15 µL (paper: 10 µL); PicoGreen (> 1 ng/µL expected); Bioanalyzer optional.
9. **Tagmentation**: 0.8 ng cDNA in 20 µL (10 µL TD buffer, 5 µL ATM), 55 °C 10 min;
   5 µL NT to stop, 5 min RT.
10. **3' enrichment PCR**: 15 µL NPM, 0.5 µL 5 µM P5NEXTPT5, 0.5 µL 5 µM N7xx, water to
    50 µL. **72 °C 3 min gap fill**, 95 °C 30 s, **13 ×** (95 °C 10 s, 55 °C 30 s, 72 °C 1 min),
    72 °C 5 min.
11. 1.0× beads, elute 20 µL; **size select 300–900 bp** on 2 % E-Gel EX (paper: 300–800 bp),
    MinElute; Bioanalyzer (> 3–5 ng/µL expected).

### SCRB-seq vs mcSCRB-seq — 🟢 (Bagnoli 2018 Supplementary Table 2, plus protocols)

| | SCRB-seq (Ziegenhain 2017 version) | mcSCRB-seq |
|---|---|---|
| cell suspension | RNAprotect | PBS |
| oligo-dT | added at RT, 1 µM | in the lysis buffer, 0.2 µM |
| RT volume / enzyme | 2 µL (after desiccation), 25 U | 10 µL, 20 U Maxima H− |
| RT enhancer | none | 7.5 % PEG 8000 |
| TSO | 5'-blocked (iso-bases), 1 µM | unblocked, 2 µM |
| pooling | Zymo column | magnetic beads |
| cDNA PCR | KAPA HiFi, 18–21 cycles (original 2014: Advantage 2) | Terra direct, 13–15 cycles |
| time / cost | 2 days, 1–2 € per cell | 1 day, 0.4–0.6 € per cell |

The library architecture (oligo-dT, handles, P5NEXTPT5, Nextera XT i7) is unchanged.

## 5. Final library — 🟡 (assembled from the oligos above; agrees with SCRB-seq Fig. S1 🟢)

Only the 3'-end fragment is amplified (top strand, 5'→3'):

```
5'- P5 · TruSeq-R1 (33) · <bc6> · <umi10> · T30 · VN · <cDNA insert, antisense to the mRNA> · ME' (19) · s7' (15) · <i7'> (8) · P7' (24) -3'
```

Written with `lib/` names: `illumina.P5` · `illumina.TRUSEQ_READ1[4:]`-overlap (the first
4 nt of Read 1 are the last 4 of P5) · BC6 · UMI10 · T30VN · insert · `nextera.ME_RC` ·
`nextera.S7_RC` · i7 · `illumina.P7_RC`. Fixed (non-insert) length: 58 + 16 + 32 + 19 + 15
+ 8 + 24 = **172 nt** 🟡 *(computed)*. With the 300–800 bp size selection, inserts of
roughly 130–630 bp 🟡.

The insert is read on the top strand starting with the poly(T), i.e. **Read 2 (from the
Nextera side) is sense to the mRNA** and Read 1 is antisense — the strand is preserved 🟡.

## 6. Sequencing / read layout

| | SCRB-seq (2014) 🟢 | mcSCRB-seq 🟢 |
|---|---|---|
| instrument | Illumina HiSeq, paired end | HiSeq 1500 high output |
| Read 1 | **17** cycles: BC6 (1–6) + UMI10 (7–16) + 1 base | **16** cycles: BC (1–6) + UMI (7–16) (`zUMIs -c 1-6 -m 7-16`) |
| Index 1 (i7) | 8 cycles | 8 cycles (only when multiplexing plates) |
| Index 2 | — | 0 |
| Read 2 | 34 cycles, cDNA | 50 cycles, cDNA |

Sequencing primers 🟡 (not named in either paper; implied by the library): Read 1 =
`illumina.TRUSEQ_READ1` (anneals to the TruSeq handle from P5NEXTPT5/E3V6NEXT); Index 1 =
`nextera.INDEX1_PRIMER` (`ME_RC + S7_RC`); Read 2 = `nextera.READ2_PRIMER` (`S7 + ME`).
A mixed TruSeq-Read-1 / Nextera-Read-2 library therefore needs a run whose primer mix
contains both, as HiSeq kits of the period did 🟡. The 2014 analysis kept read pairs whose first 6 bases
matched a designed barcode exactly (all Q ≥ 10) and whose next 10 (UMI) were all Q ≥ 30 🟢.

## 7. Check against upstream (scg_lib_structs `SCRB-seq.html`)

Agreements 🟢/🟡:

- All five oligo sequences on the upstream page match the SCRB-seq supplement and the
  protocols.io table base for base (E3V6NEXT handle, E5V6NEXT, SINGV6, P5NEXTPT5 incl. the
  phosphorothioates, Nextera N7xx).
- The final-library drawing (P5 · TruSeq R1 · 6-nt BC · 10-nt UMI · dT · cDNA · ME · s7 ·
  i7 · P7) and the three sequencing primers agree with Fig. S1 and §5–6 above.
- Upstream's account of which tagmentation products amplify ("product 7 only") agrees with
  the computed primer logic in §3 step 7.

Disagreements / imprecisions:

- **Barcode/UMI order in the E3V6NEXT line.** Upstream writes
  `…CGATCT[6-bp cell barcode][10-bp UMI]T30VN`, which matches the source; but the
  first-strand drawings render the oligo 3'→5' and look reversed — they are consistent once
  read as 3'→5'. No real disagreement 🟡.
- **Blocked TSO drawn as plain `CGC`.** Upstream's cDNA top strand begins
  `CGCACACTCTTT…` and its bottom strand ends with the complementary `…GCG-5'`, i.e. it
  shows the iso-dC/iso-dG/iso-dC as copied. The supplement's Fig. S1 draws the second
  strand as starting at `ACACTCTTT…` and notes that modified nucleotides are omitted; an
  ordinary RT does not copy iso-bases, so the first strand should end opposite the 22-nt
  handle 🟡. Harmless for the library (that end is discarded) but the drawn product is not
  the real one.
- Upstream says the library construction is "the same" for mcSCRB-seq and gives only the
  blocked E5V6NEXT; mcSCRB-seq actually uses the **unblocked** TSO
  (`ACACTCTTTCCCTACACGACGCrGrGrG`) 🟢 (protocols.io).
- Upstream omits the bulk N16 oligo-dT, the 17- vs 16-cycle Read 1, and Exonuclease I.

## 8. Open questions

- 🔴 SCRB-seq 2014: RT temperature/time and cDNA PCR cycle number/cycling with Advantage 2
  (the mcSCRB-seq paper uses "the original SCRB-seq PCR cycling conditions" for its
  polymerase screen but does not restate them).
- 🔴 mcSCRB-seq cell-barcode list (`mcSCRBseq_oligodT.txt`) — not obtained; relation to
  SCRB-seq Table S2 sets unknown.
- 🔴 Why the 17th Read 1 cycle in 2014 (probably just to read into the first T) 🟡.
- 🔴 Why phosphorothioates on P5NEXTPT5 — the exonuclease-protection reading in §3 step 7
  is an inference.
- 🟡 Whether the iso-base block was dropped in mcSCRB-seq at a cost: the authors report
  equal yield and no extra primer artefacts for the unblocked TSO (Supp. Fig. 2c–d), 🟢 for
  the claim, but the experiment is at 10 pg UHRR, not single cells.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched only the two bioRxiv full texts and the upstream page. It
**missed the SCRB-seq bioRxiv supplement** (five adjunct files at
`/highwire/filestream/827/field_highwire_adjunct_files/…`, listed on the
`003236v1.supplementary-material` page) — and that supplement holds every oligo; the
preprint main text has none. It also did not follow the mcSCRB-seq preprint to its
*Nat Commun* version or to protocols.io. All of these were fetched by hand (curl) and
converted with `tools/doctext.py`, which handled them fine. `tools/scrape_primers.py` found
all SCRB-seq oligos in `003236-1.pdf.txt`, but parsed E5V6NEXT as a primer named `iCiGi`
with sequence `CACACTCTTT…` (it took the third iso-base `iC` as a plain `C`). In the
protocols.io table, where sequences wrap across table cells, it glued the word "Biotin"
onto the sequence (`Bioti: nACACT…`), cut the unblocked TSO at the wrap
(`ACACTCTTTCCCTACACGA`, 19 nt, without `CGCrGrGrG`) and cut P5NEXTPT5 after the first
phosphorothioate. `--find` ignores wraps but not modification markers, so the unblocked
TSO (`…CGCrGrGrG`) and the protocols.io P5NEXTPT5 (`…CCG*A*T*C*T`) had to be confirmed
with grep on the de-wrapped text. Reaction order and conditions were read by hand.
