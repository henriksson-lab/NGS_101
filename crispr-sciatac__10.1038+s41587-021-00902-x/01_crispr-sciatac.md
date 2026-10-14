# CRISPR-sciATAC — research notes

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/`; the oligo
> transcriptions are pinned in `tools/selftest.py`, and `tools/selftest_source_oligos.py`
> re-finds every modelled oligo in the extracted source text.

**CRISPR-sciATAC** — Liscovitch-Brauer N, Montalbano A, Deng J, Méndez-Mancilla A,
Wessels H-H, Moss NG, Kung C-Y, Sookdeo A, Guo X, Geller E, Jaini S, Smibert P,
Sanjana NE. "Profiling the genetic determinants of chromatin accessibility with scalable
single-cell CRISPR screens." *Nature Biotechnology* 39:1270–1277 (2021).
doi:[10.1038/s41587-021-00902-x](https://doi.org/10.1038/s41587-021-00902-x). Preprint
(as "Scalable pooled CRISPR screens with single-cell chromatin accessibility profiling"):
doi:[10.1101/2020.11.20.390971](https://doi.org/10.1101/2020.11.20.390971).

Sources read (fetched into `$CHEM_DATA/sources/crispr-sciatac__10.1038+s41587-021-00902-x/`,
never committed):

| File | What | Used for |
|---|---|---|
| `s41587-021-00902-x_PMC8516442.xml(.txt)` | author manuscript, full Methods | reaction order, TnY, MEDS A/B, RT, the two PCRs, purification |
| `s41587-021-00902-x_41587_2021_902_MOESM4_ESM.xlsx` | **Supplementary Table 2** | **every oligo**: common oligos + six complete barcode plates |
| `scg_41587_2021_902_MOESM1_ESM.pdf` · `..._MOESM1_ESM.pdf` | Supplementary Information, **Supplementary Fig. 1a–e** | barcode naming, construct layout, primer design, run lengths |
| `2020.11.20.390971_v1.full.pdf` | bioRxiv v1 full text | cross-check of the Methods |
| `nmeth.4177_*` (CROP-seq) | inherited guide cassette | CROPseq-Guide-Puro map (Addgene 86708), scaffold geometry |
| `414870_*`, `j.cell.2018.11.022_*` (Perturb-ATAC) | inherited comparison method | cost/time comparison only; no sequence enters this model |

The Teichlab scg_lib_structs page for this protocol was **not** used as evidence; the
`scg_41587_2021_902_MOESM1_ESM.pdf` in that repository is the publisher's own
Supplementary Information file, mirrored, and is cited as such.

---

## 1. What it is

Pooled CRISPR screening with a single-cell ATAC read-out, built on sci-ATAC-style
combinatorial indexing. 🟢 Cells carrying a CROP-seq guide cassette are fixed with
glyoxal, lysed to nuclei, and distributed into a 96-well plate where a **barcoded
transposome** tags accessible chromatin. 🟢 In the same well the polyadenylated CROP-seq
transcript is reverse-transcribed with a primer carrying **the same barcode**. 🟢 Nuclei
are then pooled and redistributed at ~15–20 per well, and two more barcodes arrive on the
PCR primers.

🟢 "Cell barcodes consist of a unique combination of BC1, BC2, and BC3" (Supplementary
Fig. 1a legend), and "BC 1, 2, and 3 are matched for ATAC-seq and gRNA libraries, e.g.
the ATAC-seq BC 1 in well A1 … is the same as the gRNA BC 1 in well A1 …".

| | Where the barcode is installed | ATAC branch | Guide branch |
|---|---|---|---|
| BC1 | tagmentation well | barcoded T7 transposome arm (MEDS B) | barcoded RT primer |
| BC2 | PCR well | P5i ATAC PCR1 primer | P5i sgRNA PCR1 primer |
| BC3 | PCR well | P7i ATAC PCR1 primer | P7i sgRNA PCR2 primer |

🟢 All three are **8 nt** — every row of all six plates in Supplementary Table 2
decomposes into fixed flanks around an 8-nt window (🟡 *computed*: 96 BC1, 96 BC2 and 144
BC3, all unique, and the BC1/BC2/BC3 series are in the **same order** on the ATAC and the
gRNA plate of each pair).

🔴 There is **no sample index anywhere in this protocol**. Every barcode position,
including the ones that sit where i5 and i7 normally sit, is part of the cell identity.
Libraries from different experiments are kept apart by sequencing them separately
(🟢 "Libraries were separately sequenced on the MiSeq Sequencer").

### The transposase is TnY, not Tn5

🟢 "We found a range of transposon sequences that were related to the Tn5 sequence and
selected a transposon from *Vibrio parahemolyticus* (ViPar)… Two mutations were
introduced: (1) P50K… and (2) M53Q… The ViPar transposase with P50K and M53Q mutations,
henceforth referred to as **TnY**". 🟢 "The inside and outside ends (IE and OE) of the
ViPar transposon utilize the same sequence as the IE and OE of the Tn5 transposon", so the
adapters carry the **canonical 19-nt Tn5 mosaic end** and nothing else about the ends
changes. 🟢 The polymerase is likewise a home-made PfuX7, not a kit enzyme.

## 2. Oligos

🟢 Verbatim from Supplementary Table 2. `[8-nt BCn]` marks the plate window; the plates
themselves are not reproduced here.

```
pMENT                        /5Phos/CTGTCTCTTATACACATCT
T5 tagment sciATAC           TCGTCGGCAGCGTC TCCACGCAACCAACG CACCGTCTCCGCCTC AGATGTGTATAAGAGACAG
Barcoded T7 tagment sciATAC  GTCTCGTGGGCTCGG CTGTCCCTGTCC [8-nt BC1] GCGATCGAGGACGGC AGATGTGTATAAGAGACAG
RTi sgRNA                    CCGGGAGCTGCATGT [8-nt BC1] TTTTTCAAGTTGATAA
P5i ATAC (PCR1)              AATGATACGGCGACCACCGAGATCTACAC [8-nt BC2] TCGTCGGCAGCGTC
P7i ATAC (PCR1)              CAAGCAGAAGACGGCATACGAGAT [8-nt BC3] GTCTCGTGGGCTCGG
P5i sgRNA (PCR1)             CGCCGTCCAGCTGA [8-nt BC2] CCGGGAGCTGCATGT
U6 outer sgRNA PCR1          TGCATATACGATCACAAGGCTGTTAG
P5 sgRNA PCR2                AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT CGCCGTCCAGCTGA
P5 sgRNA PCR2 +1/+2/+3/+8    ...CTTCCGATCT [A | TA | GTA | NNNNNNNN] CGCCGTCCAGCTGA
P7i sgRNA (PCR2)             CAAGCAGAAGACGGCATACGAGAT [8-nt BC3] TTGTGGAAAGGACGAAACAC
Read 1 ATACseq               CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
Index 1 ATACseq              CTGTCTCTTATACACATCTGCCGTCCTCGATCGC
Read 2 ATACseq               GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
```

(Supplementary Table 2 also lists a qPCR primer/probe panel, used for the arrayed
EZH2/HOX validation by TaqMan qRT-PCR, which builds no sequencing library.)

### The four spike-in oligos

🟢 Supplementary Table 2 also prints `spike-in T5`, `spike-in T7`, `spike-in P5` and
`spike-in P7`. 🔴 **The paper never mentions them**: the word "spike" appears nowhere in
the Methods, the preprint, the Supplementary Information or any other supplement — only in
the oligo table. So this repository cannot say from its defining source what they are for,
at what ratio they were added, or at which step.

🟡 *(computed)* They are a complete, self-consistent, **cell-barcode-free** construct with
fixed indexes, not a reagent of the library above: `spike-in T7` carries a fixed 8-nt
`CGAGTAAT` where BC1 would sit, `spike-in P5`/`P7` carry fixed i5 `TCCTACGG` and i7
`CTAAGAAC`, and both transposon arms reproduce this protocol's Read 1 and Read 2 entry
points exactly — so a spike molecule is read by the same three custom sequencing primers
on the same run, while its P5/P7 handles cannot cross-prime the sciATAC library.

🟡 The set is inherited verbatim from snATAC-seq (Preissl et al. 2018): the same
`CTAAGAAC`/`TCCTACGG` indexes and the same T7 prefix appear there as `p7_spike-in`,
`p5_spike-in`, `X755_spike-in` and `X536_spike-in`, and that family's purpose is stated by
Amini et al. 2014 — a spike-in library of different transposon and PCR sequences used "to
balance the bases within each detection cycle", because Index 1 and Index 2 each hold two
barcodes around a common linker. See
[snATAC-seq](../snatac-seq__10.1038+s41593-018-0079-3/01_snatac-seq.html), which analyses
the same design. 🔴 That purpose is inference from the inherited reagent, not a claim this
paper makes.

🟡 *(computed)* Two divergences from the snATAC-seq set are worth recording, because they
mean the oligos were carried over rather than re-derived for this protocol: the T5 arm uses
`revcomp(s5)` where snATAC-seq uses `complement(s5)`, and its final 8 bases `TATAGCCT` are
snATAC-seq's barcode-slot residue, not the complement of this paper's own T5 spacer
`AACCAACG`.

How they interlock, 🟡 *(computed)* and enforced in `tools/chem.py`:

- `pMENT` **is** `nextera.ME_RC` exactly, so MEDS A and MEDS B are ordinary loaded Tn5
  adapters with only the terminal 19-nt ME duplexed.
- `T5 tagment sciATAC` = `nextera.S5` + a 15-nt **T5 spacer** + the 15-nt Read 1 entry
  point + `nextera.ME`. 🟢 Supplementary Fig. 1d labels these "Spacer on T5" and
  "MEDS A"; `Read 1 ATACseq` is exactly that last 34 nt.
- `Barcoded T7 tagment sciATAC` = `nextera.S7` + a 12-nt **T7 spacer** + BC1 + the 15-nt
  Read 2 entry point + `nextera.ME`; `Read 2 ATACseq` is again that last 34 nt.
- `Index 1 ATACseq` = `ME_RC` + `revcomp(GCGATCGAGGACGGC)`. It reads *out of* the mosaic
  end **towards P7**, i.e. through BC1, the T7 spacer, s7 and BC3 — **43 cycles**.
- The RT primer's 3′ arm `TTTTTCAAGTTGATAA` is `revcomp(TTATCAACTTGAAAAA)`, a site inside
  the **sgRNA scaffold**; so guide identity is captured from RNA, not from genomic DNA.
- `P5 sgRNA PCR2` = `illumina.P5` + `illumina.TRUSEQ_READ1` (sharing P5's 3′ `ACAC`) +
  the `CGCCGTCCAGCTGA` handle that the PCR1 primer installed.
- `P7i sgRNA PCR2`'s 3′ arm `TTGTGGAAAGGACGAAACAC` is the U6 promoter 3′ end, **nested
  inside** the U6 outer PCR1 primer: that is what trims the amplicon.

## 3. Step by step

1. 🟢 **Transduce** at MOI ~0.1 with a pool cloned into CROPseq-Guide-Puro (Addgene
   86708), "which leads to the synthesis of an RNA Pol3 transcript of the Cas9 gRNA and an
   RNA Pol2 polyadenylated transcript containing the puromycin resistance gene, a U6
   promoter, and the gRNA"; select with puromycin for one week.
2. 🟢 **Fix** with glyoxal (chosen over PFA for nucleus integrity) and lyse to nuclei.
3. 🟢 **Anneal** pMENT to the single T5 oligo ("MEDS A") and to each of the 96 barcoded T7
   oligos ("MEDS B"), mix, and assemble onto TnY.
4. 🟢 **Tagment** ~20,000 nuclei per well, 30 min at 37 °C; stop with EDTA, then quench
   the EDTA with MgCl₂ before RT. 🟡 Only the T5/T7 product can bridge P5 to P7; T5/T5 and
   T7/T7 products carry the same arm at both ends (`nextera.TAGMENTATION_OUTCOMES`).
5. 🟢 **Reverse-transcribe** 5 µl of nuclei into a plate of barcoded RT primers, "keeping
   plate orientation to match tagmentation and reverse transcription barcodes";
   RevertAid, 30 min at 37 °C.
6. 🟢 **Pool** all wells, count, and redistribute ~20 nuclei (2 µl at 10 nuclei/µl) per
   well of a new 96-well plate; **digest** with SDS/proteinase K at 65 °C and quench with
   PMSF.
7. 🟢 **PCR1** is a single multiplexed reaction: ATAC-seq primers at 0.5 µM *and*
   gRNA-PCR1 primers at 0.1 µM, PfuX7 in Phusion GC buffer, 72 °C 5 min gap fill then
   14–18 cycles. The ATAC branch is finished here (BC2 + BC3); the guide branch only gains
   BC2 and the `CGCCGTCCAGCTGA` handle.
8. 🟢 **PCR2** (guide only): 2 µl transferred "keeping plate orientation to match ATAC-seq
   and gRNA barcodes"; 20 cycles with the P5 sgRNA PCR2 and P7i sgRNA PCR2 primers, which
   add P5 + TruSeq Read 1 and BC3.
9. 🟢 **Purify**: ATAC by column + 0.9X Ampure XP; guide by column, 2 % E-gel, **~250 bp**
   band excised. 🟡 *(computed)* the modelled unstaggered guide library is 233 bp, so the
   staggered members of the pool fall at 233–241 bp.

## 4. Final libraries and read layout

🟡 *(computed, from the oligos above)*

```
ATAC   P5 - BC2 - s5 - T5spacer - Read1site - ME - [accessible gDNA] - ME_rc
          - Read2site_rc - BC1_rc - T7spacer_rc - s7_rc - BC3_rc - P7_rc
guide  P5 - TruSeqRead1 - CGCCGTCCAGCTGA - BC2 - CCGGGAGCTGCATGT - BC1
          - TTTTTCAAGTTGATAA - [scaffold_rc] - [spacer_rc] - U6end_rc - BC3_rc - P7_rc
```

🟢 Supplementary Fig. 1d,e give the run:

| Library | Read 1 | Index 1 | Index 2 | Read 2 |
|---|---|---|---|---|
| ATAC | 45 cycles, genomic | 43 cycles: BC1 + T7 spacer + s7 + BC3 | 8 cycles: BC2 | 45 cycles, genomic |
| guide | **162 cycles**, covering BC2, BC1, the scaffold, the **guide spacer**, the U6 end and BC3 | — | — | — |

Two consequences worth stating plainly:

- 🟡 The ATAC **Index 1 read carries two different cell barcodes** separated by 27 nt of
  constant sequence. It is not a sample-index read and is never used for demultiplexing.
- 🟢 The gRNA library takes **no index reads at all**; one long Read 1 does everything,
  which is why the P5 PCR2 primer pool carries 0/+1/+2/+3/+8 diversity bases ("Staggered
  P5 oligos were introduced in the library preparation to introduce sequence diversity").
- 🟡 Supplementary Table 2 lists **no i5 sequencing primer**. The ATAC runs were on a
  MiSeq and a HiSeq 2500, both forward-strand-workflow instruments, where the flow-cell P5
  oligo primes the i5 read; that is how the 8-cycle BC2 read is modelled here.
- 🟡 Guide identity is the **20-nt spacer itself**, read in Read 1 — not a separate guide
  barcode. Every gRNA in Supplementary Table 2 is 20 nt (338 of them).

## 5. What is *not* transcribed, and why

- 🟡 The 42-nt stretch of **sgRNA scaffold** between the guide spacer and the RT-primer
  site is drawn as a length-preserving placeholder. Its length is read off the
  CROPseq-Guide-Puro map (Addgene 86708) that the Methods name — the hU6 promoter ends at
  the `TTGTGGAAAGGACGAAACAC` site, the scaffold begins `GTTTTAGAGCTAGAAATAGCAAG…` and the
  RT site `TTATCAACTTGAAAAA` is 42 nt further on. 🔴 The defining paper itself never prints
  the scaffold.
- 🟡 The 1–2-nt BsmBI cloning junction between the U6 3′ end and the spacer is not drawn;
  it shifts every later Read 1 cycle by that much.
- 🔴 The **U6 outer PCR1 primer has no exact site on a CROPseq-Guide-Puro template**:
  Supplementary Table 2 prints `TGCATATACGATCACAAGGCTGTTAG`, while the vector carries
  `TGCATATACGATACAAGGCTGTTAG` (one base different, one shorter) at the corresponding
  position in the hU6 promoter. This is why `guide_pcr1` uses *regional* rather than exact
  PCR geometry. Open question: a table typo, or a different U6 variant in the lab's stock?
- 🔴 The accessible genomic insert, and the identity of the 96/144 barcodes, are plate
  data rather than construct sequence.

## 6. Inherited components

- **CROP-seq** (Datlinger et al. 2017, doi:10.1038/nmeth.4177; preprint
  doi:10.1101/083774) supplies the guide cassette: the Pol II transcript that carries the
  U6 promoter and the gRNA into the poly(A) pool and makes RT capture possible.
- **sci-ATAC-seq** (Cusanovich et al. 2015, doi:10.1126/science.aab1601) supplies the
  two-round combinatorial indexing logic.
- **Perturb-ATAC** (Rubin et al. 2019, doi:10.1016/j.cell.2018.11.022; preprint
  doi:10.1101/414870) is the method CRISPR-sciATAC is compared against on cost and time
  (Supplementary Tables 8 and 9). 🔴 None of its chemistry enters this model.

🔴 The phrase "semi-suppressive PCR" does **not** occur in this paper, its preprint or any
of its supplements; it was carried into an earlier version of this page from a secondary
schematic. What it was reaching for is the ordinary Tn5 end-selection above: only the
T5/T7 heteroduplex has one P5-capable and one P7-capable end.
