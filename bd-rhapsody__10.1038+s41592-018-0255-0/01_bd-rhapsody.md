# BD Rhapsody WTA — microwell bead capture, random-primer second strand

> **Evidence marking.** 🟢 verbatim from a primary source (the paper, or a BD vendor
> document) · 🟡 derived or inferred · 🔴 not published / not available. Relationships
> marked 🟡 *(computed)* were worked out with `lib/`; the model is `tools/rhapsody.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not verification. Its own author describes its Rhapsody bead and
> primer sequences as an educational guess reconstructed from read data — not a vendor or
> author statement, and therefore not evidence at all. An earlier revision of this note
> reproduced that reconstruction (a `CCCCCCTCTCTCT` bead head, 9-nt cell-label blocks,
> named 12/13-nt linkers, an 8-nt UMI, `(T)18`, two 16-nt "pre-amp" primers and full
> TruSeq library primers). None of it has a source, so none of it is asserted any more.
>
> **What is published.** The BD vendor documents mirrored in the scg_lib_structs
> repository under `data/BD/` are ordinary vendor originals — a filename prefix is not
> provenance, and the manifest records the mirrored and hand-fetched copies as
> md5-identical. From them, 🟢: the **Randomer** sequence, the whole RPE/RPE-PCR/index-PCR
> workflow with its conditions, and the capture-oligo **cartoon order** `Univ | CL | UMI |
> poly(T)`. Not from them, 🔴: every bead base, the cell-label and UMI **lengths**, the
> Universal Oligo, the BD AbSeq Primer and both library primers.

**BD Rhapsody** (formerly BD Resolve) is a commercial BD Biosciences platform. The
catalogue keys it to its first use in a Nature Methods paper:

Yoon S-J, Elahi LS, Pașca AM, Marton RM, Gordon A, Revah O, Miura Y, Walczak EM,
Holdgate GM, Fan HC, Huguenard JR, Geschwind DH, Pașca SP. "Reliability of human cortical
organoid generation." *Nat Methods* 16, 75–78 (2019).
doi:[10.1038/s41592-018-0255-0](https://doi.org/10.1038/s41592-018-0255-0) · PMID 30573846 ·
PMC6677388.

(H. Christina Fan, a co-author, is first author of the CytoSeq paper below.) Papers
this one cites for the platform (not fetched, not in the catalogue row): Fan HC,
Fu GK, Fodor SP, "Combinatorial labeling of single cells for gene expression cytometry",
*Science* 347, 1258367 (2015), doi:10.1126/science.1258367 (the "stochastic barcoding" /
CytoSeq origin, their ref. 15); Birey F et al., *Nature* 545, 54–59 (2017) (their ref. 12,
earlier organoid data on the same platform). The catalogue lists no other papers for
this method.

**The paper contains no oligo sequences for the platform** — it is a user of a commercial
kit. The chemistry below comes from the two BD vendor documents mirrored in the
scg_lib_structs repository's `data/BD/` directory. Where those documents stop, this note
stops: it does not substitute a reconstruction for a vendor statement.

### Sources read

All in `_data/sources/bd-rhapsody__10.1038+s41592-018-0255-0/` (never committed).
The BD files are on disk under `upstream_data_*` names because they were first fetched by
hand (curl); `catalogue/extra_sources.tsv` declares them as `scg:data/BD/*`, and the
manifest's `(duplicate) scg_…` rows record the mirrored copies as md5-identical. Their
`what` column, not their filename, is the authority: they are BD's own documents.

| File | What | Used for |
|---|---|---|
| `s41592-018-0255-0_PMC6677388.html.txt` | the paper (PMC full text) | Methods "Single-cell gene expression (BD Rhapsody system)": capture, the (older) library description, sequencing |
| `s41592-018-0255-0_41592_2018_255_MOESM1_ESM.pdf.txt` | Supplementary Information (Springer) | qPCR primers only — not Rhapsody chemistry |
| `s41592-018-0255-0_41592_2018_255_MOESM2_ESM.pdf.txt` | Reporting Summary | software only |
| `s41592-018-0255-0_41592_2018_255_MOESM3_ESM.xlsx.txt` | Supplementary Table 3 (sample sheet) | nothing chemical |
| `upstream_BD_Rhapsody.html.txt` | scg_lib_structs method page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `upstream_data_GMX_BD-Rhapsody-WTA-alpha-Protocol_UG_EN.pdf.txt` | BD "Whole transcriptome analysis alpha protocol", 23-21179-00, 12/2018 | **Randomer sequence**, RPE, RPE PCR, index PCR, cleanups |
| `upstream_data_GMX_BD-Rhapsody-Single-Cell-Analysis-System-Instrument_UG_EN.pdf.txt` | BD Rhapsody instrument user guide, Doc ID 214062 Rev. 2.0, 02/2019 | lysis, bead retrieval, RT mix, Exonuclease I; bead-oligo cartoon |
| `scg_BD_CLS1.txt`, `_CLS2.txt`, `_CLS3.txt` (on disk as `upstream_data_BD_CLS*.txt`) | three 97-member 9-mer lists kept in the scg_lib_structs `data/BD/` mirror | **provenance not established** 🔴: unlike the two PDFs these are bare sequence lists with no vendor document around them, so this note does not treat them as a BD publication and no length or base from them enters the model |

Not obtained (PMC download gate; the Springer copies above cover the same supplements
except possibly the figure DOCX):

| File | URL |
|---|---|
| Supplementary Figures (DOCX) | https://pmc.ncbi.nlm.nih.gov/articles/instance/6677388/bin/NIHMS1511424-supplement-Supplementary_Figures.docx |
| Supplementary Tables 1, 2, 4 (PDF) | https://pmc.ncbi.nlm.nih.gov/articles/instance/6677388/bin/NIHMS1511424-supplement-Supplementary_Tables_1__2__4.pdf |
| Supplementary Table 3 / Reporting Summary (PMC copies) | same `/bin/` directory |
| Fan et al. 2015 (Science; Cloudflare/paywall) | https://doi.org/10.1126/science.1258367 |
| BD's own bead-oligo bases and lengths, Universal Oligo, AbSeq Primer, Library Forward/Reverse primers and sequencing-primer sites | 🔴 unresolved in the BD guides checked so far; additional relevant kit documents remain unavailable |

---

## 1. What it is

Cells are loaded by limiting dilution into a cartridge of **>200,000 microwells**, then
barcoded **magnetic beads** are added to saturation so that each cell-containing well gets
a bead 🟢 (paper, Methods). Lysis in the well lets poly(A) RNA hybridise to the bead's
oligo-dT; beads are pulled out with a magnet and pooled, and **reverse transcription
happens on the pooled beads in one tube** 🟢. The cDNA stays covalently on the bead, so
beads can be archived and subsampled 🟢 (paper: ~67 % of beads used, rest archived; BD
alpha protocol: leftover beads stored at 4 °C up to 3 months).

The bead oligo carries, in order, a universal region, a cell label, a UMI and an oligo-dT
tail 🟢 (instrument guide cartoon: "Univ: universal oligo; CL: cell label; UMI: Unique
Molecular Identifier"). **The checked guides do not resolve its bases or lengths** 🔴: not its bases, not the
cell-label or UMI lengths, not whether the cell label is one block or several, and not the
poly(T) length. The paper says only that each molecule is tagged "with a molecular index
and cell label indicating its cell" 🟢. The model therefore draws `[UNIVERSAL OLIGO]`,
`[CELL LABEL]` and `[UMI]` as symbolic placeholders, so the drawing does not assert an unverified width.

| | Feature | Compare |
|---|---|---|
| 1 | Microwell + bead capture with a **magnet**, beads retrieved and processed in bulk | Microwell-seq, Seq-Well (also microwells; Drop-seq-type beads) |
| 2 | **Exonuclease I** after RT to clear unused bead oligo-dT | Drop-seq/Seq-Well do the same with Exo I |
| 3 | **No template switching**: the second strand is made **on the bead** by a randomer carrying a TruSeq Read 2 handle, extended by Klenow exo⁻ ("RPE") | 10x / Drop-seq use a TSO; here the cDNA 5' end is never captured |
| 4 | The library is 3'-anchored by the bead and cut at a **random** position by the randomer, so no fragmentation/tagmentation step | Tn5-based 3' kits |

Concept background: [reverse transcription](../ref/concepts/reverse-transcription.md)
(oligo-dT on a bead, random priming). Template switching
([concept](../ref/concepts/template-switching.md)) is **not** used — that is the main
difference from bead protocols that make full-length cDNA first.

### Two different library-prep descriptions — disagreement

- **The paper** (Dec 2018) says the whole-transcriptome libraries were made by
  "second strand cDNA synthesis, adaptor ligation, and universal amplification" with
  **22 PCR cycles**, followed by **random priming PCR** to enrich the 3' ends 🟢
  (paraphrased). That is a ligation-based WTA; no adaptor sequences are given 🔴.
- **The BD WTA alpha protocol** (12/2018) describes a different route: random-primer
  extension on the beads (RPE), RPE PCR (13 cycles), index PCR (8–9 cycles) 🟢. No
  ligation, no 22-cycle step.

So BD's documented commercial WTA workflow is **not** necessarily what Yoon et al. ran.
This note and the page follow the documented RPE route, and say so; the adaptor and
primer sequences of the ligation route remain unverified 🔴.

## 2. Oligos

**From the BD WTA alpha protocol** — the one Rhapsody oligo any source we hold prints, a
third-party reagent the user orders separately 🟢 (no modification is given for it):

```
Randomer   5' TCA GAC GTG TGC TCT TCC GAT CTNNNNNNNNN 3'
```

(written with codon spaces in the source; 100 µM in 10 mM Tris, 0.1 mM EDTA pH 8.0; N
hand-mixed 25 % each 🟢)

**Other ordered bases remain unverified in the checked guides 🔴.** BD's documents name the PCR reagents without
giving a sequence: **"Universal Oligo"** and **"BD AbSeq Primer"** (RPE PCR), and
**"Library Forward Primer"** and **"Library Reverse Primer (1–4)"** (index PCR) 🟢. The
capture bead is supplied pre-made; its exact sequence remains unverified 🔴. No sequencing-primer
site is given 🔴, and BD names four indexed library reverse primers without printing an
index 🔴.

So the only Rhapsody oligo this note can state is the Randomer. The model reflects exactly
that: `randomer` is the one graph input with a real sequence; `capture_beads`,
`universal_oligo`, `rpe_reverse`, `library_forward` and `library_reverse` are
`inferred=True` inputs, the RPE PCR and library PCR both carry
`PCRGeometry.unavailable(...)`, the final library is drawn as
`[P5 / READ 1 ARM] … [i7 / P7 ARM]`, and the sequencing section renders
`seqprimers.unavailable_diagram(...)` rather than a primer-binding diagram.

## 3. How the one published oligo fits — 🟡 *(computed with `lib/illumina`)*

- **Randomer** = the 3'-terminal **23 nt of `illumina.TRUSEQ_READ2`** (`TCAGACGTGTGCTCTTCCGATCT`)
  + N₉ (32 nt). So every random-primed second strand begins with a partial TruSeq Read 2
  handle, and the library's Read 2 side is a TruSeq-style end. 🟡 *(computed)*
- That is the **only** sequence relationship this protocol's sources support. The bead's
  own handle, the Universal Oligo, the library primers and the index remain unresolved in the checked guides 🔴,
  so no P5/Read 1 arm, linker arrangement, read length requirement or final-library
  length can be derived here. An earlier revision of this note derived all of those from
  the drawn page's reconstruction; those derivations are withdrawn.

## 4. Step by step

Steps 1–4 from the BD instrument user guide (Doc 214062 Rev. 2.0), steps 5–11 from the
BD WTA alpha protocol; conditions 🟢 unless marked. Strand bookkeeping 🟡.

1. **Load cells** (~10,000 in the paper 🟢; 5,000–10,000 per sample recommended for WTA
   🟢) into the cartridge microwells by limiting dilution, then load **Cell Capture
   Beads** to saturation, image with the scanner.
2. **Lyse**: 550 µL Lysis Buffer + DTT (75 µL 1 M DTT per 15 mL), room temperature
   **2 min**. Poly(A) RNA anneals to the bead `(T)18` 🟡 (the paper says only that poly(A)
   RNA hybridises to the beads).
3. **Retrieve beads** magnetically into a 5 mL tube with 4,950 µL lysis buffer; wash.
4. **RT on beads**, 200 µL cDNA mix (RT buffer, dNTP, 0.1 M DTT, Bead RT/PCR Enhancer,
   RNase inhibitor, reverse transcriptase), **37 °C 20 min, 1,200 rpm** ("shaking is
   critical"). First strand: `|--bead oligo · (T)18 · cDNA →`, the cDNA covalently on the
   bead with cell label and UMI at its 5' end 🟢 (paper: tagged on the cDNA 5' end =
   mRNA 3' end).
5. **Exonuclease I**, 200 µL, **37 °C 30 min**, then **80 °C 20 min** to inactivate.
   Removes bead oligos that captured nothing (single-stranded 3' ends) 🟡 (purpose not
   stated in the BD text).
6. **Remove the RNA**: beads 95 °C 2 min, magnet, discard supernatant. The bead now holds
   single-stranded cDNA 🟡.
7. **Random priming**: 174 µL Random Primer Mix (water, NEBuffer 2, 20 µL Randomer);
   95 °C 2 min (no shaking) → 37 °C 5 min → 25 °C 5 min (both 1,200 rpm). N₉ anneals at random positions
   along the bead-bound cDNA.
8. **Extension (RPE)**: + 26 µL (dNTP, Bead RT/PCR Enhancer, **Klenow exo⁻**), 37 °C
   **30 min**, 1,200 rpm. The randomer extends toward the bead and copies the bead oligo,
   so the new strand is `5'- R2 handle(23) · N₉-insert … · (A)18 · UMI' · CLS3' … CLS1' ·
   R1-handle' · AGAGAGAGGGGGG -3'` (complement written 5'→3') 🟡. Only extensions that
   reach the bead oligo carry the Read 1 handle.
9. **Release**: wash, resuspend in 200 µL Tris-Tween20, **95 °C 2 min**, magnet, keep the
   **supernatant = RPE product** (beads can be stored). Clean up with **1.6× SPRIselect**
   (320 µL to 200 µL), elute 30 µL.
10. **RPE PCR**: 50 µL PCR MasterMix + 10 µL Universal Oligo + 10 µL BD AbSeq Primer +
    30 µL RPE product, split in two; 95 °C 3 min; **13 ×** (95 °C 30 s, 60 °C 1 min,
    72 °C 1 min); 72 °C 2 min. 1× SPRI. Expected broad 200–2,000 bp; quantify the
    150–600 bp region. What the two sequence-unresolved primers add to each end, and therefore
    where the amplicon boundaries fall, is 🔴; the model expresses this as
    `PCRGeometry.unavailable(...)` rather than drawing a boundary.
11. **WTA Index PCR**: dilute to 2 nM (150–600 bp peak); 10 µL + Library Forward +
    one Library Reverse (1–4) in PCR MasterMix (50 µL); 95 °C 3 min; **9 cycles (1–2 nM)
    or 8 cycles (>2 nM)** of 95 °C 30 s / 60 °C 30 s / 72 °C 30 s; 72 °C 1 min.
12. **Double-sided SPRI**: + 60 µL water (110 µL), take 100 µL, + 60 µL beads (0.6×, keep
    supernatant), + 10 µL beads to the supernatant (0.7× cumulative), keep beads, elute
    30 µL 🟢 volumes / 🟡 ratios. Final libraries **~300–700 bp**, >1 ng/µL 🟢.

## 5. Final library — 🔴 not derivable

The library's **Read 2 side** is a TruSeq-style end built by the Randomer 🟡 *(computed,
§3)*, and the insert runs from wherever the randomer landed to the poly(A) junction 🟡
*(strand bookkeeping; because the top strand is cDNA, Read 2 reads the mRNA sense
strand)*.

The remaining base-level details are unresolved in the checked guides 🔴: the P5/Read 1 arm, the bead-derived barcode block
(bases, lengths and the number of cell-label segments), the index, and therefore the total
adapter length. The model draws this honestly —

```
5'- [P5 / READ 1 ARM] · [CELL LABEL] · [UMI] · poly(T)/cDNA · TCAGACGTGTGCTCTTCCGATCT · [i7 / P7 ARM] -3'
```

— with every bracketed region `inferred=True`. An earlier revision of this note printed a
complete 202-nt adapter assembly for "V1" and "Enhanced" beads, derived from the drawn
page's reconstruction, together with a computed insert-size range. Those numbers had no
source and are withdrawn.

## 6. Sequencing

- 🟢 Paper: HiSeq 2500, **101 × 2**; ~38,000 reads/cell; molecular index coverage 3.2–3.5.
- 🟢 BD alpha protocol: NextSeq / MiniSeq High or Mid output, load 1–1.2 pM with **20 %
  PhiX**; 10,000 / 50,000 / 100,000 reads per cell for shallow / moderate / deep.
- 🟢 Paper: "Cell labels and molecular indices were identified" from the reads, i.e.
  Read 1 carries the cell label and UMI and Read 2 the cDNA. Which primer site each read
  starts from is 🔴 — the checked guides provide no sequencing-primer sequence, so the page renders
  `seqprimers.unavailable_diagram(...)` and declares the run layout rather than computing
  it.
- 🔴 How many cycles the cell label and UMI occupy, and therefore where Read 1 runs into
  poly(T). An earlier revision of this note stated "≥60 cycles" and a cycle-by-cycle
  breakdown taken from the drawn page; that is withdrawn.

## 7. Open questions

- 🔴 **What the paper actually ran**: its Methods describe ligation-based WTA
  (second strand, adaptor ligation, 22-cycle universal PCR, then random-priming PCR),
  which does not match BD's 12/2018 RPE alpha protocol. The adaptor and primer sequences
  of that route remain unverified.
- 🔴 The bead-oligo sequence itself — universal region, cell-label structure and length,
  UMI length, poly(T) length, any anchor — remains unresolved in the checked guides. Other relevant BD kit documents have not been checked,
  and read-derived reconstructions are not evidence (see the preamble). The purpose of the
  "BD AbSeq Primer" in a WTA-only RPE PCR is also unexplained.
- 🔴 "Universal Oligo" sequence: BD gives only the name (one tube, so either a single
  oligo or a premixed pair).
- 🟡 Strand displacement: Klenow exo⁻ has limited strand-displacement activity, so when
  several randomers prime the same cDNA, only bead-proximal products may reach the
  Read 1 handle — this would bias inserts toward the 3' end, consistent with the
  "3' WTA" description, but nobody states it.
- 🔴 Provenance of the three 97-member CLS lists kept in the scg_lib_structs `data/BD/`
  mirror. They are bare sequence files with no vendor document around them, so this note
  does not treat them as a BD publication; if BD or a BD pipeline document publishes the
  cell-label whitelist, that document is what should be fetched and cited.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC text, JATS XML and the three Springer
supplements; the four PMC `/bin/` supplements stopped at the download gate (manual), and
the Europe PMC supplementary zip held only a placeholder GIF entry. None of these
contain Rhapsody chemistry. The BD vendor documents (alpha WTA protocol,
instrument user guide, CLS lists) were first fetched by hand with curl into the sources
directory and converted with `tools/doctext.py`; `catalogue/extra_sources.tsv` now
declares them as `scg:data/BD/*`, and because the hand-fetched copies are already on disk
`get_sources.py` records the mirrored ones as md5-identical duplicates rather than
re-saving them. `tools/scrape_primers.py` in listing mode **missed the Randomer in the BD
alpha protocol** because it is written in codon-spaced triplets (`TCA GAC GTG …`);
`--find` does locate it there. The guard
`tools/selftest_source_oligos.py` cannot yet check the Randomer against that PDF: the
file sits under an `upstream_data_*` name, which the guard excludes by prefix, and the
Randomer's N₉ tail means it has no exact whole-strand state either.
