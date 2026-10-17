# SMART-seq3 — quick reference

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not published.
>
> SMART-seq3xpress shares this architecture; its own oligo source review is
> `01_smart-seq3xpress.md`, in the `smart-seq3xpress__10.1038+s41587-022-01311-4`
> directory.

**SMART-seq3** — Hagemann-Jensen M, Ziegenhain C, Chen P, Ramsköld D, Hendriks G-J,
Larsson AJM, Faridani OR, Sandberg R. "Single-cell RNA counting at allele and isoform
resolution using Smart-seq3." *Nat Biotechnol* 2020;**38**(6):708–714.
doi:[10.1038/s41587-020-0497-0](https://doi.org/10.1038/s41587-020-0497-0)

**SMART-seq3xpress** — Hagemann-Jensen M, Ziegenhain C, Sandberg R. *Nat Biotechnol*
2022;**40**(10):1452–1457. doi:[10.1038/s41587-022-01311-4](https://doi.org/10.1038/s41587-022-01311-4)

Sources read (all primary):

- 🟢 Hagemann-Jensen et al. 2020 **Supplementary Table 1** ("Overview of Conditions"),
  Springer file `41587_2020_497_MOESM3_ESM.xlsx` — prints, per condition, the oligo-dT,
  the TSO with its UMI structure, and the forward and reverse PCR oligos. The Smart-seq3
  rows give `/5BiosG/ACGAGCATCAGCAGCATACGA(T)30VN`,
  `/5Biosg/AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG`,
  `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAATG` and `ACGAGCATCAGCAGCATACGA`.
- 🟢 The authors' **protocols.io** protocol `smart-seq3-protocol-36wgq5rjxgk5` — the ordered
  oligo table, which adds the supplier, the 5' biotin and the two terminal
  phosphorothioate linkages on each PCR primer.
- 🟢 The **preprint** Methods (bioRxiv 817924v1) restate the same four oligos inline, and
  give the `find_pattern: ATTGCGCAATG` used to call UMI-containing reads.
- 🟢 Hagemann-Jensen et al. 2022 Methods (PMC9546772) for the `WW`-spacer xpress TSO; see
  that protocol's own note.
- 🟡 Illumina's adapter-sequence documentation, via `lib/nextera.py` and `lib/illumina.py`,
  for ME, s5, s7, P5, P7 and the index primers.

Every modelled base on the SMART-seq3 page is a transcription of one of the above. No
secondary redrawing of this family is used as evidence anywhere in this protocol.

---

## 1. What changed, and why it matters

🟢 SMART-seq3 is a redesign of the oligos, not an optimisation. Four changes:

1. **Maxima H− reverse transcriptase** (an MMLV-derived RNase H− enzyme) instead of
   Superscript II; Supplementary Note I reports that swapping Superscript II for Maxima H−
   under Smart-seq2 conditions did not help, so the buffer was re-optimised with it.
2. **NaCl instead of KCl** during reverse transcription.
3. **5% PEG** as a crowding reagent, as in mcSCRB-seq.
4. **An 11-bp tag and an 8-bp UMI in the TSO**, so that a 5' read can be distinguished from
   internal reads.

Point 4 is the architectural one. SMART-seq/2 made both ends of the cDNA identical; SMART-seq3
deliberately makes them **different**, which buys three things at once:

- **PCR is no longer single-primer** — a forward and a reverse primer, each matching one end.
- **A read carrying the 11-bp tag is known to come from the transcript's 5' end**, because
  template switching only happens where the polymerase ran off the template.
- **The UMI counts molecules**, since it was attached before any amplification.

So one protocol gives **full-length coverage *and* UMI counting at allele and isoform
resolution** — which is the headline result.

🟢 **SMART-seq3xpress** retains this asymmetric library architecture while changing the
reaction volumes and using a `WW` TSO spacer. FLASH-seq is a separate SMART-handle
chemistry, documented in its own note.

## 2. Oligos

```
Smartseq3_OligodT30VN  5'- /5Biosg/ ACGAGCATCAGCAGCATACGA TTTTT...TTTTT VN -3'   (30 T)

Smartseq3_N8_TSO       5'- /5Biosg/ AGAGACAGATTGCGCAATG [8-bp UMI]       rGrGrG -3'
Smartseq3xpress_TSO    5'- /5Biosg/ AGAGACAGATTGCGCAATG [8-bp UMI] WW    rGrGrG -3'

Fwd_PCR_primer         5'- TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G -3'
Rev_PCR_primer         5'- ACGAGCATCAGCAGCATAC*G*A -3'
```

🟢 `*` is a phosphorothioate linkage (3'-exonuclease protection) and `/5Biosg/` is a
5' biotin; both are printed as such in the author protocol. The index primers are
🔴 **not** asserted here: the author protocol orders *custom* Nextera index primers and
recommends a 10-bp barcode, noting 8 bp (Illumina's length, which the page draws as an
index window) works as well, and no index member is chosen.

🟡 Note the TSO handle: `AGAGACAG` is the **3' end of the Nextera mosaic end**, and
`ATTGCGCAATG` is the 11-bp tag. So the forward PCR primer
(`s5 + ME + tag` = `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` + `ATTGCGCAATG`) rebuilds a complete
Nextera s5 arm as it amplifies — the TSO is already most of one.

SMART-seq3 and SMART-seq3xpress differ in the spacer between UMI and G-tail: none or
`WW`. FLASH-seq is a separate SMART-handle chemistry: its main branch has no UMI, while
the optional FS-UMI branch has a different `CTAAC` TSO. Its oligos are recorded in the
FLASH-seq note and must not be substituted for SMART-seq3 oligos.

## 3. Steps

🟢 From the author protocol V.3 (protocols.io DOI `10.17504/protocols.io.bcq4ivyw`,
PDF pp5–19; the version on disk is dated 19–20 Feb 2020). The published Methods are
not readable here (subscription preview), so the article's own wording of these
conditions was not compared.

| # | Step | Source conditions |
|---|---|---|
| 1 | Sort into 3 µL lysis; lyse 72 °C 10 min | 0.1% Triton X-100, 5% PEG 8000, 0.5 mM dNTP each, 0.5 µM oligo-dT30VN, RNase inhibitor (concentrations given for the 4 µL RT) |
| 2 | Add 1 µL RT mix; reverse transcribe with **Maxima H−** | 25 mM Tris-HCl pH 8.3, 30 mM NaCl, 2.5 mM MgCl2, 1 mM GTP, 8 mM DTT, 2 µM TSO, 2 U/µL RT; 42 °C 90 min, 10×(50 °C 2 min, 42 °C 2 min), 85 °C 5 min |
| 3 | Untemplated C addition and template switch onto the TSO in the same reaction | tag + UMI now attached, before amplification (🟡 CCC drawn as the conventional model) |
| 4 | Pre-amplify in 6 µL **KAPA HiFi HotStart** mix with **Fwd_PCR_primer (0.5 µM) and Rev_PCR_primer (0.1 µM)** | 98 °C 3 min; input-dependent cycles of 98 °C 20 s / 65 °C 30 s / 72 °C 4 min; 72 °C 5 min |
| 5 | Purify cDNA with 0.6:1 22% PEG / Sera-Mag beads ("preferable but optional") | elute 12 µL |
| 6 | Quantify (QuantiFluor) and normalise to 100 pg/µL | optional but recommended |
| 7 | **Tagment** 1 µL cDNA with 1 µL mix (0.08 µL Nextera XT Amplicon Tagment Mix, 1× home-made TD buffer: 10 mM Tris pH 7.5, 5 mM MgCl2, 5% DMF) | 55 °C 10 min; TDE1 stated to work as well |
| 8 | Strip Tn5 with 0.5 µL 0.2% SDS | 5 min |
| 9 | Add custom S50X/N70X index primers (0.1 µM each) and Phusion HF mix | |
| 10 | **72 °C 3 min gap fill**, then index PCR | 98 °C 3 min; 12×(98 °C 10 s, 55 °C 30 s, 72 °C 30 s); 72 °C 5 min |
| 11 | Pool wells; 0.6:1 PEG-bead cleanup | optional E-Gel size selection 550 bp–2 kb not drawn |

## 4. Two kinds of fragment

This is the point of the design, and it is why the schematic draws the final structure
twice.

**5' fragments** retained the TSO, so they carry the 11-bp tag and the UMI:

```
P5 [i5] s5 ME ATTGCGCAATG [8-bp UMI] GGG --- cDNA --- revcomp(ME) revcomp(s7) [i7] revcomp(P7)
```

**Internal fragments** did not, and are **character-for-character identical to a SMART-seq2
fragment** (asserted in `tools/selftest.py`). They carry no UMI and no positional information,
and are used for coverage only.

🟡 So the demultiplexer's job is simply: does read 1 begin with the 11-bp tag? If yes, the
read pair is 5'-derived and the following 8 nt are its UMI.

Single cells are identified by the **i5 + i7 index combination**, as in SMART-seq2 — there is
still no in-molecule cell barcode.

## 5. Back end

Identical to SMART-seq2 and to every other Nextera method — see
[`../ref/concepts/tn5-tagmentation.md`](../ref/concepts/tn5-tagmentation.md). The sequencing
primers are unchanged; only read 1's *content* differs, depending on whether that fragment
kept the TSO.


## Source recovery and index-design scope (2026-10-16)

The original publisher HTML was retrieved through `tools/get_sources.py`, but its
body explicitly labels itself a subscription preview (`s41587-020-0497-0.html.txt`,
line 627); it does not supply the published Methods. The PDF route also returns
HTML rather than the original PDF. These are concrete access restrictions; the
preprint and author protocol remain the accessible method sources.

The original author protocol is integer ID 33276, DOI
`10.17504/protocols.io.bcq4ivyw`. Its anonymous public API guidelines and its PDF
pages 3–4 give custom index-design instructions, including ten-base or eight-base
options, minimum Levenshtein distance, homopolymer/GC filtering and avoidance of
`AC` at the i5 start. This original protocol does **not identify a separate custom
index-order workbook**. The previous review's assumed unavailable Smart3 i5/i7
order list was not an identified source and has been corrected in the checklist.
The distinct xpress protocol's recovered 384-well workbook cannot establish the
original Smart3 experiment's selected pair.

The actual linked source is Costello et al., BMC Genomics 2018,
DOI `10.1186/s12864-018-4703-0` (PMC5941783). Its Additional file 3 explains the
`AC` mispriming issue and general index-design recommendations. Additional file 2
is that paper's own two sets of 96 eight-base indexes, not a Smart3 order roster;
Additional file 1 concerns index swapping. These linked originals were checked
for their actual scope. No categorical absence of a selected Smart3 pair is
claimed while the published Methods remain inaccessible.


## Bounded pass 2026-10-17 (round breadth-first-18)

Read: Supplementary Table 1 (`MOESM3` xlsx, all 96 condition rows; final Smart-seq3
conditions carry the four oligos above), Supplementary Note I (`MOESM1`, design and
additive rationale), author protocol V.3 PDF pp2–20 (oligo table, index design,
lysis, RT, pre-amplification, bead cleanup, tagmentation, index PCR, cleanup,
sequencing/bcl2fastq example `Y150N,I8,I8,Y150N`). The publisher HTML (line 627)
remains a subscription preview. The schematic graph is now declared in
`tools/build_page.py` with these reagents, conditions, both PEG-bead cleanups and
the SDS strip; molecular states still come from the shared SMART-seq helpers.
Open: the published Methods text; the acquired cycle counts of the published runs
(the bcl2fastq line is a protocol example, not a run record); which index length
and members each experiment used.
