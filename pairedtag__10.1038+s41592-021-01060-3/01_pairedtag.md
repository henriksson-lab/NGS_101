# PairedTag

## 1. What it is

PairedTag jointly profiles a chosen histone modification and the transcriptome from the same single nuclei using antibody-targeted tagmentation and shared combinatorial cell barcodes.

## 2. Authoritative source

🟢 Zhu et al., *Joint profiling of histone modifications and transcriptome in single cells from mouse brain*, Nature Methods (2021), DOI `10.1038/s41592-021-01060-3`.

The PMC full text, Extended Data Fig. 1 and Supplementary Table 1 (primer and all
12 × 96 × 96 barcode sequences) were read. The 2026-10-17 bounded pass re-read the
author-manuscript JATS (PMC7954905) Online Methods text lines 116–156 (adaptor annealing,
transposome assembly, tagmentation, RT, ligation barcoding, pre-amplification, digestion
and second-adapter tagging, indexing PCR and sequencing), Extended Data Fig. 2a legend
(line 254), all four sheets of Supplementary Table 1 and the Supplementary Notes PDF
(no chemistry). The step-by-step Protocol Exchange protocol (ref. 51,
doi:10.21203/rs.3.pex-1301/v1) is behind Cloudflare and was not read; the Paired-seq
predecessor (ref. 31) was not followed.

## 3. Construct-changing operations

1. 🟢 Antibody-bound pA–Tn5 transfers a phosphorylated adapter containing the common
   15-nt ligation arm, a modality base, four-base BC1, NotI and the 19-nt mosaic end.
   pMENTs supplies the paired mosaic-end strand and is 3′-ddC blocked.
2. 🟢 In the same first-round well, matched phosphorylated oligo(dT)16-VN and random-N6
   primers carry the common arm, RNA modality base, the same four-base BC1 and SbfI.
3. 🟢 Nuclei are pooled and split across 96 BC2 oligos. Linker-R02 bridges the BC1 arm
   to the BC2 oligo's ligating arm; Blocker-R02 is its exact reverse complement. The
   nuclei are pooled again and split across 96 R04 oligos carrying Read 2, a 10-nt UMI,
   seven-base BC3, a 0–3-nt phase block and the final ligating arm. Linker-R03 provides
   the splint and Quencher-R03 is its exact reverse complement.
4. 🟢 Ligation rounds use T4 DNA ligase (NEB M0202L); Blocker-R02 is added after round 2
   and an EDTA Termination-Solution containing the "R04 Terminator oligo" (Methods name;
   Supplementary Table 1 lists Quencher-R03) after round 3. After proteinase K lysis and
   1× SPRI, the heat-denatured DNA/cDNA is dC-tailed by NEB TdT (dCTP only); KAPA HiFi HS
   performs 16 linear Anchor-FokI-GH cycles and then 11 PA-F/PA-R cycles, followed by a
   0.2×/0.65× SPRI double size selection.
5. 🟢 The amplified material is split. DNA is cut with SbfI/FokI and receives an annealed
   P5-FokI/P5c-NNDC-FokI adapter. The source-defined phosphorylated NNDC end is the
   four-base cohesive ligation end; its complement depends on the FokI-cut molecule.
   RNA-derived material is cut with NotI and receives an N5 end by Tn5. 🟡 SbfI sites occur
   only in the RNA BC1 oligos and NotI sites only in the DNA BC1 oligos, so each digestion
   is expected to cleave the other modality; the article does not state this purpose.
6. 🟢 Both final libraries use Q5 indexing PCR (72 °C 5 min first, 9–14 cycles) with P5
   Universal (DNA) or N5 (RNA) plus P7, and were sequenced on NextSeq 550 or NovaSeq 6000 as
   PE 50 + 7 + 100 or PE 100 + 7 + 100 (Read 1 + Index 1 + Read 2). The model draws the
   100-cycle Read 1 option.

## 4. Key published oligos

- `pMENTs`: `/5Phos/CTGTCTCTTATACACATCTddC`
- `AdaptorA`: `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`
- `Linker-R02`: `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`
- `Blocker-R02`: `ATCCACGTGCTTGAGAGGCCAGAGCATTCG`
- `Linker-R03`: `GGTCTGAGTTCGCACCGAAACATCGGCCAC`
- `Quencher-R03`: `GTGGCCGATGTTTCGGTGCGAACTCAGACC`
- `P5-FokI`: `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- `P5c-NNDC-FokI`: `/5Phos/NNDCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG`
- `Anchor-FokI-GH`: `AAGCAGTGGTATCAACGCAGAGTGAAGGATGTGGGGGGGGG*H` (the table does not
  define `*`; the Methods call this oligo "Anchor-FokI-GSH-Oligo")
- `PA-F`: `CAGACGTGTGCTCTTCCGATCT`
- `PA-R`: `AAGCAGTGGTATCAACGCAGAGT`
- `P5 Universal`: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T`
  (printed with `*` before the final T; read as a phosphorothioate by convention, 🟡).

## 4a. Published discrepancies (shown on the page)

- 🟢 Online Methods (Reverse transcription) names a "barcoded T15 primer", but every
  RNA_#NN_RE member of Supplementary Table 1 (BC Plate #01 rows B1–B12) prints
  SbfI–T16–VN. The model uses the table sequence; whether "T15" is a loose name is open.
- 🟢 The Methods name "Anchor-FokI-GSH-Oligo" differs from the table entry
  "Anchor-FokI-GH". The table sequence is drawn with its printed `*`; the meaning of `*`
  and of "S" is not defined by the source.

The plate contains all 12 × 96 × 96 barcode choices; the schematic uses role-labelled
four-base BC1 and seven-base BC2/BC3 representatives rather than presenting one cell
identity as universal. A zero-phase R04 member is drawn; other wells contain one to three
phase bases at the documented position.

## 5. Final read layouts

Read 1 maps the histone-mark DNA or transcript cDNA. Read 2 traverses the 10-base UMI and
the three nuclear barcodes; the paper places the first bases of BC3, BC2 and BC1 around
positions 10–13, 47–50 and 84–87, respectively. Index 1 contains the PCR sample index.
The final model derives representative starts at cycles 11, 48 and 86 directly from the
source oligos. The six-base i7 is acquired in a seven-cycle index read; the seventh cycle
continues into the adjacent P7 sequence. The RNA final panel shows the oligo(dT)-primed
branch; random-N6-primed molecules replace that primer-derived tract.

## 6. Open questions

- The Protocol Exchange step-by-step protocol (ref. 51) is unread; it may define the
  anchor `*`, the T15/T16 naming and the R04 Terminator identity.
- The purpose of SbfI/NotI cross-modality digestion is inferred from site placement.
- Which Read 1 length (50 or 100) applies to each published dataset is not assigned.

## Native graph migration — 2026-10-18 batch04

This pass inspected three originals: the published article Methods, Supplementary Table1, and the MOESM13 numerical-data heading/initial rows (classified as data, no chemistry imported). Local native graphs now retain three independent molecular populations: histone-mark DNA, dT-primed cDNA and N6-primed cDNA. The physical pMENTs terminal ddC stays on the nontransferred oligo; an exact whole-modified-source regression covers the documented parser exemption. The printed Anchor9G tract does not establish a9C TdT tail. Tail/register/extension boundaries are unassigned while source barcode/UMI lineage remains attached. Actual Anchor-extension material feeds PA-PCR; source aliquots retain mixtures, and cleavage does not itself purify modalities. Source cleanups, sublibrary choices and the1µl aliquot remainder are explicit. Final P5/Tn5-N5 library references are conditional and cannot supply their own PCR substrates. Both published50+7+100 and100+7+100 layouts remain; the100cycle example is drawn, and index2 acquisition stays unassigned. Existing dT15-versus-tableT16 and AnchorGSH/name/*mark discrepancies remain public. Source* is retained as printed; its exact chemical meaning is not assigned.

The local curator-written review target tuple is independent of graph contents. This capped pass does not complete the missing-information audit and accepts no source-absence claim.
