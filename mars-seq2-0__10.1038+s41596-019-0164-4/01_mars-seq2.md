# MARS-seq2.0 — source note

## Active native migration, 2026-10-18

Source scope: Published2019NatureProtocols main s41596-019-0164-4.pdf.txt reagentsetup637–690, capture714–740, RT974–1043, Exo/pooling1057–1123, second-strand/IVT1170–1237, fragmentation/ligation1255–1322, RT2/PCR1324–1422, acquisition1448–1461 (Step62). Original SupplementTable1 scg_41596_2019_164_MOESM4_ESM.xlsx.txt header/A1–F1 primer rows; original SupplementTable2 scg_41596_2019_164_MOESM5_ESM.xlsx.txt header/ix1–ix10 donor/readbarcode rows. Mirrored files are publisher-original tables with matching manifest identity; protocols.io andMOESM1 not newly reread.

MARS-seq2.0 now owns its17step graph and source-specific molecules in tools/chem.py, with no2014constructor import. Published main reagent setup itself verifies RT2 and both enrichment primers. Actual RT1ssDNA undergoes ExoI before pooling, then0.9×SPRI; native second-strand synthesis,1.4×SPRI/on-bead IVT, TurboDNase,1.2×SPRI, Zn fragmentation/2×SPRI, native5P/C3blocked N5/pool4 donor ligation/1.5×SPRI, RT2/1.5×SPRI, KAPA PCR and0.7×SPRI. T7 initialGGG is restored, RNA U distinguished from DNA T, and independent originals/copies/hydrolysates carried explicitly with unresolved partition. Actual recovered ssDNA feeds PCR. Representative ix1 pool identity is4bases, N5 are adapter random bases, well7/UMI8 are separate features. Published Step62 acquisition75/15/noindices is distinct from physical molecular lengths. Historical MOESM1N3read-design discrepancy remains public and unresolved; its source was not newly reread.

This capped pass establishes no confirmed absence or complete chemistry review. Historical observations and audits below retain their original scope. Current source-owned declarations replace the earlier shared endpoint and IVT recipes; molecular hydrolysis, nonamplification and computational filtering are distinct from physical discard.

Deferred: Complete sequential first-cycle PCR geometry on actual native ssDNA; loaded sequencing-primer identities; Historical MOESM1read-design example3I.4P.50M/7W.8R.5I versus native N5adapter/75–15run; Exact molecular partition/integrity of original RNA and promoterDNA/hydrolysate through known recoveries; Other RT1/pool adapter members and robot details outside the selected representative molecular view.

Mandatory scraper: 61hits7files; main native primer strings and original RT1A1/adapterix1 modifications verified in inspected originals. No2014orupstream endpoint used as authority.


The full source audit, oligo transcription and step-by-step reading are retained in the
companion [MARS-seq / MARS-seq2.0 research note](../mars-seq__10.1126+science.1247651/01_mars-seq.md),
as historical context. The active2019 graph and molecular declarations now live independently in this protocol’s tools/chem.py.

**Defining source:** Keren-Shaul H *et al.* “MARS-seq2.0: an experimental and analytical
pipeline for indexed sorting combined with single-cell RNA sequencing.” *Nature
Protocols* 14, 1841–1862 (2019).
doi:[10.1038/s41596-019-0164-4](https://doi.org/10.1038/s41596-019-0164-4).

## Where the schematic's bases come from

> **The scg_lib_structs drawn page is not a source for this protocol.** It is inspiration
> only: no sequence, length or reaction below is taken from it, and agreement with it is
> not verification. Its copy on disk is `upstream_MARS-seq.html`, and `scg_adapters_MARS_seq.fa`
> is that page's own derived FASTA; neither is cited here.

Every concrete base the MARS-seq2.0 schematic draws is either a canonical `lib/illumina.py`
constant or a transcription from the authors' own supplementary tables and detailed
[protocols.io procedure](https://doi.org/10.17504/protocols.io.7hkhj4w).

| Modelled material | Bases | Authoritative document |
|---|---|---|
| RT1: `CGATTGAGGCCGG` · T7 promoter · `GC` · `GACGTGTGCTCTTCCGATCT` · 7-nt well barcode · N8 · T20 · N | 🟢 verbatim | **MARS-seq2.0 Supplementary Table 1**, sheet `MARS-seq2.0_RT1 primers` — all 384 rows match this form |
| pool ligation adapter: `/5Phos/` · 4-nt pool barcode · N5 · `AGATCGGAAGAGCGTCGTGTAG` · `/3SpC3/` | 🟢 verbatim | **MARS-seq2.0 Supplementary Table 2**, sheet `MARS-seq2.0_Ligation adapter` — all 32 rows |
| RT2: `CTACACGACGCTCTTCCGATCT` | 🟢 verbatim | **Keren-Shaul et al., protocols.io (2019), reagent setup p. 29**, “Second RT primer (RT2)”; local PDF text line 1416 |
| P5-Read1 enrichment primer: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` | 🟢 verbatim | **Author protocols.io procedure, reagent setup p. 29**, “Enrichment primers”; local PDF text line 1420 |
| P7-Read2 enrichment primer: `CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` | 🟢 verbatim | **Author protocols.io procedure, reagent setup p. 29**, “Enrichment primers”; local PDF text line 1421 |
| read design `3I.4P.50M` / `7W.8R.5I` | 🟢 verbatim | MARS-seq2.0 Supplementary Information (`MOESM1`), pipeline `seq_batches.txt` field description |
| MARS-seq (v1) RT1, ligation adapters `lig_NNNX4_ix1..8`, second RT primer `TCTAGCCTTCTCGCAGCACATC`, `P5_Rd1 PCR`, `P7_Rd2 PCR` | 🟢 verbatim | **Jaitin *et al.* 2014, *Science* 343:776, Supplementary Material**, Tables S7–S8 — the predecessor this protocol versions; not used by the 2.0 graph. The printed v1 RT2 is the exact reversal of the 2.0 RT2 (an apparent 3′→5′ listing, shown on the 2014 page) |
| P5, P7, TruSeq Read 1 / Read 2 arms | 🟢 canonical | `lib/illumina.py` (Illumina adapter document), imported, never retyped |

Supplementary Tables 1 and 2 are the publisher's own ESM files
`41596_2019_164_MOESM4_ESM.xlsx` / `MOESM5_ESM.xlsx`. They are on disk in two copies —
fetched from `media.springernature.com` and mirrored in the scg_lib_structs repository
under `data/MARS-seq/` — and the manifest records them as **md5-identical**
(`9e4f86ee2a009492cb706cc8be6906ee`, `e584f05dc8ca8f0593e6bfd28f5bacba`). A filename
prefix is not provenance: these are Springer originals regardless of which copy is read.
The same holds for `scg_jaitin-sm.pdf`, which is the Science supplement's own PDF (its
first page reads “Supplementary Material for … Published 14 February 2014, Science 343,
776 (2014), DOI: 10.1126/science.1247651”).

## Remaining source limits

🟢 The Nature Protocols main text was not available in the local source bundle, but the
same authors published a [detailed 54-page laboratory procedure](https://www.protocols.io/view/mars-seq2-0-an-experimental-and-analytical-pipelin-7hkhj4w.pdf)
on protocols.io. Its reagent setup, p. 29, prints RT2 and both enrichment primers, and
the procedure gives reaction conditions. The original PDF is fetched through the declared
`protocolsio:` source route, and its `.txt` twin was checked with the oligo scraper.
The previous conclusion
that these three oligos were unpublished was wrong because it checked only the paper
supplements.

🟡 The final library's P5 and P7 flanks use the v2 enrichment primers. The pre-PCR
template is derived from the published RT2 and RT1 handles, and both enrichment primers
are checked for inward-facing annealing before the full library is drawn. The authors do
not separately name the instrument sequencing-primer oligos; TruSeq read primer landing
sites are computed from the published library flanks.

## Bounded pass 2026-10-17 (round 12): author procedure Steps 31–80

Read in the protocols.io PDF text (its Steps, Box 1 and Reagent setup follow the
*Nature Protocols* numbering; the paywalled journal original itself was not compared):

- 🟢 Steps 29–35: RT 2 min 42 °C, 50 min 50 °C, 5 min 85 °C; Exo I added to each
  384-well (30 min 37 °C, 10 min 80 °C) **before** the two-step robotic pooling
  (Steps 36–39) and 0.9× AMPure cleanup (Step 40, QC1). The graph's Exo I-then-pool
  order already matched.
- 🟢 Steps 42–49: second strand (NEBNext module, enzyme 1:8, 2 h 16 °C), 1.4× AMPure
  with IVT on beads, T7 IVT overnight 37 °C, Turbo DNase I 15 min 37 °C, 1.2× AMPure.
- 🟢 Steps 51–55: Zn²⁺ fragmentation 3 min 70 °C (median ~300 nt), stop solution, 2× AMPure.
- 🟢 Steps 56–61: plate adapter + RNA 3 min 70 °C, T4 RNA ligase 1 with DMSO, ATP and PEG,
  2 h 22 °C, 1.5× AMPure.
- 🟢 Steps 62–67: RT2 "complementary to the ligated adaptor" annealed 3 min 72 °C;
  AffinityScript 2 min 42 °C, 50 min 45 °C, 5 min 85 °C; 1.5× AMPure, QC2.
- 🟢 Steps 69–72: KAPA HiFi HotStart ReadyMix + enrichment primer mix; example 2 cycles
  at 57 °C + 12 at 65 °C (typically 11–17); 0.7× AMPure, QC3.
- 🟢 Steps 79–80: 20–40 plates per NextSeq 500 High Output v2 flow cell;
  **rd1 = 75 bp, rd2 = 15 bp, no indices.**

Model changes: named reagents/conditions on the fragmentation, ligation, RT2 and PCR
inputs; acquired read lengths 75/15 declared in the sequencing plan; Read 2 drawn as
well barcode + UMI only. The MOESM1 `seq_batches.txt` example (3I.4P.50M / 7W.8R.5I,
57 + 20 cycles) is shown on the page as a different configuration, not the Step 80 run.
The AMPure cleanups are not yet separate graph steps (deferred).
