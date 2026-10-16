# LIANTI — linear whole-genome amplification via a T7-promoter transposon

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/`; the model is `tools/lianti.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length, kit or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence. **No LIANTI oligo is published anywhere we can
> read** 🔴: the paper's main text prints none, and the Supplementary Materials and
> Methods (.docx) sits behind the PMC download gate. Earlier revisions of this note and
> of the model carried the drawn page's reconstruction — a 10-nt spacer `GAACAGAATT`, an
> 8-nt UMI, a 6-nt i7, and the NEBNext Ultra hairpin-adaptor/USER library prep. None of
> that has a source, so none of it is asserted any more; the model now draws placeholders
> and an unspecified indexed Illumina adapter instead.

**LIANTI** — Chen C, Xing D, Tan L, Li H, Zhou G, Huang L, Xie XS. "Single-cell whole-genome
analyses by Linear Amplification via Transposon Insertion (LIANTI)." *Science* 2017;
356:189–194. doi:[10.1126/science.aak9787](https://doi.org/10.1126/science.aak9787)
(PMID 28408603, PMC5538131). Data: SRA SRP102259. No other papers for this method are
listed in `catalogue/scg_lib_structs.tsv`.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/lianti__10.1126+science.aak9787/`, never committed):

| File | What | Used for |
|---|---|---|
| `science.aak9787_PMC5538131.html(.txt)` | PMC author manuscript, full main text + figure legends | concept, transposon architecture (Fig. 1B/C legend), yields, coverage |
| `upstream_LIANTI.html(.txt)` | scg_lib_structs method page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `lib/nextera.py`, `lib/illumina.py` (in-repo, from Illumina's adapter document) | canonical Tn5 mosaic end; canonical single-index Illumina arms | the only real bases the model draws, each one created `inferred=True` |

Not fetched (needs a manual download):

| URL | What | Why it matters |
|---|---|---|
| https://pmc.ncbi.nlm.nih.gov/articles/instance/5538131/bin/NIHMS883352-supplement-Supplementary.docx | Supplementary Materials (6.1 MB docx; PMC download gate) | Materials and Methods: the transposon and second-strand primer **as published**, buffers, temperatures, IVT time, RT enzyme, RNase, library kit, sequencer |
| https://www.science.org/doi/10.1126/science.aak9787 | publisher version + supplement | same supplement, alternative route |

---

## 1. What it is

A single-cell **whole-genome amplification (WGA)** method that replaces exponential
amplification (PCR, MDA, MALBAC) with **linear** amplification by **T7 in vitro
transcription (IVT)**. 🟢 (main text) Tn5 fragments a single cell's genome and in the same
event attaches a **T7 promoter** to every fragment end; IVT then makes thousands of RNA
copies from the original fragments, never from copies of copies, so errors and size bias
do not compound. 🟢 (main text, Fig. 1A argument)

Builds on [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) and
[reverse transcription](../ref/concepts/reverse-transcription.md).

| | New thing here | Consequence |
|---|---|---|
| 1 | **One-oligo hairpin transposon**: a 19-bp double-stranded Tn5 binding site closed by a **single-stranded T7-promoter loop** 🟢 (Fig. 1B legend) | the Tn5 homodimer carries the same adaptor on both monomers 🟡 *(computed from the one-oligo design)*, so **every** cut is amplifiable, not only the s5/s7 hetero-ended fragments of ordinary Nextera 🟡 |
| 2 | **Linear amplification by IVT** of the tagmented genome 🟢 | CV of read depth lowest of the WGA methods compared at all bin sizes (Fig. 1E) 🟢 |
| 3 | **RNA self-priming at its 3' end** for RT 🟢 (Fig. 1C legend) | no random or degenerate primer anywhere before library prep |
| 4 | **UMI in the second-strand primer** 🟢 ("tagged with unique molecular barcodes", Fig. 1C legend) plus "same ends = same fragment" digital counting (Fig. 2A) 🟢 | micro-CNV detection at ~10 kb 🟢 |

No cell barcode: one cell per tube, identity carried by the library index. 🟡 The
authors note that cell barcodes could be added to the transposon and primer. 🟢

## 2. Oligos — 🔴 none is published

The paper describes the transposon in words and draws it as a cartoon; it prints no
sequence, and its oligo table is in the supplement we cannot read. What it does fix:

> "LIANTI transposon consists of a 19-bp double-stranded transposase binding site and a
> single-stranded T7 promoter loop. Equal molar of LIANTI transposon and Tn5 transposase
> are mixed and dimerized to form LIANTI transposome." — Fig. 1B legend 🟢
> ("… double-stranded LIANTI amplicons tagged with unique molecular barcodes … ready for
> DNA library preparation and next generation sequencing" — Fig. 1C legend 🟢)

### What the model draws, and why each piece is allowed

| Segment | Drawn as | Basis |
|---|---|---|
| mosaic end / its reverse complement | `nextera.ME`, `nextera.ME_RC` | the paper's "19-bp double-stranded transposase binding site" is the Tn5 mosaic end; the bases are canonical 🟡 (identification), 🟢 (bases, Illumina adapter document) |
| T7 promoter | `TAATACGACTCACTATA` | canonical T7 class III promoter for the paper's "T7 promoter" 🟡 |
| T7 start | `GGG` | the +1…+3 of that promoter; the paper's transcripts must begin with them 🟡 *(computed)* |
| **spacer between ME-rc and the promoter** | `N` × 10 placeholder, `inferred=True` | 🔴 **neither bases nor length published.** The 10 drawn positions are a drawing extent, not a measurement — the loop must merely be long enough to be single-stranded |
| UMI | `N` × 8 placeholder | role 🟢 ("unique molecular barcodes"), width 🔴 |
| i7 index | `I` × 6 placeholder | that the library was indexed 🟡, width 🔴 |
| library arms (P5, Read 1, index-read site, P7′) | `illumina.P5`, `TRUSEQ_READ1`, `INDEX1_PRIMER`, `P7_RC`, all `inferred=True` | 🔴 the paper names **no** library-prep kit, adapter or primer; a single-index TruSeq-style shell is the model's assumption and is drawn entirely dotted |

### How the transposon folds — 🟡 *(computed)*

The transposon's total length is **not** derivable, because the spacer's is not; the
order is:

| Segment | Length | Identity |
|---|---|---|
| `CTGTCTCTTATACACATCT` | 19 | reverse complement of the mosaic end (`nextera.ME_RC`) 🟡 |
| `NNNN…` | 🔴 unknown (drawn as 10) | the unverified spacer 🔴 |
| `TAATACGACTCACTATA` | 17 | T7 promoter (class III consensus, −17 to −1) 🟡 |
| `GGG` | 3 | T7 +1 to +3; first transcribed bases 🟡 |
| `AGATGTGTATAAGAGACAG` | 19 | mosaic end (ME), 3′ end = the transferred end 🟡 |

- ME pairs with ME-rc: a **19-bp double-stranded Tn5 binding site**, with the spacer,
  promoter and `GGG` forming the **single-stranded loop** — exactly the paper's
  description 🟢, assembled from canonical parts 🟡. `tools/lianti.py` builds this with
  `HairpinTopology`, so the stem is checked rather than drawn by hand.
- The promoter points **into** the genomic insert from the ME side 🟡; there is no s5/s7
  arm at all — the Nextera entry points are replaced by the loop 🟡.
- The **second-strand primer = UMI + `GGG` + ME**, i.e. the 5′ end of every transcript, so
  it anneals to the 3′ end of every first-strand cDNA 🟡 *(computed)*.
- Whether the transposon's 5′ end is phosphorylated, and whether the hairpin is joined to
  the genome after gap fill or displaced, is not stated 🔴.

### Library prep — 🔴 unestablished

The paper ends at "ready for DNA library preparation and next generation sequencing". The available main article names no kit, adapter, primer, cycle count, read length or index set, and the supplement is unavailable. The model therefore uses an **unspecified indexed Illumina adapter** with
`AdapterLigationGeometry.unavailable(...)`, a library PCR with `PCRGeometry.unavailable(...)`,
and a library whose every adapter-derived base is inferred. Ligation-based prep (rather
than tagmentation) is the model's one structural assumption, on the ground that tagmenting
the amplicon would destroy the UMI end that the method exists to install 🟡.

## 3. Step by step

Conditions in the paper's main text are sparse; all buffer, time and temperature
details are in the supplement we could not fetch. 🔴

1. **Transposome assembly**: equimolar LIANTI transposon and Tn5 transposase, dimerised
   🟢 (Fig. 1B legend). Both monomers carry the identical hairpin adaptor (homodimer). 🟡
2. **Single cell picking and lysis** (BJ fibroblasts; how cells were collected is in the
   unavailable supplement 🔴). Optional **UDG** treatment of the lysate removes deaminated
   cytosines (C→U), the source of the C→T false positives 🟢 (Fig. 3B, main text).
3. **Tagmentation**: random fragmentation, average ~400 bp 🟢. Each fragment end carries
   the transposon's ME (3' end joined to the genomic 5' end), the hairpin and the usual
   9-nt gap on the opposite strand 🟡 ([Tn5](../ref/concepts/tn5-tagmentation.md)).
4. **Gap extension** with a DNA polymerase (not named in the main text 🔴) converts the
   single-stranded T7 loops into **double-stranded T7 promoters on both ends** 🟢 (Fig. 1C
   legend). Computed product (top strand, 🟡):
   transposon · insert · reverse complement of the transposon.
5. **In vitro transcription, overnight** 🟢. Both ends have a promoter pointing inward, so
   both strands of the fragment are transcribed. Computed transcript (5'→3', 🟡):
   `GGG` · ME · <insert> · ME-rc · `CCC` · promoter reverse complement · spacer reverse
   complement (unverified) · ME — 22 nt before the insert. What happens where the two
   convergent transcription events meet is not addressed by the paper 🔴.
6. **Self-primed reverse transcription** 🟢 ("capable of self-priming on the 3' end").
   Computed mechanism 🟡: the transcript's 3'-terminal ME pairs with its internal ME-rc
   copy (the same 19-bp stem as the transposon, now in RNA, with the unverified
   spacer/promoter interval as the turn), and the 3' OH primes RT back across the insert
   to the transcript's 5' end. First strand
   (after RNA removal): rc(insert) · `CTGTCTCTTATACACATCT` · `CCC`.
7. **RNase digestion** removes the RNA 🟢 (enzyme not named in the main text 🔴).
8. **Second-strand synthesis** with the UMI primer 🟢 (concept) / 🔴 (sequence) →
   double-stranded LIANTI amplicon: `<UMI> · GGG · ME · <insert>`, with the far end bare
   genomic sequence (the other Tn5 insertion point). Yield ~20 ng per cell 🟢.
9. **Library prep**: 🔴 the paper says only "ready for DNA library preparation and next
   generation sequencing". No kit, adapter, primer or cycle count is named in the available main article. The model
   draws end repair → dA-tailing → ligation of an *unspecified indexed Illumina adapter*
   (ligation geometry `unavailable`) → PCR (geometry `unavailable`).

## 4. Final library — 🟡 assumed shell, 🔴 no published base

Both orientations occur, because the dA-tailed amplicon is ligated on both ends. The
shell below is the model's assumption (a single-index Illumina library), not a published
layout; on the page every one of its bases is dotted.

Orientation 1 (UMI on the P5 side):

```
5'- P5 · TruSeq Read 1 · <UMI> · GGG · ME · <genomic insert> · A · index-read site · <i7> · P7' -3'
```

Orientation 2 (UMI on the P7 side):

```
5'- P5 · TruSeq Read 1 · <genomic insert> · ME-rc · CCC · <UMI, rc> · A · index-read site · <i7> · P7' -3'
```

## 5. Read layout / sequencing

All 🟡 — consequences of the assumed shell, not source-verified read-layout facts. 🔴 No read length,
index width or sequencer is given anywhere we hold.

- Read 1 (TruSeq Read 1 primer): in orientation 1, UMI + `GGG` + 19-nt ME before the
  genome; in orientation 2, genome from cycle 1. 🟡 *(computed)*
- Index 1: the i7 index, drawn 6 nt wide as a placeholder 🔴.
- Read 2 (TruSeq Read 2 primer): the mirror of Read 1. 🟡
- Both read starts are Tn5 insertion points, which is what makes "same start and end =
  same original fragment" digital counting possible (Fig. 2A) 🟢 / 🟡.
- Platform, read lengths: not in the main text 🔴. BJ cells sequenced to ~30× 🟢; 97 %
  genome coverage, 17 % ADO 🟢 (Table S1, as quoted in main text).

## 6. What the paper fixes and what it leaves open

| Point | Science 2017 main text / figure legends | Status in the model |
|---|---|---|
| Transposon architecture | 19-bp ds transposase binding site + ss T7 promoter loop 🟢 | `HairpinTopology` stem of canonical ME; loop = placeholder spacer + canonical promoter 🟡 |
| Transposome | equimolar transposon + Tn5, dimerised 🟢 | homodimer; both monomers carry the same adaptor 🟡 |
| Gap extension | converts the ss loops into ds T7 promoters on both ends 🟢 | drawn, all bases inferred |
| IVT | overnight, linear 🟢 | drawn; the "one strand dominates" idea has no source and is not asserted |
| Self-primed RT, RNase, second strand with UMI 🟢 | drawn | widths placeholders |
| Transposon spacer, UMI width, i7 width | not given 🔴 | placeholders |
| Library prep, PCR primers, cycles, sequencer, read lengths | not given 🔴 | unavailable geometries, inferred shell |

Detailed conditions may be in the unavailable supplement 🔴; the available main article does not establish them, and nothing in the model depends on them.

## 7. Open questions

- 🔴 The supplement's exact oligo sequences (the spacer and its length, the UMI width,
  modifications) and the transposon annealing protocol. The Science supplement
  (`https://www.science.org/doi/10.1126/science.aak9787`) and the PMC author-manuscript
  DOCX are both gated; the manifest carries both URLs for a manual download. A second
  lead: the acknowledgments disclose a Harvard patent application covering LIANTI 🟢,
  and a granted patent would print the transposon — the application number is not given
  in the paper, so no source row can be added for it yet.
- 🔴 Role of the spacer, and of a transposon 5′ phosphate (is the hairpin ligated to the
  genome after gap fill, or displaced?).
- 🔴 Which RT, which RNase, and whether IVT also yields RNA-templated RNA products (refs
  51–53 in the paper concern T7 RNA-dependent synthesis, suggesting the supplement
  discusses it).
- 🔴 Sequencer, read lengths, PCR cycles.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text but not the supplement (.docx behind
PMC's download gate; listed as `(manual)` in MANIFEST.tsv). `tools/scrape_primers.py`
found **no oligo in any source that counts** — the paper contains none — which is the
whole finding for this protocol. Transposon architecture and every derived structure were
computed with `lib/` (`nextera.ME`, `illumina.*`, `chemdraw.revcomp`,
`topology.HairpinTopology`).
