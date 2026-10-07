# LIANTI — linear whole-genome amplification via a T7-promoter transposon

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available to us. Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note; they are not yet asserted in a self-test (no
> `tools/` module for this protocol). **Every oligo sequence below is 🟡**: the paper's
> main text contains no sequences, the Supplementary Materials (.docx) could not be
> fetched, and the only sequence source we hold is the scg_lib_structs page (a secondary
> source).

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
| `upstream_LIANTI.html(.txt)` | scg_lib_structs method page (Teichmann lab) | **all oligo sequences**, step-by-step, NEB library prep, read layout |

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
| 1 | **One-oligo hairpin transposon**: a 19-bp double-stranded Tn5 binding site closed by a **single-stranded T7-promoter loop** 🟢 (Fig. 1B legend) | the Tn5 homodimer carries the same adaptor on both monomers, so **every** cut is amplifiable, not only s5/s7 hetero-ended fragments (upstream page) 🟡 |
| 2 | **Linear amplification by IVT** of the tagmented genome 🟢 | CV of read depth lowest of the WGA methods compared at all bin sizes (Fig. 1E) 🟢 |
| 3 | **RNA self-priming at its 3' end** for RT 🟢 (Fig. 1C legend) | no random or degenerate primer anywhere before library prep |
| 4 | **UMI in the second-strand primer** 🟢 ("tagged with unique molecular barcodes", Fig. 1C legend) plus "same ends = same fragment" digital counting (Fig. 2A) 🟢 | micro-CNV detection at ~10 kb 🟢 |

No cell barcode: one cell per tube, identity carried by the library index. 🟡 The
authors note that cell barcodes could be added to the transposon and primer. 🟢

## 2. Oligos

🟡 As written on the scg_lib_structs page (secondary; not checkable against the
supplement, which we do not have). Modifications as the page gives them: `/Phos/`,
`/5Phos/` = 5' phosphate, `(dU)` = deoxyuridine, `*` = phosphorothioate.

```
LIANTI transposon DNA    5'- /Phos/CTGTCTCTTATACACATCTGAACAGAATTTAATACGACTCACTATAGGGAGATGTGTATAAGAGACAG -3'
Second strand primer     5'- [8-bp UMI]GGGAGATGTGTATAAGAGACAG -3'

NEBNext Hairpin Adaptor  5'- /5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC(dU)ACACTCTTTCCCTACACGACGCTCTTCCGATC*T -3'
NEBNext Universal PCR Primer for Illumina
                         5'- AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'
NEBNext Index 1 Primer for Illumina
                         5'- CAAGCAGAAGACGGCATACGAGAT[i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -3'

Illumina TruSeq Read 1 primer   5'- ACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'
Illumina TruSeq Read 2 primer   5'- GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -3'
Sample index sequencing primer  5'- GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -3'
```

The NEBNext oligos are commercial kit components (NEBNext Ultra DNA Library Prep, with
Index Primers); that the paper used this kit is the upstream page's statement, "as
suggested from the LIANTI paper". 🟡

### How the transposon is built — 🟡 (computed)

The 68-nt transposon is one oligo that folds on itself:

| nt | Segment | Identity |
|---|---|---|
| 1–19 | `CTGTCTCTTATACACATCT` | **reverse complement of the mosaic end** (`nextera.ME`) |
| 20–29 | `GAACAGAATT` | 10-nt spacer, function not given 🔴 |
| 30–46 | `TAATACGACTCACTATA` | **T7 promoter** (class III consensus, −17 to −1) |
| 47–49 | `GGG` | T7 +1 to +3; first transcribed bases |
| 50–68 | `AGATGTGTATAAGAGACAG` | **mosaic end (ME)**, 3' end = the transferred end |

- ME (50–68) pairs with ME-rc (1–19): a **19-bp double-stranded Tn5 binding site**, and
  nts 20–49 form a **30-nt single-stranded loop** carrying the promoter. This matches the
  paper's description exactly (19-bp ds binding site + ss T7 promoter loop, Fig. 1B). 🟡
  computed, 🟢 description.
- The promoter is written in the **sense** direction of the transposon's 5'→3', i.e. it
  points **into** the genomic insert from the ME side. 🟡
- There is no s5/s7 sequence at all: the Nextera s5/s7 arms are replaced by the loop.
  🟡
- The **second-strand primer = 8-nt UMI + `GGG` + ME** (22 nt after the UMI), i.e. the
  5' end of every transcript, so it anneals to the 3' end of every first-strand cDNA. 🟡
- The 5' phosphate on the transposon presumably lets the hairpin's 5' end be joined
  after gap fill; the paper text we hold does not say. 🔴

### The library-prep oligos — 🟡 (computed against `lib/illumina.py`)

- The hairpin adaptor equals `illumina.NEBNEXT_HAIRPIN` (65 nt, U at position 32):
  `NEBNEXT_ARM_READ2` + U + `TRUSEQ_READ1`.
- Universal PCR primer = `illumina.P5` + `TRUSEQ_READ1` overlapping by 4 nt (`ACAC`),
  58 nt (= `illumina.NEBNEXT_UNIVERSAL_PRIMER`).
- Index 1 primer = `illumina.P7` + i7 (6 nt) + `illumina.TRUSEQ_READ2`; the i7 is
  embedded as the reverse complement of what the index read reports (standard NEBNext). 🟡

## 3. Step by step

Conditions in the paper's main text are sparse; all buffer, time and temperature
details are in the supplement we could not fetch. 🔴

1. **Transposome assembly**: equimolar LIANTI transposon and Tn5 transposase, dimerised
   🟢 (Fig. 1B legend). Both monomers carry the identical hairpin adaptor (homodimer). 🟡
2. **Single cell picking and lysis** (BJ fibroblasts; "sort or handpick" per the
   upstream page 🟡). Optional **UDG** treatment of the lysate removes deaminated
   cytosines (C→U), the source of the C→T false positives 🟢 (Fig. 3B, main text).
3. **Tagmentation**: random fragmentation, average ~400 bp 🟢. Each fragment end carries
   the transposon's ME (3' end joined to the genomic 5' end), the hairpin and the usual
   9-nt gap on the opposite strand 🟡 ([Tn5](../ref/concepts/tn5-tagmentation.md)).
4. **Gap extension** with a DNA polymerase (Q5 per upstream 🟡) converts the
   single-stranded T7 loops into **double-stranded T7 promoters on both ends** 🟢 (Fig. 1C
   legend). Computed product (top strand, 🟡):
   transposon (68) · insert · reverse complement of the transposon (68).
5. **In vitro transcription, overnight** 🟢. Both ends have a promoter pointing inward, so
   both strands of the fragment are transcribed. Computed transcript (5'→3', 🟡):
   `GGG` · ME · <insert> · ME-rc · `CCC` · `TATAGTGAGTCGTATTA` (promoter rc) ·
   `AATTCTGTTC` · ME — 22 nt before the insert, 68 nt after.
   Upstream says the two convergent transcription events clash and one strand dominates
   (upstream claim, not in the main text) 🟡.
6. **Self-primed reverse transcription** 🟢 ("capable of self-priming on the 3' end").
   Computed mechanism 🟡: the transcript's 3'-terminal ME pairs with the ME-rc 30 nt
   upstream (the same 19-bp stem as the transposon, now in RNA, with a 30-nt loop), and
   the 3' OH primes RT back across the insert to the transcript's 5' end. First strand
   (after RNA removal): rc(insert) · `CTGTCTCTTATACACATCT` · `CCC`.
7. **RNase digestion** removes the RNA 🟢 (enzyme not named in the main text 🔴).
8. **Second-strand synthesis** with the UMI primer 🟢 (concept) / 🟡 (sequence) →
   double-stranded LIANTI amplicon: `<UMI8> · GGG · ME · <insert>`, with the far end bare
   genomic sequence (the other Tn5 insertion point). Yield ~20 ng per cell 🟢.
9. **Library prep** (NEBNext Ultra, per upstream 🟡): end repair, dA-tailing, ligation of
   the NEBNext hairpin adaptor, USER cleavage at the dU to open the hairpin, PCR with the
   Universal primer and an Index 1 primer. Cycle numbers 🔴.

## 4. Final library — 🟡 (assembled from the oligos above)

Both orientations occur, because the A-tailed amplicon is ligated on both ends.

Orientation 1 (UMI on the P5 side), 152 nt + insert (computed):

```
5'- P5 · TruSeq Read 1 (overlapping P5 by 4 nt) · <UMI 8> · GGG · ME · <genomic insert> · AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC · <i7 6, as rc> · P7' -3'
```

Orientation 2 (UMI on the P7 side):

```
5'- P5 · TruSeq Read 1 · <genomic insert> · ME-rc · CCC · <UMI 8 rc> · AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC · <i7 6, as rc> · P7' -3'
```

The upstream page draws the same two orientations. 🟡 agreement.

## 5. Read layout / sequencing

- Read 1 (TruSeq Read 1 primer): in orientation 1, **8 nt UMI + 22 nt `GGG`+ME, then
  genome at cycle 31**; in orientation 2, genome from cycle 1. 🟡 computed.
- Index 1 (sample index primer): 6-nt i7. 🟡
- Read 2 (TruSeq Read 2 primer): the mirror of Read 1. 🟡
- Both read starts are Tn5 insertion points, which is what makes "same start and end =
  same original fragment" digital counting possible (Fig. 2A) 🟢 / 🟡.
- Platform, read lengths: not in the main text 🔴. BJ cells sequenced to ~30× 🟢; 97 %
  genome coverage, 17 % ADO 🟢 (Table S1, as quoted in main text).

## 6. Upstream page vs the paper — agreements and disagreements

- **Agree**: transposon = 19-bp ds Tn5 site + ss T7 loop (computed from the upstream
  sequence, matches Fig. 1B legend); gap extension → ds promoters on both ends; overnight
  IVT; self-priming; RNase; second-strand synthesis with UMI-tagged primer.
- **Disagreement (internal to upstream)** 🟡: upstream's IVT product (its step 4) begins
  `GAGAUGUG…`, i.e. with only one G before ME, and its cDNA (step 6) correspondingly ends
  in `…CTC`. With the class III promoter in the transposon, transcription starts at the
  first G of `GGG`, so the RNA begins `GGGAGAUGUG…` — which is also what the upstream's own
  second-strand primer (`…GGGAGATG…`) assumes. Treat the one-G drawing as a simplification.
- **Not checkable**: every sequence, the Q5 polymerase, the NEB kit, the "one strand
  dominates" claim — all only on the upstream page until the supplement is read.

## 7. Open questions

- 🔴 The supplement's exact oligo sequences (including whether the UMI is 8 nt and how
  it is written), modifications, and the transposon annealing protocol.
- 🔴 Role of the 10-nt spacer `GAACAGAATT` and of the transposon's 5' phosphate (is the
  hairpin ligated to the genome after gap fill, or displaced?).
- 🔴 Which RT, which RNase, and whether IVT also yields RNA-templated RNA products (refs
  51–53 in the paper concern T7 RNA-dependent synthesis, suggesting the supplement
  discusses it).
- 🔴 Sequencer, read lengths, PCR cycles.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text and the upstream page but not the
supplement (.docx behind PMC's download gate; listed as `(manual)` in MANIFEST.tsv).
`tools/scrape_primers.py` found every oligo on the upstream page and labelled ME, ME-rc
and the Illumina/NEBNext parts correctly; it also picked up fragments of the upstream's
ASCII duplex drawings (e.g. `CGACTCACTATAGGG`, `CGAGAAGGCTAG`) as if they were oligos,
and named some by neighbouring words ("UMI", "dU", "i7"). It found nothing in the paper,
which contains no sequences. `scrape_primers.py --find` does not match the hairpin
adaptor across its inline `(dU)` (with either U or T at that position); it was checked by
hand against the upstream text instead. Transposon architecture and every derived structure were
computed with `lib/` (`nextera.ME`, `illumina.*`, `chemdraw.revcomp`).
