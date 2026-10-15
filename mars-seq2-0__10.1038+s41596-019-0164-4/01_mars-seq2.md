# MARS-seq2.0 — source note

The full source audit, oligo transcription and step-by-step reading are retained in the
companion [MARS-seq / MARS-seq2.0 research note](../mars-seq__10.1126+science.1247651/01_mars-seq.md),
because one model (`../mars-seq__10.1126+science.1247651/tools/`) renders both versions.

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
constant or a transcription from one of two **authors' own documents**, both on disk:

| Modelled material | Bases | Authoritative document |
|---|---|---|
| RT1: `CGATTGAGGCCGG` · T7 promoter · `GC` · `GACGTGTGCTCTTCCGATCT` · 7-nt well barcode · N8 · T20 · N | 🟢 verbatim | **MARS-seq2.0 Supplementary Table 1**, sheet `MARS-seq2.0_RT1 primers` — all 384 rows match this form |
| pool ligation adapter: `/5Phos/` · 4-nt pool barcode · N5 · `AGATCGGAAGAGCGTCGTGTAG` · `/3SpC3/` | 🟢 verbatim | **MARS-seq2.0 Supplementary Table 2**, sheet `MARS-seq2.0_Ligation adapter` — all 32 rows |
| read design `3I.4P.50M` / `7W.8R.5I` | 🟢 verbatim | MARS-seq2.0 Supplementary Information (`MOESM1`), pipeline `seq_batches.txt` field description |
| MARS-seq (v1) RT1, ligation adapters `lig_NNNX4_ix1..8`, second RT primer `TCTAGCCTTCTCGCAGCACATC`, `P5_Rd1 PCR`, `P7_Rd2 PCR` | 🟢 verbatim | **Jaitin *et al.* 2014, *Science* 343:776, Supplementary Material**, Tables S7–S8 — the predecessor this protocol versions |
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

## What remains unpublished

🔴 The Nature Protocols **main text is paywalled**, so for MARS-seq2.0 itself the second
RT primer (RT2) and the two library-PCR primer sequences were never read. The model does
not guess them: the MARS-seq2.0 final library is drawn with `[P5-SIDE ARM]` and
`[P7-SIDE ARM]` inferred placeholders, `rt2`, `pcr_forward` and `pcr_reverse` are
`inferred=True` graph inputs, the library PCR carries `PCRGeometry.unavailable(...)`, the
sequencing section renders `seqprimers.unavailable_diagram(...)`, and the page states the
boundary in its caveat.

🟡 Whether v2 reuses the v1 primers printed in Jaitin Table S7 is **not** established. The
v1 page may use them (they are published for v1); the v2 page may not, and does not.

🔴 All v2 reaction conditions (lysis buffer, RT enzyme and volumes, exonuclease,
second-strand mix, IVT time, fragmentation, ligation, RT2, PCR cycles) are in the
paywalled main text.
