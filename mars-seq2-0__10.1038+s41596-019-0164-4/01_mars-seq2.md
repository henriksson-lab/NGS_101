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
