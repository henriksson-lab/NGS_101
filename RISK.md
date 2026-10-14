# Risk review of the CRISPR / break-mapping protocol group

An automatic approval review once rejected primer-level work on the 13 protocols listed
below and left their oligo source audit pending. That audit is now complete: every one
of the 13 has a reviewed row in `catalogue/protocol_oligo_audit.tsv`, which remains the
authoritative status for every protocol.

## Why the group was flagged, and what it actually contains

Twelve of the 13 are **host-genome** assays. They map double-strand breaks, repair
outcomes, integration sites or integrated screen cassettes in human DNA; several exist
specifically to *measure* off-target editing, which is safety-side work. They amplify no
pathogen, and their oligos are tag duplexes, adapters, barcodes and locus primers.

The thirteenth, ARTIC tiled-amplicon sequencing, does sequence a viral genome from
clinical samples. It is a diagnostic and genomic-surveillance protocol whose primer
schemes the authors published in full and the ARTIC Network maintains openly for that
purpose; what this repository documents is the library chemistry downstream of the two
primer pools.

So the flag tracked vocabulary — CRISPR, off-target, double-strand break, amplicon — not
the chemistry. Each of the 13 now carries an `## Assay target` section in its notes
stating what is actually sequenced, so the same inference does not have to be re-derived.
Grep for that heading.

| Protocol | Directory | Audit |
|---|---|---|
| CRISPR-MIP | `crispr-mip__10.1101+2024.03.28.587082` | reviewed |
| CRISPR-UMI (Schmierer) | `crispr-umi-schmierer__10.15252+msb.20177834` | reviewed |
| ARTIC tiled-amplicon sequencing | `artic-tiled-amplicon-sequencing__10.1038+nprot.2017.066` | reviewed |
| GUIDE-seq | `guide-seq__10.1038+nbt.3117` | reviewed |
| CRISPR CIRCLE-seq | `crispr-circle-seq__10.1038+nmeth.4278` | reviewed |
| CRISPR-sciATAC | `crispr-sciatac__10.1038+s41587-021-00902-x` | reviewed |
| DSBCapture | `dsbcapture__10.1038+nmeth.3960` | limited |
| UDiTaS | `uditas__10.1186+s12864-018-4561-9` | reviewed |
| SITE-Seq | `site-seq__10.1038+nmeth.4284` | reviewed |
| Digenome-seq | `digenome-seq__10.1038+nmeth.3284` | limited |
| PEM-seq | `pem-seq__10.1038+s41421-019-0088-8` | reviewed |
| CHANGE-seq | `change-seq__10.1038+s41587-020-0555-7` | limited |
| DISCOVER-Seq | `discover-seq__10.1126+science.aav9023` | reviewed |

## What the audit found

Missing primers were not the problem. Nine of the 13 publish every structural oligo
verbatim; the three `limited` rows are limited by commercial kit oligos or an
inaccessible article Methods, not by withheld sequences. Two protocols (PEM-seq, UDiTaS)
have a *per-experiment by design* locus primer with published exemplars, which the
ledger records as such rather than as unpublished.

What the audit did find is **model fidelity**: pages that disagree with sources our own
`_data/sources/` already holds. Fixed during the audit:

- SITE-Seq `SITE_A2_R` was transcribed one base short (`ACTACTAATACGAC` for the published
  `ACTACTAATACGAC*T`), with the phosphorothioate rendered on the wrong base. The dropped
  T is the 3' overhang that T/A-ligates to the dA-tailed end, so Adapter 2 was drawn as a
  blunt duplex with no ligation chemistry. The selftest pinned the truncated string.
- CRISPR-UMI: two hand-typed insert displays ended in 7 T where the source has 6 — the
  error the same page announces as already fixed.
- CHANGE-seq: the circularization statement named CIRCLE-seq inside the CHANGE-seq graph.
- DSBCapture: the flow-cell claim cited main Figure 1d; it is Supplementary Figure 1d.
- ARTIC: the paper's own primer-concentration and denaturation figures disagree between
  its prose and its reaction table; the notes now record both.
- CRISPR-UMI: Michlits (`10.1038/nmeth.4466`), the defining source for note 03, had no
  `catalogue/extra_sources.tsv` row and was reachable only from a gitignored local PDF.

Outstanding, and **not** risk items:

- **CRISPR-sciATAC** needs a real rewrite. Supplementary Table 2 publishes every oligo as
  complete 96/144-member plates and the model transcribes none of it: a hardcoded 12-nt
  barcode that the source contradicts (all three are 8 nt) and that the selftest pins;
  two cell barcodes modelled where there are three; BC2/BC3 mislabelled as sample indexes
  when the protocol has no sample index; wrong barcode placement; both Tn5 arm spacers
  absent; canonical Nextera/TruSeq sequencing primers where the source has three custom
  ones; no RT step although guide identity comes from a reverse-transcribed transcript;
  and TnY described as Tn5. `catalogue/properties.tsv` repeats several of these.
- **GUIDE-seq** and **CIRCLE-seq** declare adapter geometry `unavailable` that their
  sources publish in full. GUIDE-seq additionally draws a pure TruSeq library where the
  real one is a hybrid, and places the `NNWNNWNN` molecular index inside Read 1 instead of
  the index-2 read.
- **PEM-seq**'s inline index is 4 nt, not the 5 the model and notes assert.
- **Digenome-seq** marks steps green against follow-up papers and Illumina documentation
  that have no manifest row, and renders an inherited TruSeq shell as established.
- **DSBCapture**'s oligo evidence is figure artwork only; the Protocol Exchange entry its
  Methods cites (printed as a literal "(Insert URL)") is not in `extra_sources.tsv`.
- **CHANGE-seq**'s defining article Methods was never fetched; both PMC routes are gated.
- `catalogue/properties.tsv` claims `exact oligos public` for Digenome-seq, which
  publishes no structural oligo at all.
