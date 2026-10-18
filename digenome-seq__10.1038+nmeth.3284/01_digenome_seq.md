# Digenome-seq — research notes

Kim et al., *Nature Methods* (2015), [10.1038/nmeth.3284](https://doi.org/10.1038/nmeth.3284).

## Source scope and version boundary

This capped migration read two originals: the defining2015 published supplement
`41592_2015_BFnmeth3284_MOESM426_ESM.pdf` (identity, contents, supplementary figure
legends and Supplementary Note1 instrument scope), and the authors’2016 multiplex
Digenome-seq companion article, [10.1101/gr.199588.115](https://doi.org/10.1101/gr.199588.115),
PMC4772022 Discussion114 and Methods132–136. The defining2015 main Methods were
not inspected. No complete source audit or publication-absence conclusion follows.
Mandatory bundle-scraper hits were reviewed; targeted-validation primers in the
supplement do not establish whole-genome-library adapters or index geometry.

🟢 The2016 source describes RNaseA50µg/ml followed by DNeasy repurification,
Covaris400–500bp fragmentation,1µg DNA ligated with adapters using a TruSeq DNA
library prep kit, and HiSeqXTen WGS at30–40×depth. Its Discussion explicitly says
no PCR amplification precedes WGS. These are qualified as2016 companion-method
statements, not proof of the original2015 catalogue member or run recipe.
Supplementary Note1’s HiSeqX10 statement concerns knockout-clone WGS and is not
silently reassigned to the cleaved-genome samples.

## Source-owned model

The local native graph represents cell-free Cas9 cleavage, RNase/repurification,
unenriched fragmentation and a source-limited WGS finishing boundary. Independent
Cas9/shear and random-shear fragment classes occupy separate Scene rows.
`batch_ngs.truseq_library`, `protocol_model.da_tailed` and the unused legacy
workflow are removed.

The former explicit end-repair, size-selection, dA-tail and cohesive ligation
choices were generic kit assumptions. They now remain unresolved inside the named
library-finishing boundary, rather than assigning uninspected source chemistry.
The final construct has variable insert and unknown physical adapter/index ends;
no six-base index, single/dual-index format, dA junction, canonical arm or acquired
read mask is added. The2016 no-PCR statement remains explicit without choosing a
PCR-Free catalogue product. This corrects unsupported molecular exactness rather
than establishing that the sources omit these details.

## Deferred evidence

The original2015 main Methods, exact historical kit member and physical adapter
ends, detailed library processing, index architecture/member/width and acquired
run cycles remain deferred. Previously consulted2017 PCR-Free and2025 adapter
vendor references are historical candidate family sources; they are not used to
supply study-specific defaults in this migration. The2021 detailed protocol and
its applicability/version identity were not investigated. No new retrieval or
reference chase was performed.
