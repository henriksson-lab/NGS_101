# Source manifest

Originals live in `$CHEM_DATA/sources/chip-seq__10.1038+nmeth1068/` (see its `MANIFEST.tsv`).

- Robertson et al. 2007, *Nature Methods* 4:651–657, DOI `10.1038/nmeth1068` (no dot;
  the dotted `nmeth.1068` does not resolve). Main article PDF is subscription-only: download
  it by hand from `https://doi.org/10.1038/nmeth1068` on an institutional network and save
  it as `nmeth1068.pdf` in the source directory (done 2026-10-17). Obtain the
  Springer supplements with `python3 tools/get_sources.py --doi 10.1038/nmeth1068`
  (Supplementary Information `MOESM305` PDF; Supplementary Data 1 `MOESM306` ZIP,
  whose README states the 27 bp read length).
- Schmidt et al. 2009, DOI `10.1016/j.ymeth.2009.03.001`, PMC4052679: author-manuscript
  JATS XML via `python3 tools/get_sources.py --doi 10.1016/j.ymeth.2009.03.001`.
- Illumina, *Preparing Samples for ChIP Sequencing of DNA*, part 11257047 Rev. A (2007):
  `https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/cluster_station/11257047_ChIP_Sample_Prep.pdf`
  (fetched by hand; no `extra_sources.tsv` row yet).
- Illumina oligonucleotide/adapter sequence document covering Adapter Oligo Mix 1000521
  and PCR primers 1000537/1000538: not obtained.
