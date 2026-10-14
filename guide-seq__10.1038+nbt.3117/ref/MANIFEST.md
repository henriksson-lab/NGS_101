# Source manifest

Obtain everything with `python3 tools/get_sources.py --doi 10.1038/nbt.3117`; the files
land in `$CHEM_DATA/sources/guide-seq__10.1038+nbt.3117/` with `.txt` twins.

- `nbt.3117_PMC4320685.xml` / `.pdf` — Tsai et al. 2015, DOI `10.1038/nbt.3117`,
  PMC4320685 author manuscript: main text and Online Methods (the dsODN strands).
- `nbt.3117_41587_2015_BFnbt3117_MOESM21_ESM.pdf` — **the authoritative file**, Springer
  supplementary, corrected 25 June 2015: Supplementary Fig. 1 (library cartoon),
  Supplementary Table 4 (every common oligo) and the Supplementary Methods
  "Protocol for GUIDE-seq Library Preparation". The correction un-swaps the
  `Nuclease_off_+_GSP1` / `Nuclease_off_-_GSP1` labels, so an uncorrected copy of this
  file must not be used.
- `nbt.3117_41587_2015_BFnbt3117_MOESM22_ESM.xlsx` — off-target site tables; no oligos.

`python3 tools/selftest_source_oligos.py guide-seq__10.1038+nbt.3117` re-checks every
exact oligo in `tools/guide_seq.py` against the extracted text of these files, and skips
rather than fails when they are absent.
