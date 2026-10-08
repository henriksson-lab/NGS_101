# CROP-seq — third-party reference material

No third-party source files are committed here. Re-obtain the sources used for the note and
model from Datlinger *et al.*, doi:[10.1038/nmeth.4177](https://doi.org/10.1038/nmeth.4177):

- PMC full text, PMCID **PMC5334791**;
- Supplementary Protocol (cloning, Drop-seq, oligo table and sequencing configuration);
- Supplementary Table 1 (vector-construction and validation oligos); and
- Supplementary Data, the **CROPseq-Guide-Puro / Addgene #86708** GenBank map.

`python3 tools/get_sources.py "CROP-seq"` stores fetched copies under `$CHEM_DATA/sources/`.
The repository model contains only small sequence facts transcribed from those sources.
