# CUT&Tag

## 1. What it is

CUT&Tag maps chromatin-bound proteins by recruiting adapter-loaded Protein A–Tn5 through a specific antibody, so target-proximal cleavage and sequencing-adapter transfer are the same molecular event.

## 2. Sources and evidence

- 🟢 Kaya-Okur HS et al. *CUT&Tag for efficient epigenomic profiling of small samples and single cells.* Nature Communications 10:1930 (2019). [doi:10.1038/s41467-019-09982-5](https://doi.org/10.1038/s41467-019-09982-5).
- 🟢 The Methods specify an equimolar preannealed Tn5MEDS-A/Tn5MEDS-B mixture loaded onto pA–Tn5, antibody binding in permeabilized cells, magnesium activation, gap-filling extension and indexed PCR.
- 🟡 The completed page uses canonical Nextera sequences from `lib/nextera.py`; the paper names the standard Tn5 adapters rather than reprinting every base in its main text.

## 3. Molecular path

Antibody binds the chromatin epitope, Protein A recruits a loaded Tn5 dimer, and washing removes unbound transposase. Magnesium activates target-proximal tagmentation. Gap fill and indexed PCR select heterologous S5/S7 products and complete the Illumina library.
