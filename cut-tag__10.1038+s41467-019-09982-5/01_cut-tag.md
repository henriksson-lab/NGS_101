# CUT&Tag

## 1. What it is

CUT&Tag maps chromatin-bound proteins by recruiting adapter-loaded Protein A–Tn5 through a specific antibody, so target-proximal cleavage and sequencing-adapter transfer are the same molecular event.

## 2. Sources and evidence

- 🟢 Kaya-Okur HS et al. *CUT&Tag for efficient epigenomic profiling of small samples and single cells.* Nature Communications 10:1930 (2019). [doi:10.1038/s41467-019-09982-5](https://doi.org/10.1038/s41467-019-09982-5).
- 🟢 The Methods specify an equimolar preannealed Tn5MEDS-A/Tn5MEDS-B mixture loaded onto pA–Tn5, antibody binding in permeabilized cells, magnesium activation, gap-filling extension and indexed PCR.
- 🟢 The CUT&Tag Methods cite Picelli et al., [doi:10.1101/gr.177881.114](https://doi.org/10.1101/gr.177881.114), for the preannealed Tn5MEDS-A/B loaded onto pA–Tn5. Picelli's Methods print Tn5ME-A `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, Tn5ME-B `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`, and the shared reverse mosaic-end strand `5′-/5Phos/CTGTCTCTTATACACATCT-3′`. The A/B bodies match the canonical segments in `lib/nextera.py`.
- 🟡 The defining paper names the standard adapter pair rather than reprinting these bases. Its own Supplementary Files 1–3 do not supply a different adapter sequence; the sequence assignment follows its cited construction source.
- 🔴 The defining methods do not specify which particular Nextera N5/N7 sample-index members were used for each library; the named pair establishes adapter geometry, not a sample barcode choice.

## 3. Molecular path

Antibody binds the chromatin epitope, Protein A recruits a loaded Tn5 dimer, and washing removes unbound transposase. Magnesium activates target-proximal tagmentation. Gap fill and indexed PCR select heterologous S5/S7 products and complete the Illumina library.
