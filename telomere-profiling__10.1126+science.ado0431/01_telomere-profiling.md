# Telomere Profiling

## 1. What it is

Telomere Profiling combines restriction digestion, a phased biotinylated TeloTag, streptavidin enrichment and Oxford Nanopore sequencing to concentrate intact human chromosome ends.

## 2. Sources read

- 🟢 Karimian et al. (2024), [doi:10.1126/science.ado0431](https://doi.org/10.1126/science.ado0431), author-manuscript article text and method description.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison.

## 3. Construct and reaction order

1. 🟢 Genomic DNA is restriction-digested while native telomere ends are retained.
2. 🟢 Six telomeric splints cover repeat phase and guide ligation of the barcoded, biotinylated TeloTag.
3. 🟢 Streptavidin capture and washing enrich tagged chromosome ends.
4. 🟢 Restriction cleavage at the adapter release site elutes selected fragments.
5. 🟢 Oxford Nanopore ligation library preparation adds the motor adapter.

The supplementary ordered bases were not available from the automated source endpoint. Named structural regions are retained without guessed sequence.

## 4. Final library and sequencing

The motor-loaded Oxford Nanopore adapter, rather than a sequencing primer, controls strand entry through the pore.
