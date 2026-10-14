# Telomere Profiling

## 1. What it is

Telomere Profiling combines restriction digestion, a phased biotinylated TeloTag, streptavidin enrichment and Oxford Nanopore sequencing to concentrate intact human chromosome ends.

## 2. Sources read

- 🟢 Karimian et al. (2024), [doi:10.1126/science.ado0431](https://doi.org/10.1126/science.ado0431), author-manuscript article text and method description.
- 🟢 Karimian et al. bioRxiv v2, [doi:10.1101/2023.12.21.572870](https://doi.org/10.1101/2023.12.21.572870), Extended Data Table 2 oligonucleotide sequences linked to the published article.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), workflow comparison.

## 3. Construct and reaction order

1. 🟢 Genomic DNA is restriction-digested while native telomere ends are retained.
2. 🟢 Six telomeric splints cover repeat phase and guide ligation of the barcoded, biotinylated TeloTag.
3. 🟢 Streptavidin capture and washing enrich tagged chromosome ends.
4. 🟢 Restriction cleavage at the adapter release site elutes selected fragments.
5. 🟢 Oxford Nanopore ligation library preparation adds the motor adapter.

🟢 Extended Data Table 2 publishes six complete 99-nt TeloTag strands with 5′ phosphate and 3′ BioTEG, plus six 105-nt phase splints for each tag. The molecular model stores all six TeloTags and derives each matching splint family by reverse-complementing the paired 87-nt tag body and appending the appropriate 18-nt telomere-repeat phase. NB50 and all six of its phase splints are shown as the representative duplex.

## 4. Final library and sequencing

The motor-loaded Oxford Nanopore adapter, rather than a sequencing primer, controls strand entry through the pore.

🔴 The defining workflow does not print the complete historical motor-bearing ONT adapter assembly. Current ONT adapter oligonucleotides are not substituted into this earlier library without evidence that the kit revision is identical.
