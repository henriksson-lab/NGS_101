# GAGE-seq — source notes

Defining preprint: Zhou et al. (2023), [doi:10.1101/2023.07.20.549578](https://doi.org/10.1101/2023.07.20.549578). Published version: Zhou et al. (2024), [doi:10.1038/s41588-024-01745-3](https://doi.org/10.1038/s41588-024-01745-3).

- 🟢 The full preprint, published paper, supplementary protocol and the published Supplementary Table 1 workbook were read.
- 🟢 RNA is reverse-transcribed with biotinylated poly(dT) or random-hexamer primers; genomic DNA undergoes in-situ contact ligation.
- 🟢 DNA and cDNA receive two rounds of 96-well ligation barcodes using analyte-specific adapter sets. Streptavidin beads then separate cDNA from DNA.
- 🟢 Both branches are preamplified, Nextera-tagmented, indexed and sequenced PE150.
- 🟢 The version-of-record workbook supplies the material that was missing from
  the preprint bundle: every 96-well RNA and Hi-C adapter plate, strand-specific
  modifications, RNA splints and blockers, the post-TdT primers, Hi-C pre-amplification
  primers, and terminal index-primer sets.
- 🟢 The two RNA ligation rounds can be reconstructed exactly. Linker 1 bridges
  the RT-primer handle to Round1_01; the RNA linker 2 bridges Round1_01 to Round2_01.
  Both complementary arms are checked by the shared splinted-ligation constructor.
- 🟢 The terminal shells differ by branch, as Supplementary Figure 1 indicates:
  chromatin-contact molecules use Nextera P5 with TruSeq P7, whereas RNA molecules use
  TruSeq P5 with Nextera P7. Representative A1/index primers are modeled exactly;
  alternate wells change only their printed whitelist sequences.
