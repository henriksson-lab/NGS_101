# scCARE-seq — source notes

Defining paper: Deng et al., *Nature Structural & Molecular Biology* (2023), [doi:10.1038/s41594-023-01066-9](https://doi.org/10.1038/s41594-023-01066-9).

- 🟢 The paper, Extended Data Figure 1 and Supplementary Table 1 were read. The
  source files are the publisher's article XML, figure and supplementary workbook.
- 🟢 Supplementary Table 1 prints Tn5-N5, Tn5-N7, phosphorylated Tn5-Rev,
  carrier ssDNA, the common/DNA/RNA pre-amplification primers, all eight
  UMI–barcode–oligo(dT) RT primers, and the complete N501–N516, N701–N724 and
  T501–T516 indexed PCR-primer sets.
- 🟢 A biotin-dT capture-primer form labels RNA; the method combines RNA expression with single-cell chromatin-contact reads.
- 🟢 Extended Data Figure 1 assigns N5xx/N7xx to the DNA-contact branch and
  T5xx/N7xx to the RNA branch. Thus the RNA endpoint is a hybrid library: TruSeq Read 1
  and i5 geometry at one end, Nextera Read 2 and i7 geometry at the other.
- 🟡 The page shows Oligo-dT-1, N501, N701 and T501 as exact representatives;
  the remaining printed members change only the documented whitelist sequence.
