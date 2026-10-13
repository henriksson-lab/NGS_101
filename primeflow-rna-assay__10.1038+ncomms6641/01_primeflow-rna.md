# PrimeFlow RNA Assay

## 1. What it is

A commercial branched-DNA fluorescence in situ hybridization assay for measuring selected RNA targets in individual suspension cells by flow cytometry, optionally together with surface and intracellular proteins.

## 2. Sources and evidence

- 🟢 [Thermo Fisher PrimeFlow RNA Assay Kit User Guide, MAN0019788 Rev B.0](https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0019788_PrimeFlowRNAAssay_UG.pdf), dated 22 February 2022.
- 🟢 Porichis F, Hart MG, Griesbeck M, et al. *High-throughput detection of miRNAs and gene-specific mRNA at the single-cell level by flow cytometry.* Nature Communications 5:5641 (2014). [doi:10.1038/ncomms6641](https://doi.org/10.1038/ncomms6641).
- 🟢 [Supplementary Table 1](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms6641/MediaObjects/41467_2014_BFncomms6641_MOESM1192_ESM.pdf) gives the target binding intervals and target-region sequences used for the IL-2, IL-21, IFN-γ, IDO, HMBS and miRNA probe sets in the paper.
- 🟢 The manual specifies 20–40 adjacent target-probe pairs, then preamplifier, amplifier and fluorescent label-probe hybridizations. One complete tree has 400 label sites and the full target set gives 8,000–16,000-fold theoretical signal amplification.
- 🟢 The current assay supports up to four targets using Types 1, 4, 6 and 10, followed by standard flow-cytometric analysis.
- 🟡 Those published target regions determine the reverse-complementary target-binding arms. They are not complete orderable probe sequences: the universal tails attached to the arms are absent from the table.
- 🔴 The universal probe tails (including their lengths and docking bases) and the
  preamplifier, amplifier and fluorescent label-probe sequences remain undisclosed.

## 3. Molecular path

After optional antibody staining, suspension cells are fixed and permeabilized. Gene-specific probe pairs bind adjacent regions of each target RNA. Preamplifier and amplifier hybridizations build a branched scaffold, fluorescent label probes occupy its terminal sites, and a flow cytometer measures fluorescence in each retained cell. The target RNA is never amplified.

## 4. Public sequence example

For human IL-2 (NM_021803), the first adjacent target regions are bases 2–24, `CCTTGGTCTCGTTTTCACTTCAG`, and bases 25–52, `TCTCATAAGTACCAACAGTAGAGCTAGA`. The schematic derives the two complementary probe arms from these published target sequences and leaves each proprietary universal tail visibly unresolved.
