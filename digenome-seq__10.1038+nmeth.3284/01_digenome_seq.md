# Digenome-seq — research notes

Kim et al., *Nature Methods* (2015), [doi:10.1038/nmeth.3284](https://doi.org/10.1038/nmeth.3284).

## Assay target

The sequenced material is purified human genomic DNA cleaved in vitro by Cas9, sequenced as an ordinary whole genome. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a stock TruSeq whole-genome library kit, with no cut-end-specific oligo). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining article abstract and integrated supplementary text/figures, including the 5′-end plots and Supplementary Notes.
- 🟢 Kim et al., *Genome Research* (2016), [doi:10.1101/gr.199588.115](https://doi.org/10.1101/gr.199588.115), the same group's open-access multiplex Digenome-seq paper, for explicit library details: Covaris 400–500-bp fragmentation, 1 µg input to TruSeq DNA library preparation, HiSeq X Ten WGS, and 30–40× depth.
- 🟢 Kim, Kang and Kim, *Nature Protocols* (2021), [doi:10.1038/s41596-020-00453-6](https://doi.org/10.1038/s41596-020-00453-6), for the protocol-level statement that Digenome-seq uses about 30× WGS without cut-end enrichment by PCR.
- 🟢 Illumina TruSeq DNA PCR-Free documentation and the repository's canonical Illumina adapter-sequence source for the commercial WGS shell.

## Construct logic

1. 🟢 Purified genomic DNA is cleaved in vitro with Cas9–sgRNA. Digenome-seq therefore measures biochemical cleavage of naked DNA, not editing inside a cell.
2. 🟢 RNase A removes guide RNA, protein is removed, and the DNA is repurified. No barcode or affinity handle is attached to a cut end.
3. 🟢 The entire digested genome is sheared by Covaris. The resulting pool contains a small class of fragments bounded by a Cas9 cut and a random shear, among the much larger population bounded only by random shearing.
4. 🟢 The fragments enter an ordinary Illumina whole-genome library preparation. The defining implementation used HiSeq X Ten WGS; the same authors' detailed follow-up specifies TruSeq DNA library prep and 30–40× depth, and the protocol description explicitly distinguishes Digenome-seq from PCR-enriched end-capture assays.
5. 🟢 After alignment, a cleavage site produces forward- and reverse-strand reads sharing precise 5′ endpoints. Random shearing distributes endpoints, so it does not make the same local pileup.

## What is—and is not—encoded in the molecule

There is no Digenome-specific oligo, UMI, cut-site barcode, or selective adapter. The standard sample index identifies a library, but it does not identify a cut end. The assay signal is the genomic coordinate of the physical insert boundary, interpreted across many ordinary WGS molecules.

The page consequently shows the standard single-index TruSeq PCR-free construct and declares the genomic insert boundary as the measurement. It does not add a fictional proximal adapter or enrichment step.
