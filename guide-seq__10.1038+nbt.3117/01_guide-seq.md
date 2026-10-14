# GUIDE-seq — research notes

## Assay target

The sequenced material is human genomic DNA flanking nuclease-induced double-strand breaks in living cells. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (a non-coding 34-bp tag duplex, a single-tail adapter and nested Illumina primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Tsai et al. 2015, DOI [10.1038/nbt.3117](https://doi.org/10.1038/nbt.3117), PMC4320685, main text, Online Methods and supplement.

## Molecular order

🟢 A blunt, 5-prime-phosphorylated 34-bp dsODN with phosphorothioate-protected ends is
introduced with the programmable nuclease and captured by NHEJ at cellular DSBs. Genomic
DNA is randomly sheared and ligated to a single-tail adapter. Two mirrored STAT-PCRs,
each using a dsODN-strand primer and adapter primer, recover the two genomic flanks.
Nested amplification installs an 8-bp random molecular barcode and Illumina arms.

## Boundary

GUIDE-seq measures breaks that capture the dsODN in living cells. It is distinct from
CIRCLE-seq, where purified genomic circles are cut by Cas9 in vitro.
