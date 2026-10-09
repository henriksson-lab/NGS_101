# Human NET-seq

## 1. What it is

Human NET-seq purifies engaged RNA polymerase complexes and sequences the native 3′ ends of their nascent RNAs, so each aligned end reports polymerase position.

## 2. Authoritative sources

🟢 Mayer et al., *Native Elongating Transcript Sequencing Reveals Human Transcriptional Activity at Nucleotide Resolution*, Cell (2015), DOI `10.1016/j.cell.2015.03.010`.

🟢 Mayer and Churchman, detailed human NET-seq protocol, Nature Protocols (2016), DOI `10.1038/nprot.2016.086`.

## 3. Construct-changing operations

🟢 A pre-adenylated 3′ RNA linker carries a random molecular barcode. Reverse transcription, depletion of abundant unwanted cDNAs and circularization/PCR produce the sequencing library.

## 4. Oligos and topology

🟢 The detailed protocol prints the pre-adenylated DNA linker as `/5rApp/(N)6CTGTAGGCACCATCAAT/3ddC`, reverse PCR primer oNTI231 as `CAAGCAGAAGACGGCATACGA`, and custom sequencing primer oLSC006 as `TCCGACGATCATTGATGGTGCCTACAG`.

🟢 The linker is ligated to the nascent RNA 3′-OH before alkaline fragmentation. First-strand cDNA is size-selected, circularized with CircLigase, depleted with biotinylated probes, and amplified with an Illumina Index forward primer plus oNTI231.

🟡 The named commercial index-forward primer is not printed in the protocol, so the public construct keeps that outer arm as a role token. The custom read-primer site, N6 position and RNA-facing read direction are directly supported.
