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

🟢 Its Table 1 also prints the complete RT primer oLSC007 as `/5Phos/ATCTCGTATGCCGTCTTCTGCTTG/iSp18/CACTCA/iSp18/TCCGACGATCATTGATGGTGCCTACAG`. The two internal 18-atom spacers are chemical links, not nucleotide runs; the 5′ phosphate permits CircLigase closure. The table prints a separate RNA ligation control, oGAB11, `agucacuuagcgauguacacugacugug`. Table 2 gives all twenty 5′-biotin-TEG depletion oligos, which remove mature-RNA cDNAs rather than contributing to the sequenced construct.

🟢 Table 2 depletion oligos, each with a **5′ biotin-TEG** modification (all sequences 5′→3′):

| Target | Sequence after 5′ biotin-TEG |
|---|---|
| RNVU1-1 | `CAGGGGAAAGCGCGAACGCAGTCCCCCACTACC` |
| RNVU1-7 | `CAGGGGAAAGCGCGGACGCAGTCCCCCACT` |
| RNU2-1 | `TGCACCGTTCCTGGAGGTACTGCAATACCAGGTCGATG` |
| RNU4-1 | `GCCTACAGCAGTCTCCGTAGAGACTGTCAAAAATTGCCA` |
| RNU5B-1 | `AGCCTTGTCGGAACAAGGCCTCAAAAAATTAGCTTAAGAC` |
| RNU5E-1 | `AGCCTTGCCAAAGCAAGGCCTCAAAAAATTGGGTTAAG` |
| SNORD3D | `CCACTCAGACCGCGTTCTCTCCCTCTCACTCCCCAA` |
| SNORD12B | `CAGACAGAAACTGGCTTAGAAAAGTTAGATCAACATAGTCGATC` |
| SNORD12C | `CTAACTGGCAAAATATAAGACGTCAGCATTGTCGATCTGATG` |
| SNORD27 | `CTCAGTAGTAAGATGACATCACTTGAAAGTTCAGCCATATGC` |
| SNORD29 | `CTCAGGTGTTCATGTATTTTCACTGTCGGTCATAGTGAGC` |
| SNORD31 | `TACCTTTCAGTCACACATTGATCAGACTGGGGCGGTAT` |
| SNORD80 | `CATCAGATAGGAGCGAAAGACTTAATATTGCTCATCAGCG` |
| SNORD81 | `TTGTCATCAAGTAATCAGTGAGAGAGTTCAAGTTGG` |
| SNORD118 | `GGAGCAATCAGGGTGTTGCAAGACCTGATTACGCAGAG` |
| RNA5S1 | `CACCCGGTATTCCCAGGCGGTCTCCCATCCAAGTAC` |
| RNA5-8S5 | `AGCGACGCTCAGACAGGCGTAGCCCCGGGAGG` |
| RNA28S5 | `ACAAACCCTTGTGTCGAGGGCTGACTTTCAATAGATCGCA` |
| MT-TM | `GTACGGGAAGGGTATAACCAACATTTTCGGGGTATGGG` |
| MT-TV | `TGGTCAGAGCGGTCAAGTTAAGTTGAAATCTCCTAAGTG` |

🟢 The linker is ligated to the nascent RNA 3′-OH before alkaline fragmentation. First-strand cDNA is size-selected, circularized with CircLigase, depleted with biotinylated probes, and amplified with an Illumina Index forward primer plus oNTI231.

🟡 The named commercial index-forward primer is not printed in the protocol, so the public construct keeps that outer arm as a role token. The custom read-primer site, N6 position and RNA-facing read direction are directly supported.
