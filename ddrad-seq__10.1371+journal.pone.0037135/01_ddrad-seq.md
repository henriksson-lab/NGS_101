# ddRAD-seq

## 1. What it is

ddRAD-seq is a reduced-representation DNA protocol in which two restriction enzymes define fragment ends, compatible adapters introduce sample identities, and size selection chooses a reproducible subset of the genome.

## 2. Authoritative source

🟢 Peterson, Weber, Kay, Fisher and Hoekstra, *Double Digest RADseq: An Inexpensive Method for De Novo SNP Discovery and Genotyping in Model and Non-Model Species*, PLOS ONE (2012), DOI `10.1371/journal.pone.0037135`.

## 3. Construct-changing operations

🟢 A double restriction digest creates two distinguishable cohesive ends. Adapter ligation is end-selective; pooled fragments are size-selected and PCR completes the sequencing library.

## 4. Published oligos transcribed

🟢 Supplementary Protocol S1 and Table S1 give P1.1 as `ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATG`, P1.2 as `/5Phos/AATTCATGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`, P2.1 as `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, and P2.2 as `/5Phos/CGAGATCGGAAGAGCGAGAACAA`. GCATG is one of the published five-base inline barcodes.

🟢 The PCR2 example named ATCACG carries `CGTGAT` in the primer as the reverse complement; the Index 1 read reports ATCACG.

## 5. Construct interpretation

🟡 The final schematic expresses the published oligos as a conventional single-i7 TruSeq-compatible duplex. EcoRI and MspI overhangs are derived from their cut coordinates, and adapter compatibility is enforced when the model is constructed.
