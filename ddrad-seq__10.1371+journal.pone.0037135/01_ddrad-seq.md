# ddRAD-seq

## 1. What it is

ddRAD-seq is a reduced-representation DNA protocol in which two restriction enzymes define fragment ends, compatible adapters introduce sample identities, and size selection chooses a reproducible subset of the genome.

## 2. Authoritative source

🟢 Peterson, Weber, Kay, Fisher and Hoekstra, *Double Digest RADseq: An Inexpensive Method for De Novo SNP Discovery and Genotyping in Model and Non-Model Species*, PLOS ONE (2012), DOI `10.1371/journal.pone.0037135`.

## 3. Construct-changing operations

🟢 A double restriction digest creates two distinguishable cohesive ends. Adapter ligation is end-selective; pooled fragments are size-selected and PCR completes the sequencing library.

## 4. Published oligos transcribed

🟢 Supplementary Protocol S1 and Sequences S1 give P1.1 as `ACACTCTTTCCCTACACGACGCTCTTCCGATCTGCATG`, P1.2 as `/5Phos/AATTCATGCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`, P2.1 as `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`, and P2.2 as `/5Phos/CGAGATCGGAAGAGCGAGAACAA`. GCATG is one of the published five-base inline barcodes.

🟢 The PCR2 example named ATCACG carries `CGTGAT` in the primer as the reverse complement; the Index 1 read reports ATCACG.
🟢 The same supplementary spreadsheet prints PCR1 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACG` and the complete PCR2_Idx_1_ATCACG oligo `CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC`. Its `PCR_primers` sheet lists all twelve PCR2 index alternatives; its `adapter_P1-EcoRI` sheet lists all 48 five-base inline barcode pairs. The representative schematic uses one published member of each family.

## 5. Construct interpretation

🟡 The final schematic expresses the published oligos as a conventional single-i7 TruSeq-compatible duplex. EcoRI and MspI overhangs are derived from their cut coordinates, and adapter compatibility is enforced when the model is constructed.

## 6. Bounded source pass — 2026-10-17

Read the published PMC3365034 JATS Methods (lines 57–67 and combinatorial
indexing at 113–122), Protocol S1 `pone.0037135.s001.doc.txt` lines 8–12,
18–47 and 53–115, and Sequences S1 `pone.0037135.s002.xls.txt` representative
EcoRI/MspI rows 3, 51, 101–102 and PCR sheet 205–219. Manifest identity and
article DOI agree; the oligo scraper was run over the bundle and relevant
modified hits checked in these passages.

🟢 Protocol S1 specifies T4 DNA ligase and Phusion PCR. Both P1.2 and P2.2
carry the printed `/5Phos/` modification; the schematic now displays these
ends and the four adapter ordering strands.

🟢 Protocol S1 lines 90–93 explicitly says P1/P1 products can amplify but do
not cluster. Its biotin-P2/streptavidin removal is optional and occurs after
size selection. The graph follows the route without this optional capture.
The previous immediate discard of same-end digest fragments was a repository
model error: no separation happens at digestion. Same-end ligation products
now have an unresolved downstream fate rather than an invented physical
removal. Their molecular states and subsequent recovery remain unrepresented.

🟢 Article lines 117–122 permit single- or paired-end runs with or without
an index read. The schematic's exact primer sites establish compatibility;
this generic example does not assign one study's acquired cycle recipe.

Deferred: Analysis S1 (`s003.pdf`, not read in this capped pass); the linked
live author protocol, design/calculator documents and vendor procedures;
explicit post-digest, pooled and post-PCR AMPure recovery states; complete
adapter annealing and pre-PCR strand geometry; same-end product propagation
and a selected experimental run. Historical source readings remain historical.
This is a bounded pass, not a new independent absence review or complete
experimental validation.
