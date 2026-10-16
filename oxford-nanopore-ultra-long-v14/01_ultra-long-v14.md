# Oxford Nanopore Ultra-Long DNA Sequencing Kit V14

## 1. What it is

SQK-ULK114 begins with ultra-high-molecular-weight genomic DNA, deliberately uses a very
dilute fragmentation mix to create few cuts, attaches Rapid Adapter, and recovers the
extremely long adapted DNA by precipitation rather than bead manipulation. The protocol
targets read-length N50 above 50 kb and is PCR-free.

Commercial defining source: Oxford Nanopore, “Ultra-Long DNA Sequencing Kit V14
(SQK-ULK114),” ULK_9177_v114.

## 2. Source read

- 🟢 [Oxford Nanopore protocol ULK_9177_v114](https://nanoporetech.com/document/ultra-long-dna-sequencing-kit-sqk-ulk114), including extraction, tagmentation, Rapid Adapter attachment, precipitation cleanup and loading.
- 🟢 Oxford Nanopore, *Chemistry Technical Document*, CHTD_500 v1 rev AV, ultra-long kit description and reagent table (text lines 421–423 and 2124–2125), plus adapter table (lines 2081–2087).

## 3. Molecular workflow

🟢 Six million cells are extracted to uHMW DNA. The protocol mixes 750 µl DNA with
250 µl Fragmentation Mix diluted 6:244 in FRA Dilution Buffer. It incubates 10 min at
room temperature, 10 min at 75 °C and at least 10 min on ice. Five microlitres of Rapid
Adapter are added and incubated 30 min at room temperature.

🟢 The library is recovered with Precipitation Buffer and a Precipitation Star and eluted
overnight. Wide-bore tips and slow mixing preserve molecule length. This cleanup is part
of why the ultra-long workflow differs from ordinary rapid preparation.

🟢 The published Rapid Adapter strand is `CCTGTACTTCGTTCAGTTACGTATTGCT`.

🟢 The technical document says the transposase fragments and tags DNA, then the
motor-loaded Rapid Adapter attaches to those tags without an enzyme. 🔴 The tag bases,
exact tag–adapter junction and complete motor-bearing assembly remain unresolved from
the reviewed vendor documents. Those regions remain role tokens.

## 4. Sequencing entry

🟢 No sequencing primer binds. A motor in the Rapid Adapter meters the attached strand
through an R10.4.1 nanopore. The model marks the tag–adapter boundary as an attachment,
matching the vendor's enzyme-free description.
