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

## 3. Molecular workflow

🟢 Six million cells are extracted to uHMW DNA. The protocol mixes 750 µl DNA with
250 µl Fragmentation Mix diluted 6:244 in FRA Dilution Buffer. It incubates 10 min at
room temperature, 10 min at 75 °C and at least 10 min on ice. Five microlitres of Rapid
Adapter are added and incubated 30 min at room temperature.

🟢 The library is recovered with Precipitation Buffer and a Precipitation Star and eluted
overnight. Wide-bore tips and slow mixing preserve molecule length. This cleanup is part
of why the ultra-long workflow differs from ordinary rapid preparation.

🔴 Fragmentation tag, Rapid Adapter, attachment chemistry, motor-bearing strand and exact
adapter bases are proprietary. The diagram therefore uses role tokens, not nucleotide
sequences borrowed from other Nanopore kits.

## 4. Sequencing entry

🟢 No sequencing primer binds. A motor in the Rapid Adapter meters the attached strand
through an R10.4.1 nanopore. The model marks the tag–adapter boundary as an attachment,
not as a ligation, because no ligase is used in library construction.
