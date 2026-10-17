# Oxford Nanopore Ultra-Long DNA Sequencing Kit V14

## 1. What it is

SQK-ULK114 begins with ultra-high-molecular-weight genomic DNA, deliberately uses a very
dilute fragmentation mix to create few cuts, attaches Rapid Adapter, and recovers the
extremely long adapted DNA by precipitation rather than bead manipulation. The protocol
targets read-length N50 above 50 kb and is PCR-free.

Commercial defining source: Oxford Nanopore, “Ultra-Long DNA Sequencing Kit V14
(SQK-ULK114),” ULK_9177_v114 Rev Q, 6 November 2025.

## 2. Source read

- 🟢 [Oxford Nanopore protocol ULK_9177_v114](https://nanoporetech.com/document/ultra-long-dna-sequencing-kit-sqk-ulk114), including extraction, tagmentation, Rapid Adapter attachment, precipitation cleanup and loading.
- 🟢 Oxford Nanopore, *Chemistry Technical Document*, CHTD_500 v1 rev AV, ultra-long kit description and reagent table (text lines 421–423 and 2124–2125), plus adapter table (lines 2081–2087).

Bounded re-read 2026-10-17 (round breadth-first-14): protocol lines 164–297
(identity, workflow overview), 1166–1290 (tagmentation and RA attachment),
1291–1403 (precipitation clean-up) and 1453–1480 (library loading mix); CHTD lines
248–251, 379–423 and 2068–2126.

## 3. Molecular workflow

🟢 Six million cells are extracted to uHMW DNA. The protocol mixes 750 µl DNA with
250 µl Fragmentation Mix diluted 6:244 in FRA Dilution Buffer. It incubates 10 min at
room temperature, 10 min at 75 °C and at least 10 min on ice. Five microlitres of Rapid
Adapter are added and incubated 30 min at room temperature.

🟢 The reaction must be cooled on ice before Rapid Adapter is added "to prevent
denaturing the enzyme" (protocol note after the incubation table); the enzyme is not
named. FRA Dilution Buffer (FDB) is modelled as the dilution condition.

🟢 The library is recovered by adding a Precipitation Star and 500 µl Precipitation
Buffer, rotating 20 min, aspirating the supernatant and eluting in 300 µl Elution
Buffer overnight (at least 12 h). 🟡 Precipitation is not an adapter-specific selection,
so the model keeps unadapted tagmentation products as an unresolved carried fate
rather than a physical discard; only the aspirated supernatant is discarded. Wide-bore tips and slow mixing preserve molecule length. This cleanup is part
of why the ultra-long workflow differs from ordinary rapid preparation.

🟢 The published Rapid Adapter strand is `CCTGTACTTCGTTCAGTTACGTATTGCT`.

🟢 The technical document says the transposase fragments and tags DNA, then the
motor-loaded Rapid Adapter attaches to those tags without an enzyme. The tag bases,
exact tag–adapter junction and complete motor-bearing assembly are open model
questions, not confirmed publication absences. Those regions remain role tokens, and
the published RA strand is a separate graph input rather than being placed across an
unverified attachment junction (model correction 2026-10-17: the earlier endpoint
drew the RA strand directly 5′ of the tag).

## 4. Sequencing entry

🟢 No sequencing primer binds. A motor in the Rapid Adapter meters the attached strand
through an R10.4.1 nanopore. The endpoint is a component-level relationship; the
loading mix (100 µl Sequencing Buffer UL, 10 µl Loading Solution UL, 90 µl library,
30 min at room temperature) is not modelled as a graph step.
