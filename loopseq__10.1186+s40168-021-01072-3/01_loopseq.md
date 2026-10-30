# LoopSeq — intramolecular barcode distribution

Callahan et al. (2021), *Ultra-accurate microbial amplicon sequencing with synthetic long reads*, [10.1186/s40168-021-01072-3](https://doi.org/10.1186/s40168-021-01072-3).

## Inspected sources and version

Two originals: published PMC8179091 article, Introduction and **Sequencing library preparation** with its full-length 16S subsection; and its directly linked vendor manual **LoopSeq 16S Microbiome SSC 24-Plex Kit / LoopSeq 16S & 18S Microbiome SSC 24-Plex Kit, With Single Sample Calibration**, printed **Version 1.7**, revision history **July 2019**. The archived original is [the article-linked v1.7 manual](https://github.com/benjjneb/LoopManuscript/blob/master/Docs/LoopSeq_16S_16S18S_Microbiome_SSC_Manual_v1p7.pdf), fetched through the standard source tool. Inspected manual: cover, revision identity and Sections I–III (pp. 7–8); Section IV recipe order/conditions (p. 9) were checked and represented operationally, without complete polymer geometry. The standard full-bundle oligo scraper was run. No current vendor revision replaces this historical original. Retrieved superseded preprint and unrelated cited protocols are not active chemistry evidence.

This is one technology entry; the paper's 16S, 18S-ITS and bacterial-genome applications are not counted as three preparation methods. The local graph is scoped to the historically linked **SSC microbiome v1.7** procedure, not every kit evaluated in the article.

## Source operation and distinction

🟢 Published Methods describes attaching a parent-molecule UMI and a sample-specific tag, amplification/multiplexing, intramolecular distribution of each UMI to random positions, and fragmentation at these positions. UMI-tagged fragments average about 400 bp and undergo Illumina PE150 sequencing; assembly groups reads by sample then UMI. The computation reconstructs a long sequence, not a physically reassembled long DNA molecule.

🟢 This differs from existing TELL-seq/stLFR bead co-barcoding of linked fragments. Intramolecular barcode redistribution is stated by this original; the model does not substitute MuA, Tn5, nicking or circularization as its proprietary mechanism. Barcode widths and kit oligo/mix identities remain unassembled, not confirmed absent. Structured molecule-UMI and sample-tag meanings are displayed as separate components without invented bases or joined strands.

## Historical v1.7 workflow

🟢 Section 1.2 adds 5 µL diluted genomic DNA to 15 µL predispensed Barcoding Adapters, then performs two assignment-PCR cycles. Section 1.3 adds Inactivation Enzyme SSC; Section 1.4 performs SPRIselect purification. The graph treats the kit PCR as a coarse operational transformation because the internal primer species and exact molecular geometry have not been assembled. It does not invent two independent primers or a symmetric primer to satisfy a structural PCR checker.

🟢 A separate aliquot is used for Section II calibration qPCR. Calibration determines the amount carried into the next reaction; its analytical PCR products do not enter the sequencing library. Remaining material and untransferred aliquots have explicit unresolved fates rather than assumed disposal.

🟢 Section 3.1 amplifies calibrated individual samples for 22 cycles. Section 3.2 physically pools **5 µL from each of 24 wells**, performs SPRIselect cleanup and elutes into 15 µL. Two symbolic sample classes show the pool without inventing experimental identities or implying a two-sample kit.

🟢 Section 3.3 adds Distribution Mix and Distribution Enzyme, incubates at 20°C for 15 min and then 75°C for 5 min. Section 3.4 adds Activation Mix/Enzyme, incubates at 20°C for 2 h, adds Neutralization Enzyme and incubates at 37°C for 15 min. Section 3.5 performs SPRIselect cleanup. These source labels do not establish an atom-level enzyme mechanism. Residual reagents/products remain with reaction material until the physical cleanup; unbound supernatant/washes are removed there, and DNA-containing eluate is retained.

🟢 Section 4.1 subsequently specifies combined fragmentation/end repair/A-tailing (32°C for 5 min, 65°C for 30 min); Section 4.2 performs adapter ligation at 20°C for 15 min. Section 4.3 performs SPRIselect cleanup, Section 4.4 performs 12-cycle library-index PCR with one of P1–P4, and Section 4.5 performs SPRIselect cleanup. These operational steps are represented through final purification, while full polymer/end geometry and exact final readout remain deferred. The model does not infer an exact adapter T overhang from the source's A-tailing label. The proprietary index-primer mix is not modeled as an invented pair or symmetric primer species.

## Lineage and deferred review

The article cites BAsE-Seq (Hong et al. 2014, [10.1186/s13059-014-0517-9](https://doi.org/10.1186/s13059-014-0517-9)) and Stapleton et al. (2016), *Haplotype-phased synthetic long reads from short-read sequencing*, [10.1371/journal.pone.0147229](https://doi.org/10.1371/journal.pone.0147229), as earlier academic work. This is a documented broad conceptual ancestry statement, **not** evidence of exact shared reagents, constructs or recipes. Those originals remain unreviewed; no constructor inheritance or exact component lineage is declared.

Deferred: historical manual's exact kit polymers/mix identities and end/ligation geometry, Section IV molecular assembly/read plan, other application-specific kit manuals, chemistry-relevant directly cited originals, remaining main/SI sections and independent review. This capped pass does not establish complete source review or accepted missing information. No finished schematic is published.
