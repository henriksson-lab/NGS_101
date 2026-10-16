# Human NET-seq

Human NET-seq recovers chromatin-associated nascent RNA and sequences its native 3′ end. The first aligned base after removal of the molecular barcode reports polymerase position.

## Original sources and selected version

🟢 [Mayer et al., Cell (2015)](https://doi.org/10.1016/j.cell.2015.03.010) and its [author-hosted original with Extended Experimental Procedures](https://www.cecerelab.com/_files/ugd/594009_59c7fec32bcc42cf8c824804675fe319.pdf) define the displayed human implementation. Current extracted methods are at lines 1075–1100.

🟢 That source explicitly inherits library construction from [Churchman and Weissman (2012)](https://doi.org/10.1002/0471142727.mb0414s98), recovered as an [original PDF from the official author website](https://churchman.med.harvard.edu/pdfs/2012-native-elongating-transcript-sequencing-net-seq.pdf). Its yeast immunoprecipitation front end is not assigned to human fractionation. The source also cites [Churchman and Weissman (2011)](https://doi.org/10.1038/nature09652).

🟢 [Mayer and Churchman (2016)](https://doi.org/10.1038/nprot.2016.047), PMC5326677, is a later detailed human protocol, not the defining 2015 experiment. Earlier repository citations incorrectly linked its title to `nprot.2016.086`, which is a different protocol. This was a repository citation error.

## Core oligos

The 2015 original prints these ordering strings with modifications:

- Linker: `AppNNNNNNCTGTAGGCACCATCAAT/3ddC`.
- oLSC007: `Phos/ATCTCGTATGCCGTCTTCTGCTTG/iSp18/CACTCA/iSp18/TCCGACGATCATTGATGGTGCCTACAG`.
- oLSC006: `TCCGACGATCATTGATGGTGCCTACAG`.

The 2012 inherited original Table 4.14.1 prints reverse PCR primer oNTI231, `CAAGCAGAAGACGGCATACGA`, full oNTI230 and indexed oMHL001/002 alternatives. These are positive published primer data. They do not identify the exact forward member selected in 2015 after its RT-primer changes. The declaration keeps that outer arm and selected index unassigned; this is model data awaiting assignment, not accepted publication absence.

🟡 Circularization retains the first three bases of the 2015 RT primer before the reverse-primer-complement tract. The earlier endpoint silently omitted them. The final construct now preserves the complete source-known RT5 tract across the closure junction. Biological inserts are variable-length roles; the 28-nt oGAB11 sequence is a separate control, not the length of every nascent RNA.

## Source-specific reaction route

🟢 The 2015 route ligates the preadenylated, blocked N6 linker before partial alkaline fragmentation and selects 35–100-nt RNA. It then reverse-transcribes, circularizes, depletes nineteen abundant mature-RNA targets by biotinylated probe hybridization and streptavidin capture, retains the supernatant and amplifies.

The inherited 2012 library procedure supplies truncated T4 RNA ligase 2, SuperScript III, post-RT NaOH removal of RNA, cDNA PAGE, CircLigase, Phusion PCR and post-PCR gel selection. Its cDNA selection is described relative to its own RT-primer/control bands. The 85–160-nt cDNA interval belongs to the 2016 successor; it is not assigned to the selected 2015 branch. Earlier repository prose attributed this numerical interval to 2015 incorrectly.

The inherited procedure warns that RT-primer carryover can yield a separate PCR band. The graph represents possible primer-derived background and its gel-exclusion route without assigning an exact 2015 molecular structure or complete removal. The later 2016 protocol explicitly identifies empty-circle background near 120 nt; that exact size is not transferred to the 2015 experiment.

## Readout and filtering

🟢 The 2015 original reports NextSeq 500 single-read 75 nt and HiSeq 2000 single-read 50 nt, both with oLSC006. The displayed run selects the former. The first six acquired bases are the molecular barcode, followed by sequence copied from the nascent RNA end. Index acquisition and the selected sample-index member are not assigned.

The defining article's processing section computationally removes barcode-matched mispriming, PCR duplicates and reads exactly at selected splice-intermediate positions. These analysis filters do not imply that the corresponding physical library molecules were discarded during preparation.

## Remaining source obligations

The 2015 Table S1 is explicitly a separate original, PMC SD6 / `NIHMS680477-supplement-6.pdf` (65.9 KB). Current PMC download routes are gated and the PMC Cloud TDM record has no attachments. The recovered author PDF cites the table but does not contain it. The twenty depletion oligos in the 2016 successor are not substituted for the nineteen-target 2015 set. This is an access blocker, not confirmed absence.

The source-review checklist preserves exact located readings and named predecessor, supplement, data and vendor scope still awaiting inspection. No missing-information claim is accepted from that incomplete inventory.
