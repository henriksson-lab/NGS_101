# R-ChIP — in-vivo R-loop capture with inactive RNASEH1

## Sources read

- 🟢 Chen *et al.*, “R-ChIP Using Inactive RNase H Reveals Dynamic Coupling of R-loops with Transcriptional Pausing at Gene Promoters,” *Molecular Cell* (2017), DOI `10.1016/j.molcel.2017.10.008` (defining paper; PMCID: PMC5957070). STAR Methods “Strand-specific R-ChIP Library Construction” read in the PMC author-manuscript XML text (`j.molcel.2017.10.008_PMC5957070.xml.txt`, lines 299–303) and “R-ChIP Data Analyses” (lines 309–311). The PMC supplement DOCX (NIHMS964699) is behind the PMC download gate and was not read.
- 🟢 Chen *et al.*, “R-ChIP for genome-wide mapping of R-loops by using catalytically inactive RNASEH1,” *Nature Protocols* (2019), DOI `10.1038/s41596-019-0154-6` (authors' stepwise protocol; PMCID: PMC6604627). Read in `PMC6604627_PMC6604627.xml.txt`: experimental-design read-length passage (line 139), primer extension, dA tailing, adaptor preparation/ligation, PCR and library recovery, sequencing (lines 604–766), Table 1 (lines 1475–1485).

## Molecular interpretation

- 🟢 Cells stably express V5-tagged catalytic mutant RNASEH1 D210N; the binding-defective WKKD mutant supplies the negative control. Cells are fixed in 1% formaldehyde (10–15 min, RT; glycine quench), nuclei are lysed and chromatin is sonicated to 250–600 bp before overnight anti-V5 bead capture and TSE I/II/III and TE washes (2017 Methods).
- 🟢 Elution in SDS buffer, overnight 65 °C de-cross-linking, RNase A then proteinase K, phenol and phenol:chloroform:isoamyl alcohol extraction, and ethanol precipitation recover the “precipitated hybrid fragment” (2017 Methods).
- 🟢 Captured DNA is heat-denatured (95 °C 5 min, 25 °C 5 min) with the tailed N9 primer `5′-/invddt/CAAGCAGAAGACGGCATACGAGNNNNNNNNN-3′` and extended with phi29 DNA polymerase (30 °C 20 min; 65 °C 10 min) (2019 protocol; the 2017 Methods prints the same primer).
- 🟢 Klenow fragment (3′→5′ exo−) with dATP (37 °C 30 min) adds the A. The 2017 Methods states the adapter ligates only to one end because the other end carries the 5′ overhang introduced by the N9 primer; PCR therefore preserves strand orientation.
- 🟢 Adapter: oligo A `5′-/Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-3′` and oligo B `5′-AGACGTGTGCTCTTCCGATCT-3′`, annealed 95 °C → 25 °C at 0.1 °C/s; ligated with T4 DNA ligase for 2 h at RT. 🟡 Their 12-base `GATCGGAAGAGC` overlap leaves oligo B's 3′-T over the dA, and oligo A's 5′ phosphate joins the copy strand's 3′-A. Whether oligo B's 3′-T is sealed to the template's 5′ end depends on that end's phosphorylation, which is not stated.
- 🟢 PCR with Phusion, the PCR primer `5′-CAAGCAGAAGACGGCATACGAG-3′` and the barcode primer `5′-AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′` (`NNNNN` = flexible barcode). 🟡 The barcode primer is P5 + 5-nt barcode + the TruSeq Read 1 sequence (checked against `lib/illumina.py` in the module); its 3′ 33 nt anneal to the reverse complement of dA + oligo A, and the PCR primer repeats the N9-primer tail.
- 🟢 Products are gel-isolated and sequenced on the HiSeq 2500 with the Illumina Read 1 sequencing primer, single-end, 40 nt (2019 protocol line 139 and sequencing section). The protocol's trimming command `cutadapt –a CTCGTATGCCGTCTTCTGCTTG` (line 762) removes the reverse complement of the N9-primer tail from the 3′ end of reads, consistent with Read 1 starting at the ligated adapter.

## Published discrepancy shown on the page

- 2017 Methods: 16–18 PCR cycles, 130–350 bp gel-isolated. 2019 protocol: test PCR at 14 and 16 cycles (18 if weak), 140–350 bp excised. The schematic shows both and chooses neither.

## Final readout

🟡 The final library (copy-strand orientation) is tailed-N9 handle (P7-side) – N9 – insert – A – oligo A – 5-nt barcode – P5. The library is reconstructed from the printed oligos; seqprimers places Read 1 on dA + oligo A, so Read 1 reads into the insert from the ligated end. This gives the sequence of the captured template strand. The sequenced fragment is the DNA component of an in-vivo RNASEH1-bound R-loop. It is a strand-specific, one-ended library and is not interchangeable with S9.6 DRIP-seq or RNA-component DRIPc-seq.

## Open questions

- How the 5-nt barcode, located in the i5 position between P5 and Read 1, was read on a single-end HiSeq 2500 run is not stated in the inspected passages.
- The PMC supplement DOCX of the 2017 paper (gated) and the cited ChIP predecessor (Ji *et al.* 2013) were not read in the capped pass.
