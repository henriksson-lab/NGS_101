# DRIPc-seq — strand-specific sequencing of R-loop RNA

## Sources read

- 🟢 Sanz *et al.*, “Prevalent, dynamic, and conserved R-loop structures associate with specific epigenomic signatures in mammals,” *Molecular Cell* (2016), DOI `10.1016/j.molcel.2016.05.032` (defining paper and open author manuscript, PMCID: PMC4955522).
- 🟢 Sanz and Chédin, *Nature Protocols* (2019), DOI `10.1038/s41596-019-0159-1` (stepwise protocol, PMCID: PMC6615061).
- 🟢 Its open XML text (`PMC6615061_PMC6615061.xml.txt`) was checked at the reagents list (lines 606–621) and the example adapter-trimming FASTA (lines 2310–2325).
- 🟢 Round 2026-10-17-breadth-first-09 bounded pass: defining author manuscript JATS (`j.molcel.2016.05.032_PMC4955522.xml.txt`) Results line 61 and Experimental Procedures "DRIP-seq and DRIPc-seq" line 137; 2019 protocol procedure overview lines 403–446, reagents 528–621, steps 26–31 (lines 1020–1066), 39–58 (lines 1100–1330), Table 1 (qPCR primers only, lines 2104–2143). The 2016 Supplemental Experimental Procedures PDF remains behind the PMC download gate.

## Molecular interpretation

- 🟢 Gentle extraction, five-enzyme restriction digestion and S9.6 capture are shared with DRIP-seq.
- 🟢 2016 article (line 137): eluted DRIP material is treated with 4 U DNase I (NEB) for 35 min at 37 °C; RNA is ethanol-precipitated and reverse-transcribed with iScript RT supermix (Bio-Rad) using a mix of poly(T) and random hexamers; second-strand synthesis uses dUTP instead of dTTP; cDNA is sonicated to ~200 bp; a UDG step precedes PCR; libraries are sequenced on an Illumina HiSeq 2000. Hybrids must be >50 bp to be detected.
- 🟢 2019 protocol: DNase I 6 U (3 µL in 100 µL), 45 min 37 °C, then EDTA and 75 °C 15 min (step 26–27); glycogen/NaOAc/ethanol precipitation (28–30); iScript 25 °C 5 min, 42 °C 30 min, 85 °C 5 min (31); AMPure cleanup. Second strand (39–40): dUTP-for-dTTP dNTPs, NAD, RNase H, E. coli DNA ligase and DNA polymerase I, 16 °C 30 min — nick-translation replacement of the RNA, no added second-strand oligo. Bioruptor sonication (43), NEBNext End Repair (44), Klenow exo− dA-tailing (46), NEB Quick ligation of Illumina index adapters (48), AmpErase UNG 20 min 37 °C (50), Phusion Flash PCR with P5/P7 primers, 10–15 cycles (51–52), two-step AMPure selection 0.65× then added beads for 200–500 bp (53–56), qPCR and BioAnalyzer QC (57–58).
- 🟡 DNase I dose/time differ between the article (4 U, 35 min) and the 2019 protocol (6 U, 45 min); the page shows both and the graph does not assign a dose. The library-construction steps on the page are the 2019 protocol's; the 2016 article does not itemize them.
- 🟢 Because only the first cDNA strand survives amplification, read orientation identifies the RNA strand that participated in the original hybrid.
- 🟢 The detailed protocol specifies Illumina TruSeq Single index adapters and gives library-amplification primers P5 `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′` and P7 `5′-CAAGCAGAAGACGGCATACGAGAT-3′`. Its analysis example gives adapter-trimming motifs `ACACTCTTTCCCTACACGACGCTGTTCCATCT` and `CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT`; these are trimming references, not evidence for the complete indexed ligation oligos or their physical end configuration.
- 🟡 Read length for the 2016 HiSeq 2000 data is not stated in the article text read; the 2019 protocol mentions single-end 100 bp only as an example case (line 420). GEO GSE70189 run metadata was not inspected.
- 🔴 The exact selected TruSeq Single index oligos and the instrument's sequencing primer are not printed in the defining paper or detailed protocol; the page therefore still shows the indexed ligation ends by role.

## Final readout

The library sequences RNA-derived cDNA rather than the captured genomic DNA. That change provides higher-resolution, strand-specific R-loop mapping and makes DRIPc-seq molecularly distinct from DRIP-seq.
