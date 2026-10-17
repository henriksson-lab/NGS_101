# TELL-seq — research notes

## Defining source

Chen Z, Pham L, Wu T-C, et al. *Ultralow-input single-tube linked-read library method enables short-read second-generation sequencing systems to routinely generate highly accurate and economical long-range sequencing information.* Genome Research 2020;30:898–909. [doi:10.1101/gr.260380.119](https://doi.org/10.1101/gr.260380.119) (open full text PMC7370886).

## Sources read

- 🟢 Bounded pass 2026-10-17 (round breadth-first-20) used the published version. The PMC7370886 JATS and PDF were fetched with `python3 tools/get_sources.py --slugs tell-seq__10.1101+gr.260380.119`. Sections read: Results "Clonal barcoding with simultaneous capture and strand transfer reactions" and the Figure 1 legend (JATS lines 80–88); Figure 1A/B, inspected in the rendered PDF on p. 899; Methods "TELL-seq library construction and sequencing" and "Primary sequencing data processing" (JATS lines 178–188).
- 🟢 Supplemental Figure S1 was inspected in the rendered original. So was the Supplemental Code Tell-Read user guide, section "Sample index names" (pp. 11–12).
- The bioRxiv 852947 v1 preprint is superseded and is not used for current claims.
- Not re-read in this pass: the Universal Sequencing Technology TELL-Seq WGS Library Prep and Library Sequencing user guides. They are not in the source bundle.

## Evidence record

- 🟢 0.1–5 ng genomic DNA, 3–10 million barcoded TELL beads, non-barcoded beads and transpososomes are mixed in one 0.2-mL PCR tube (Figure 1 legend). Methods give 0.1 or 0.5 ng microbial DNA with about 3 million beads, and 5 ng human DNA with about 8 million beads.
- 🟢 MuA strand-transfer complexes (STCs) are captured onto bead barcode oligos while strand transfer is still happening. The Results describe this as strand transfer plus hybridization. Non-barcoded spacer beads and higher viscosity keep the barcoding clonal.
- 🟡 Supplemental Fig. S1A Method 1 lists transpososomes + barcoded beads + ligase + genomic DNA in one reaction. The Results text maps this simultaneous method to TELL-seq. The model therefore lists a ligase in the barcoding step, but the ligated junction is not described.
- 🟢 A second transpososome tags DNA between STCs and introduces the second amplification priming site. STCs are then broken and the magnetic beads washed. Methods: "TELL beads with captured barcoded DNA were amplified": 13–14 cycles for microbial samples, 8 for human. P5 and P7 are added in this PCR.
- 🟢 TELL beads are 3 µm and carry about 50,000 barcode templates each. About 1% of templates share a barcode with another bead. The molecular barcode is 18 degenerate bases (maximum homopolymer 6) drawn from a pool of more than 2.4 billion.
- 🟢 Sequencing: 2 × 146 paired-end reads, an 18-cycle Index 1 read (TELL barcode) and an 8-cycle Index 2 read (sample index), run on MiSeq/NextSeq (microbial) or NovaSeq (human) "based on the manufacturer's protocols". Reads are demultiplexed by I2, and I1 reads are the TELL barcodes.
- 🟢 The Tell-Read guide lists Index 2 primer names and 8-base sequences T500–T524, e.g. T501 `TGAACCTT`, T502 `TGCTAAGT`, and T500 `NNNNNNNN`. It does not say whether these are written in read orientation or primer orientation. The model keeps the sample index as a placeholder.
- 🔴 Neither the paper nor the inspected supplements print the bead adaptor, the transposon ends, the PCR primers or the Read 1/Read 2 region sequences. The schematic does not imply standard TruSeq or Nextera primer binding.
- ⚪ Earlier notes said that custom R1, R2 and I1 primers are required and that I2 needs one only on some instruments. That came from the vendor *TELL-Seq Library Sequencing User Guide* Rev B, which is not in the bundle and was not re-read in this pass.

## Final construct

Figure 1B, read in the rendered original, gives the order P5 · 8-base sample index (Index 2) · Read 1 region · 100–400 bp insert · Read 2 region · 18-base TELL barcode (Index 1) · P7. In the arrows, Index 1 reads from the Read 2 region into the barcode and Index 2 reads from the Read 1 region into the sample index. In Figure 1A the bead barcode oligo is joined to the first (blue) transposon end, which is the Read 2 region. The second (red) transpososome supplies the Read 1 region.

**Model correction (2026-10-17):** the earlier model drew the 8-nt sample index on the P7 side after the barcode, with separate invented Index 1/Index 2 primer sites. It also swapped the Read 1/Read 2 roles of the two transposon ends. Both now follow Figure 1B.

Other changes in the same pass:
- After STC disruption the fragments stay on the beads. The bead wash is now the represented separation step.
- Fragments lacking one of the two different priming ends are now an unresolved candidate, not a discard.

Reads that share an I1 barcode are linked back to one long input molecule. The short genomic fragments themselves remain ordinary paired-end inserts.

## Open questions

- STC-disruption condition, wash buffer and PCR polymerase (kit reagents). These need the TELL-Seq WGS Library Prep user guide (deferred).
- Which ligase is used, and the exact junction it makes on the bead (Supplemental Fig. S1 only names "ligase").
- Bead adaptor, transposon-end and PCR-primer sequences, and whether each read uses a custom primer.
