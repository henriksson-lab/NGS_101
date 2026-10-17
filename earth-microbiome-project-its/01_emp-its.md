# Earth Microbiome Project ITS amplicon sequencing

## 1. What it is

The Earth Microbiome Project fungal ITS workflow amplifies ITS1 with primers that already contain complete Illumina arms and a twelve-base Golay sample barcode. Unlike the two-stage Illumina 16S workflow, this is a fusion-primer library and uses locus-specific custom sequencing primers extending beyond the PCR-primer boundary.

## 2. Authoritative source

🟢 Earth Microbiome Project, *ITS Illumina Amplicon Protocol* (EMP.ITSkabir), and the linked November 2016 `ITS1f_ITS2` oligo workbook. The workbook publishes the fusion constructs and all three custom sequencing primers.

## 3. Published fusion primers

🟢 Forward ITS1f construct, 5′→3′:

`AATGATACGGCGACCACCGAGATCTACAC GG CTTGGTCATTTAGAGGAAGTAA`

This is the P5-side Illumina arm, a `GG` linker and ITS1f.

🟢 Reverse ITS2 construct, 5′→3′:

`CAAGCAGAAGACGGCATACGAGAT NNNNNNNNNNNN CG GCTGCGTTCTTCATCGATGC`

This is the P7-side arm, a twelve-base Golay barcode, a `CG` linker and ITS2. Spaces above only expose the published components.

🟢 **Source discrepancy.** The [EMP ITS page](https://earthmicrobiome.org/protocols-and-standards/its/), reverse-primer construct, and its [2018 protocols.io version](https://www.protocols.io/view/emp-its-illumina-amplicon-protocol-pa7dihn.pdf), p. 4, show `NNNNNNNNNN` (10 bases). The linked [November 2016 primer workbook](https://earthmicrobiome.org/wp-content/uploads/2022/06/ITS1f_ITS2_Nov2016.xlsx), *Reverse Primer Constructs* sheet, lines 11–971 of the extracted text, lists 960 twelve-base Golay members. [Smith and Peay 2014](https://web.stanford.edu/~kpeay/papers/Smith%26Peay_PLOS_2014.pdf), *Molecular methods*, p. 3, also specifies a 12-base error-correcting Golay barcode and a 12-base index read. The model uses the ordering workbook’s 12 bases; whether the 10-N placeholder is shorthand or a distinct unlisted construct remains unknown.

## 4. Read geometry

🟢 Read 1: `TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC`.

🟢 Read 2: `CGTTCTTCATCGATGCVAGARCCAAGAGATC` (`V` and `R` are degenerate positions).

🟢 Index 1: `TCTCGCATCGATGAAGAACGCAGCCG`.

These produce paired-end reads plus one custom index read over the twelve-base Golay barcode; there is no i5 read. For an exact computed landing, the schematic uses the valid Read 2 pool member `V=G, R=A` and labels that choice rather than treating it as inferred chemistry.

🟢 The cited [Caporaso et al. 2012 Supplementary Methods](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fismej.2012.8/MediaObjects/41396_2012_BFismej20128_MOESM84_ESM.doc), *MiSeq preparation and sequencing protocol*, describes a **historical 16S** run: 12 index cycles, custom Read 1/index/Read 2 primers in cartridge reservoirs 12/13/14, and 151×12×151 cycles (extracted lines 229–287). It recommends 15–30% PhiX for its 16S amplicons (lines 208–218); the later EMP ITS page advises 5–10% for ITS. The 2012 instrument settings are an example of how the custom-primer architecture was run, not a documented cycle/PhiX selection for every EMP ITS run.

## 4a. Amplification and pooling

🟢 EMP ITS page (retrieved 2026-10-17), *PCR reaction mixtures*, *Thermocycler conditions* and *Amplification protocol*; the same steps appear in the 2018 [protocols.io version](https://dx.doi.org/10.17504/protocols.io.pa7dihn), pp. 6 and 10–11:

- 25 µl reactions: 13 µl water, 10 µl Platinum Hot Start PCR Master Mix (2x; 0.8x final), 0.5 µl each 10 µM primer (0.2 µM final), 1 µl template.
- 94 °C 1 min; 35 cycles of 94 °C 30 s, 52 °C 30 s, 68 °C 30 s; 68 °C 10 min; 4 °C hold (current web page table; the protocols.io thermocycler table is an image and was not transcribed).
- Each sample is amplified in triplicate and the replicates are pooled to 75 µl; samples are checked on an agarose gel, quantified with PicoGreen, and combined at 240 ng each.
- The pool is cleaned with the MoBio UltraClean PCR Clean-Up Kit; optionally half the pool is gel-extracted when spurious bands were seen. Its A260/A280 should be 1.8–2.0. The pool is then sent for sequencing with the three custom primers, and 5–10% PhiX is added.

🟢 **Source discrepancy (amplicon size).** The current EMP page gives "Amplicon size: ~250–600 bp (Bokulich & Mills, 2013; Hoggart et al., 2018)" under *Thermocycler conditions*, but "Expected band size for ITS1f-ITS2 is ~230 bp" in *Amplification protocol* step 3. The 2018 protocols.io version gives ~230 bp in both places. The model keeps the fungal insert as a variable placeholder and chooses neither value.

## 5. Scope

This page documents the fusion-primer EMP implementation, not every possible ITS region or primer pair. Primer choice remains taxon- and study-dependent.

## Bounded pass 2026-10-17 (round breadth-first-09)

Re-fetched and read the EMP ITS web page and the 2018 protocols.io PDF; re-checked the workbook *Forward and Sequencing Primers* sheet rows 4–9 and the first *Reverse Primer Constructs* rows. The graph now names the Platinum Hot Start master mix and the cycling program, the triplicate PCR, equal-mass pooling and the UltraClean clean-up. The page shows the amplicon-size discrepancy. Workbook sheet labels were corrected; they were previously cited as non-existent "Read1/Read2/Index" sheets. Open questions: acquired cycle counts for EMP ITS runs (the page defers to Caporaso 2012, whose settings are for a historical 16S run); the 10-N versus 12-nt barcode discrepancy; Smith & Peay 2014 Figure S1 and Walters et al. 2016 Table S1 were not re-checked in this pass.
