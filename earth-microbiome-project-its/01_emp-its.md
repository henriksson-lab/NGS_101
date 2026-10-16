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

## 5. Scope

This page documents the fusion-primer EMP implementation, not every possible ITS region or primer pair. Primer choice remains taxon- and study-dependent.
