# SHAPE-Seq

**Defining source:** Lucks *et al.*, *PNAS* (2011),
doi:10.1073/pnas.1106501108.

## Evidence used

- 🟢 Each designed RNA carries a DNA barcode at its 3′ end before pooling.
- 🟢 SHAPE chemistry acylates flexible 2′ hydroxyls; primer-extension stops encode the
  modification sites relative to a minus-reagent control.
- 🟢 A second adapter is ligated to cDNA and paired-end sequencing connects each stop
  with the barcode identifying the RNA template.
- 🟢 The original supplementary methods were recovered from the CaltechAUTHORS mirror.
  They print the v1 A_adapter_t, A_adapter_b, treated/control handle patterns and PCR
  construction; the shared model now uses those sequences.

## Published adapter and primer architecture

🟢 `A_adapter_t` / `PE_forward` is
`AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
It is the P5 arm followed by the Read 1 site and forms the 5′ tail of each RT primer.

🟢 The RT primers then carry a four-base treatment handle—`RRRY` for the 1M7-treated
branch or `YYYR` for the solvent control—followed by a target-specific RT-primer sequence.
There is therefore no single universal terminal RT-primer sequence to assert.

🟢 The single-stranded ligation oligo `A_adapter_b` is
`/5Phos/AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG/3SpC3/`.
Its 5′ phosphate permits ligation to cDNA and its 3′ spacer prevents concatemerization.

🟢 The defining experiment used 2×50 cycles on an Illumina Genome Analyzer IIx.
Read 1 begins with the treatment handle and RNA-identity end; the paired read reports the
SHAPE-dependent primer-extension stop from the opposite end.

🔴 The supplement does not identify a separate, immutable ordering sequence for the
historical GAIIx paired-end sequencing-primer reagent. The exact adapter landing regions are
shown, but that instrument reagent is not silently replaced with a modern TruSeq primer.

## Construct logic

The original protocol is a multiplexed in-vitro accessibility assay for defined RNA
pools. It is distinct from SHAPE-MaP, where modifications are read as mutations rather
than termination sites.
