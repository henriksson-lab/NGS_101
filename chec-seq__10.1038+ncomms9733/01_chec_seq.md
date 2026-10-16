# ChEC-seq — research notes

Zentner et al., *Nature Communications* (2015),
[doi:10.1038/ncomms9733](https://doi.org/10.1038/ncomms9733).

## Sources read

- 🟢 The defining article and supplementary information (PMCID: PMC4618392).
- 🟢 Zentner et al., *Genome-wide Mapping of Protein-DNA Interactions with ChEC-seq
  in Saccharomyces cerevisiae* (PMCID: PMC5608238), for the detailed experimental
  sequence and controls.
- 🟢 Orsi et al., *Mapping regulatory factors by immunoprecipitation from native
  chromatin* (PMCID: PMC4410783), Basic Protocol 3, for the exact paired-end adapter
  and enrichment-primer oligos in a separately documented Henikoff-laboratory
  low-input library variant. The defining ChEC paper does not assign a concrete
  sample-index member to each published library.
- 🟢 GEO series GSE67453, including GSM1647311, for the extraction, size-selection and
  25-base paired-end run metadata of the defining data set.

## Molecular logic

1. 🟢 A yeast chromatin protein is genetically fused to 3×FLAG-MNase. Nuclear calcium is
   normally too low for strong MNase activity, so the fusion can bind chromatin before
   the timed cleavage reaction.
2. 🟢 Cells are washed and permeabilized for 5 min at 30 °C in 0.1% digitonin. Calcium is
   then added to 2 mM, activating tethered MNase in situ. Aliquots are withdrawn across
   short time points and stopped with EDTA, EGTA, high salt and SDS.
3. 🟢 Rapid cuts are enriched around sites occupied by the fused factor. Two nearby cuts
   release a small DNA fragment containing the protected binding region. Later cleavage
   increasingly includes transient encounters and bulk accessibility.
4. 🟢 Nucleic acid is phenol/chloroform-extracted, ethanol-precipitated and RNase-treated.
   A 2.5:1 AMPure-bead:sample ratio binds large genomic DNA; the unbound fraction carrying
   small ChEC fragments is recovered for library preparation.
5. 🟢 Ends are repaired with T4 DNA polymerase, Klenow and T4 PNK, then Klenow exo-minus
   adds one 3′ dA. Annealed single-index TruSeq adapters are ligated and the library is
   enriched with KAPA polymerase. This end repair removes the biochemical overhang but
   preserves the genomic coordinate of each MNase boundary.

## Exact library oligos

The cited low-input procedure publishes the ordering forms, including modifications:

- Universal Adapter:
  `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`
- Indexed Adapter:
  `5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′`
- P5 enrichment primer: `5′-AATGATACGGCGACCACCGAG-3′`
- P7 enrichment primer: `5′-CAAGCAGAAGACGGCATACGAG-3′`

🟢 `*` is the terminal phosphorothioate linkage, `/5Phos/` is the indexed adapter's
5′ phosphate, and `NNNNNN` is the six-base i7 sample barcode. The defining construct is
single-indexed; an i5 read is not part of this library.

## Read interpretation and controls

🟢 The defining data set used 25-base paired-end reads on HiSeq 2500. Both mate starts
mark repaired MNase boundaries around the selected short fragment. The paired interval,
not merely local coverage, therefore represents the two-cut ChEC product.

🟢 A soluble nuclear MNase strain under the matched promoter measures MNase sequence bias
and accessible-chromatin background. No-calcium and early-time samples distinguish basal
damage and rapid tether-proximal cleavage. ChEC-seq is antibody-free and unfixed, but it
requires an expressed fusion whose residence and linker geometry may affect cleavage.
