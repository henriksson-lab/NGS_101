# YAMAT-seq

**Defining source:** Shigematsu *et al.*, *Nucleic Acids Research* (2017),
doi:10.1093/nar/gkx005. The PMC full text and supplement were read.

## Evidence used

- 🟢 Methods lines 52–60 publish all five Y-adapter oligos, their DNA/RNA composition,
  end modifications, Rnl2 annealing/ligation, TruSeq Small RNA RT-PCR and 100-nt
  single-read sequencing.
- 🟡 **Source discrepancy:** Methods, “Y-shaped adapter design” (PMC JATS line 52),
  prints Y-3′-AD as `5′-/5Phos/GTATCCAGTTGGAATTCTCGGGTGCCAAGG/3ddC/-3′`,
  with a blocked 3′ end. Results, “Y-AD” description (line 88), calls the
  same strand 3′-OH. The model follows the explicit ordered oligo `/3ddC/`;
  the physical oligo used cannot be verified from these conflicting statements.
- 🟢 Lines 80–88 establish deacylation, recognition of the mature-tRNA 5′-NCCA-3′
  overhang and sealing of both acceptor-stem nicks.
- 🟢 The paper names the TruSeq Small RNA kit but does not reproduce its RT/PCR
  oligos. Illumina's [TruSeq Small RNA adapter-sequence sheet](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm)
  gives RTP as 5′-GCCTTGGCACCCGAGAATTCCA-3′, RP1 as
  5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA-3′, and the
  indexed PCR primer family as
  5′-CAAGCAGAAGACGGCATACGAGAT[6-base index]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′.
  The same sheet lists each index sequence. The paper does not say which index was
  chosen for each sample, so a specific six-base index remains unspecified.
- 🟡 The current page shows the final flow-cell ends by role. Its exact sequencing
  primer placement remains unmodeled; the Illumina sheet establishes RT and PCR
  oligos but does not identify the HiSeq 2500 run primer by sequence.
- 🟢 The additional nucleotide sequences in the paper's methods (lines 64–68 of the
  PMC full-text extraction) are a synthetic tRNA substrate and RT/qPCR primers for
  **ligation-efficiency experiments**. They are not the library amplification primers
  used for the sequencing runs; that step explicitly uses primers supplied with the
  TruSeq Small RNA kit (line 60). This distinction prevents the assay-control primers
  from being substituted for the unresolved sequencing-primer landing site.

## Construct logic

The four Y-5′ variants differ only in their terminal ribodinucleotide and collectively
pair with all discriminator-base/NCCA ends. Y-3′-AD is 5′ phosphorylated and 3′ blocked.
Annealing creates two Rnl2 substrates; amplification is possible only after both relevant
adapter–tRNA boundaries have been sealed.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5605243 JATS46–60, deacylation, adapter design, ligation, RT and sequencing.

Specified deacylation at pH9/37°C/40min, 37°C1h plus overnight4°C ligation, SuperScriptIII55°C60min and 11-cycle library PCR; recorded single-read100. Printed /3ddC/ and the existing Methods/Results discrepancy notice retained.

Deferred: Original Results discrepancy reread, selected kit primer geometry, exact acceptor-stem annealing and gel-selection intermediates; not a new publication error finding. This pass does not establish a complete source review or publication-wide absence.

Corrected the adapter-set state to five separate oligo strands instead of one concatenation, fixing a real source-guard failure.
