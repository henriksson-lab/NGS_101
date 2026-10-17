# RiboMethSeq

**Defining source:** Birkedal *et al.*, *Angewandte Chemie* (2015),
doi:10.1002/anie.201408362. **Illumina implementation:** Marchand *et al.*, *Nucleic
Acids Research* (2016), doi:10.1093/nar/gkw547. The open full method text was read.

## Evidence used

- 🟢 The defining method establishes protection of the phosphodiester bond immediately
  3′ of a 2′-O-methylated ribose and inference from depleted fragment endpoints.
- 🟢 The Illumina method compares alkaline, MgCl₂ and ZnCl₂ fragmentation and selects
  alkaline hydrolysis for its precision.
- 🟢 It specifies 3′ dephosphorylation, 5′ phosphorylation, NEBNext Small RNA library
  preparation, indexed amplification, paired-end MiSeq and SR50 HiSeq 1000 alternatives.
- 🟡 Only the public NEBNext 3′ adapter is written base-by-base; other kit oligos remain roles.

## Oligo source boundary

🔴 The **defining** Birkedal method used an Ion Torrent library with a specially prepared 5′ adapter carrying a 2′,3′-cyclic phosphate and a mutant RNA ligase, according to Marchand *et al.*'s comparison (PMC5027498, “Comparison of fragmentation approaches”). The original Wiley Supporting Information is listed as `anie_201408362_sm_miscellaneous_information.pdf`, but the available source fetcher cannot retrieve it; its original adapter and primer sequences therefore remain unverified. Do not transfer the later NEBNext Illumina oligos back onto that defining Ion Torrent chemistry.

🟢 Marchand *et al.*'s **later Illumina implementation** uses direct 3′ and 5′ ligation with the NEBNext Small RNA Library kit after end repair (PMC5027498, “Comparison of fragmentation approaches”). Its supplementary methods also print `5′-GTTGAATTTCTTCACTTTGACATTCAGAGCACTGGGCAGAAATCACATTG-3′` as a 25S-rRNA capture oligo for the validation assay. That capture oligo is not part of the sequenced library construct.

## Construct logic

The informative absence is local: protected bonds produce fewer 5′ and 3′ fragment ends.
End repair must occur after hydrolysis, and direct adapter ligation is essential because
tailing or template-switching would blur the nucleotide-resolution cleavage coordinate.

## Bounded curation round 4 (2026-10-16)

Inspected: Related published Marchand2016 PMC5027498 JATS1–29 (identity),62–78 (fragmentation, end repair, library and reads).

Made explicit that the displayed NEBNext/Illumina recipe is Marchand2016, not a verified reconstruction of the catalogue-defining Birkedal2015 protocol. Named Antarctic-phosphatase3-prime repair followed by T4PNK/ATP5-prime phosphorylation.

Deferred: Defining 2015 article/supplements remain unavailable in this bundle; protocol-boundary/catalogue correction and exact version-specific oligos require a separate pass. This pass does not establish a complete source review or publication-wide absence.
