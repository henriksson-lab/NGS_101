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

## Construct logic

The informative absence is local: protected bonds produce fewer 5′ and 3′ fragment ends.
End repair must occur after hydrolysis, and direct adapter ligation is essential because
tailing or template-switching would blur the nucleotide-resolution cleavage coordinate.
