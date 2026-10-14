# TeloPCR-seq

## 1. What it is

TeloPCR-seq attaches an anchor to native *S. pombe* chromosome termini, amplifies each end with a subtelomeric primer, and reads the resulting telomere amplicons as PacBio circular libraries.

## 2. Sources read

- 🟢 Bennett et al. (2016), [doi:10.1002/1873-3468.12444](https://doi.org/10.1002/1873-3468.12444), PMC5561429 full text, `_data/sources/doi_10.1002+1873-3468.12444/1873-3468.12444_PMC5561429.xml.txt` line 48, Materials and Methods. The article itself prints the anchor, its terminal modifications, and both subtelomeric primers.
- 🟡 Makarova et al. (2025), [doi:10.3389/fmolb.2025.1725112](https://doi.org/10.3389/fmolb.2025.1725112), method overview.

## 3. Construct and reaction order

1. 🟢 CircLigase attaches a single-stranded anchor at the native telomere terminus.
2. 🟢 A chromosome-end-specific subtelomeric primer and the anchor primer amplify the complete telomere tract.
3. 🟢 The amplicon is converted into a circular PacBio template.

The defining article gives these ordered oligos 🟢 (5′→3′):

| Role | Oligo as printed |
|---|---|
| Ligated anchor | `5′-[Phos]AGGGTTAATAAGCGGCCGCGTCG[AmC3]-3′` |
| Subtelomeric primer 1 | `5′-CGAGGCTGCGGGTTACAAGGTTACG-3′` |
| Subtelomeric primer 2 | `5′-GGGAATTTAGGAAGTGCGGTAAGTTG-3′` |
| Anchor PCR primer | `5′-CGACGCGGCCGCTTATTAACCCT-3′` |

The anchor PCR primer is the reverse complement of the unmodified anchor bases 🟡 (computed). `[Phos]` is the ligating 5′ phosphate; `[AmC3]` blocks the anchor 3′ end. The model shows primer 1 as one representative subtelomeric boundary. The PacBio hairpin oligos and sequencing-primer site remain undisclosed 🔴; a modern SMRTbell kit sequence cannot establish which historical kit revision was used.

## 4. Final library and sequencing

The selected telomere amplicon lies between two hairpin adapters. A PacBio sequencing primer binds a hairpin primer site and the polymerase makes repeated passes around the closed template.
