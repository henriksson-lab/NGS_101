# TruePrime

**Source.** Picher et al., *Nature Communications* 7, 13296 (2016), [doi:10.1038/ncomms13296](https://doi.org/10.1038/ncomms13296), PMCID PMC5141293.

- 🟢 Figure 9 and Methods: denatured DNA is primed by TthPrimPol with principal 7–9-nt DNA primer products; phi29 DNA polymerase extends them and displaces downstream strands; new exposed ssDNA receives further primers.
- 🟢 No synthetic random primer is added to TruePrime.
- 🟢 Figure 3a uses a synthetic 3′-TTTTTTTTTTGTCCTTTTTTTTTTTTTTT-5′ template to assay primase initiation at `GTCC`. It is an assay substrate, not an amplification primer or part of the WGA workflow. The paper also gives two mycoplasma-detection PCR primers in the cell-culture Methods; those are quality-control reagents, not TruePrime primers.
- 🟢 Whole-genome amplification of purified DNA: Buffer D alkaline denaturation for 3 min, Buffer N neutralization, then 3 h at 30 °C followed by 10 min at 65 °C. This page models that purified-DNA branch; the intact-cell L2-lysis branch has its own TruePrime Single Cell WGA page.
- 🟡 A later [4basebio WGA handbook](https://www.4basebio.com/s/3x0xxx_TruePrime-WGA-Kit-1.pdf), p. 5, corroborates the purified-DNA Buffer D/N schedule and 42.5 µl amplification mix. It cannot establish the exact kit revision used for the 2016 experiments.
- 🟢 Sequencing Methods: amplified DNA was Covaris-fragmented, prepared with “NebNext,” and paired-end sequenced on HiSeq 2500.
- 🟢 The paper-linked [ENA study SRP085855](https://www.ebi.ac.uk/ena/browser/api/xml/SRP085855) and its four experiment records distinguish TruePrime from comparator WGA runs but give no downstream library-construction protocol or NEBNext kit name.
- 🔴 The NEBNext kit generation and adapter oligos are not identified. The schematic ends the named protocol at amplified DNA and does not invent primer sites.

TruePrime is whole-genome amplification, not the downstream sequencing-library preparation.
