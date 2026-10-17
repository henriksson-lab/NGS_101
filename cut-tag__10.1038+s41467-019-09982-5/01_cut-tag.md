# CUT&Tag: original 2019 bulk method

CUT&Tag recruits adapter-loaded Protein A–Tn5 to chromatin through antibodies, coupling target-proximal cleavage with adapter transfer. This page represents the original **bulk** method. Its universal-i5/barcoded-i7 PCR differs from the paper's separate single-cell method with 72 i5 and 72 i7 indices.

## Sources and version boundaries

- 🟢 [Kaya-Okur et al., Nature Communications (2019)](https://doi.org/10.1038/s41467-019-09982-5), Methods: transposome preparation, bulk CUT&Tag and sequencing. The original article specifies NEBNext HiFi PCR and AMPure DNA extraction before PCR.
- 🟢 Cited [Picelli et al., Genome Research (2014)](https://doi.org/10.1101/gr.177881.114), Methods, “Transposome assembly in solution”: the complete three loaded oligos below. Its supplemental inventory and reaction-condition table describe Picelli's own variants; those conditions are not CUT&Tag selections.
- 🟢 Historical author protocol [wnufdew, Version 1](https://doi.org/10.17504/protocols.io.wnufdew), March 6, 2019, integer ID 18868. The article links this version. Its PCR footnote permits Nextera primers **or** indexed primers from Buenrostro2015.
- 🟢 [z6hf9b6, Version 2](https://doi.org/10.17504/protocols.io.z6hf9b6), April 21, 2019, ID 22441, is a separate revision. The current default export is Version 3, February 25, 2020, bcuhiwt6, ID 33385. A redirect to that export does not make it the original cited method.
- 🟢 [Buenrostro et al., Nature (2015)](https://doi.org/10.1038/nature14590), original MOESM36 workbook: 92 indexed i5 and 96 indexed i7 PCR oligos. The two leading NoMX rows are short C1 PCR arms, not complete P5/P7 library primers.

## Exact loaded adapter pool

Picelli's source strands are written 5′→3′ with modifications retained:

ME-A: `{{= cut_tag.ME_A }}`

ME-B: `{{= cut_tag.ME_B }}`

ME reverse: `/5Phos/{{= cut_tag.ME_REVERSE }}`

Each long strand anneals through its mosaic end to the same phosphorylated reverse strand. The source specifies an equimolar preannealed A/B pool. That alone does not establish the adapter pairing of every loaded dimer or the frequency of each resulting fragment class.

## Molecular path and represented uncertainty

The graph declares antibodies, the source adapter pool and pA–Tn5 locally. Washing precedes magnesium activation. It then represents a conditional A/B fragment and explicitly shows inferred A/A and B/B possibilities. Their physical fates remain unresolved; the model does not turn the absence of two different handles into an absolute claim that single-primer amplification is impossible.

The published bulk route extracts DNA before gap filling and PCR. In the article's bulk Methods, 21 µL DNA, 2 µL of universal i5 plus uniquely barcoded i7 primer and 25 µL NEBNext HiFi 2× master mix are cycled 72 °C 5 min (gap filling), 98 °C 30 s, 14 cycles of 98 °C 10 s / 63 °C 30 s, then 72 °C 1 min; the PCR product is cleaned with 1.1× AMPure XP, two 80% ethanol washes and 30 µL Tris elution. The graph represents that post-PCR clean-up and its discarded supernatant; which molecules it removes is not characterized. **The published article and cited Version 1 describe different extraction routes:** the article's bulk Methods use AMPure; Version 1's DNA-extraction section uses phenol/chloroform and ethanol precipitation. The page follows the article. These are documented source/version differences; the sources do not establish which route every individual experiment used.

The exact universal i5 oligo is not assigned in this model. A conditional endpoint uses the explicitly permitted Buenrostro i7 family, with its selected member and sample assignment left unassigned. This is not evidence that every experiment used that family, a variable i5 barcode, or a particular eight-cycle dual-index run. Exact loaded sequencing primers are also unassigned. The article describes paired-end Illumina sequencing and an economical 25×25 option, rather than one universal acquired run.

The article also discusses pA–Tn5-associated *E. coli* DNA used for calibration. That source branch is outside the current molecular graph; no mechanism or fate for that carryover is inferred here.

## Source audit still in progress

Accessible Picelli originals have replaced the earlier retrieval gap. Historical protocols.io versions have been recovered separately. The review retains the precise scopes of both readers and pending linked-source checks. Unassigned primer choices and sample mapping are model questions, not accepted claims that publication-wide information is absent.
