# ChRO-seq — research notes

Defining source: Chu *et al.*, *Nature Genetics* (2018), [10.1038/s41588-018-0244-3](https://doi.org/10.1038/s41588-018-0244-3) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6204104/)).

## 1. What it is

Chromatin run-on sequencing (ChRO-seq) maps active RNA polymerase from isolated chromatin. Unlike nuclear run-on, the input can be solid tissue and can contain substantially degraded RNA.

## Sources read

- 🟢 The defining article's ChRO-seq/leChRO-seq Methods and supplementary-method description were read in the open manuscript.
- 🟢 Its Methods cite Mahat *et al.*, *Nature Protocols* (2016), DOI `10.1038/nprot.2016.086`, as the library-preparation method. Table 1 of that open protocol (`nprot.2016.086_PMC5502525.xml.txt`, lines 2870–2905) was checked for the prior PRO-seq oligos. The citation establishes procedural ancestry, not that every ChRO-seq sample used an identical synthesis order.
- The repository source fetch records the publisher copy as manual: `python3 tools/get_sources.py --doi 10.1038/s41588-018-0244-3`.

## Construct-changing steps

1. 🟢 Chromatin is isolated in NUN buffer; engaged polymerases remain attached to template and nascent RNA.
2. 🟢 Run-on is performed in the presence of biotinylated ribonucleotides and sarkosyl. The newly extended 3′ end becomes the affinity handle.
3. 🟢 RNA is fragmented, biotin-labelled nascent RNA is captured on streptavidin, and strand-specific RNA adapters are installed.
4. 🟢 Reverse transcription and PCR complete the Illumina-compatible library. Repeated biotin selection suppresses unlabelled RNA background.

The biotin-NTP mark is not part of the sequenced cDNA: it selects the RNA whose endpoint is interpreted as polymerase position.

## Final library and read layout

- 🟢 The genomic alignment is strand-specific and the read end adjacent to the run-on-labelled RNA end is the positional measurement.
- 🟢 The cited PRO-seq method prints 3′ RNA adapter VRA3 `5′-pGAUCGUCGGACUGUAGAACUCUGAAC-/Inverted dT/-3′`, 5′ RNA adapter VRA5 `5′-CCUUGGCACCCGAGAAUUCCA-3′` (unphosphorylated 5′ end), RT/PCR primer RP1 `5′-AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA-3′`, and indexed PCR primer RPI-n `5′-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA-3′` (Table 1; `NNNNNN` is the sample barcode).
- 🟡 ChRO-seq explicitly follows that prior workflow for 3′ and 5′ adapter ligation, reverse transcription and amplification, so these are strong candidate oligos for implementations that did not substitute adapters. The defining paper does not name VRA3/VRA5/RP1/RPI-n or document which samples carry six-base UMIs; their assignment to every published ChRO-seq library remains an inference. The schematic retains role-level arms until that version boundary is resolved.

The public page lists those four oligos as the **cited predecessor's** sequences. It does not label them as confirmed ChRO reagent orders or assign their bases to the ChRO-specific final construct.

## Boundary and uncertainty

leChRO-seq is a named low-input extension described in the same paper. This page shows ChRO-seq and does not silently merge the extra linear-amplification branch into it.
