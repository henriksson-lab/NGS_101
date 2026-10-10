# ChRO-seq — research notes

Defining source: Chu *et al.*, *Nature Genetics* (2018), [10.1038/s41588-018-0244-3](https://doi.org/10.1038/s41588-018-0244-3) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6204104/)).

## 1. What it is

Chromatin run-on sequencing (ChRO-seq) maps active RNA polymerase from isolated chromatin. Unlike nuclear run-on, the input can be solid tissue and can contain substantially degraded RNA.

## Sources read

- 🟢 The defining article's ChRO-seq/leChRO-seq Methods and supplementary-method description were read in the open manuscript.
- The repository source fetch records the publisher copy as manual: `python3 tools/get_sources.py --doi 10.1038/s41588-018-0244-3`.

## Construct-changing steps

1. 🟢 Chromatin is isolated in NUN buffer; engaged polymerases remain attached to template and nascent RNA.
2. 🟢 Run-on is performed in the presence of biotinylated ribonucleotides and sarkosyl. The newly extended 3′ end becomes the affinity handle.
3. 🟢 RNA is fragmented, biotin-labelled nascent RNA is captured on streptavidin, and strand-specific RNA adapters are installed.
4. 🟢 Reverse transcription and PCR complete the Illumina-compatible library. Repeated biotin selection suppresses unlabelled RNA background.

The biotin-NTP mark is not part of the sequenced cDNA: it selects the RNA whose endpoint is interpreted as polymerase position.

## Final library and read layout

- 🟢 The genomic alignment is strand-specific and the read end adjacent to the run-on-labelled RNA end is the positional measurement.
- 🟡 Adapter and sequencing-primer arms are represented by molecular role because the defining article does not prescribe one invariant primer sequence across ChRO-seq implementations.

## Boundary and uncertainty

leChRO-seq is a named low-input extension described in the same paper. This page shows ChRO-seq and does not silently merge the extra linear-amplification branch into it.
