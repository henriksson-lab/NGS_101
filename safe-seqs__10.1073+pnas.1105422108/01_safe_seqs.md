# Safe-SeqS — research notes

**Defining paper:** Kinde *et al.*, “Detection and quantification of rare mutations with massively parallel sequencing,” *PNAS* (2011), DOI [10.1073/pnas.1105422108](https://doi.org/10.1073/pnas.1105422108).

## Sources read

- 🟢 Defining article and Methods, including Fig. 1 and the exogenous-UID workflow (PMCID PMC3111315).
- 🟢 The paper describes a random sequence such as 14 N placed 5′ of one gene-specific primer, limited-cycle UID assignment, universal amplification, and the ≥95% family criterion.
- 🟡 A US Patent Trial and Appeal Board exhibit reproduces the original PNAS article and ten-page supplement, including Table S4. The supplement's methods were readable, but the table's sequence cells were not captured by text extraction and the PDF could not be fetched locally for visual transcription. The publisher and PMC copies remain gated. This is an open source-reading task, not evidence that the oligos were unpublished.
- 🔴 The exact ordered oligos and modifications in Table S4 have not yet passed two independent source readings. The cited Herman *et al.* capture predecessor also needs full inspection.

## Construct-changing steps

1. 🟢 One target-specific primer carries a 5′ universal handle followed by a random UID (the paper's example is 14 nt); the other target primer carries the opposite universal handle.
2. 🟢 Two UID-assignment PCR cycles tag template strands before bulk amplification.
3. 🟢 The original supplemental method digests excess assignment primers with Exonuclease I before universal amplification; the resulting primer fragments are discarded.
4. 🟢 Universal PCR generates redundant descendants sharing the UID.
5. 🟢 Illumina sequencing groups reads by UID; a mutation is called a “supermutant” only when at least 95% of the UID family supports it.

## Boundary of the schematic

🟡 Safe-SeqS also describes endogenous identifiers from randomly sheared endpoints and a capture workflow. The public schematic represents the targeted exogenous-UID implementation. It uses a representative 14-nt UID and an inferred TruSeq-compatible PCR shell; Table S4 lists actual ordered oligos that still need transcription and reconciliation. The depicted read orientation and primer landing sites are therefore conditional, while target-specific arms are shown by role.
