# Safe-SeqS — research notes

**Defining paper:** Kinde *et al.*, “Detection and quantification of rare mutations with massively parallel sequencing,” *PNAS* (2011), DOI [10.1073/pnas.1105422108](https://doi.org/10.1073/pnas.1105422108).

## Sources read

- 🟢 Defining article and Methods, including Fig. 1 and the exogenous-UID workflow (PMCID PMC3111315).
- 🟢 The paper describes a random sequence such as 14 N placed 5′ of one gene-specific primer, limited-cycle UID assignment, universal amplification, and the ≥95% family criterion.

## Construct-changing steps

1. 🟢 One target-specific primer carries a 5′ universal handle followed by a random UID (the paper's example is 14 nt); the other target primer carries the opposite universal handle.
2. 🟢 Two UID-assignment PCR cycles tag template strands before bulk amplification.
3. 🟢 Universal PCR generates redundant descendants sharing the UID.
4. 🟢 Illumina sequencing groups reads by UID; a mutation is called a “supermutant” only when at least 95% of the UID family supports it.

## Boundary of the schematic

🟡 Safe-SeqS also describes endogenous identifiers from randomly sheared endpoints and a capture workflow. The public schematic represents the defining targeted exogenous-UID implementation, because it is the explicit primer-installed molecular-barcode construction. Universal arms are shown by role where the defining article does not make a single immutable oligo set part of the method.
