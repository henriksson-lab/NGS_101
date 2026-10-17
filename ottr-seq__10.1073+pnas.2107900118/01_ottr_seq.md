# OTTR-seq

**Defining source:** Upton *et al.*, *PNAS* (2021), doi:10.1073/pnas.2107900118.
The published PMC article defines the active workflow. Earlier notes copied oligos
from bioRxiv v1; those sequences have been withdrawn pending verification against
published SI Table S1. The locally downloaded supplementary ZIP contains article
images, not the SI Appendix. This is an unresolved retrieval/checking obligation,
not a claim that the published oligos are absent.

## Evidence used

- 🟢 Published Methods specify ddA/ddG tailing, shrimp alkaline phosphatase,
  the primer duplex with +1T/+1C overhangs, terminating template and the two
  library branches.
- 🟢 Universal adapters use Q5 PCR (4–8 cycles); full-length adapters use
  RNase A/H treatment and omit PCR.
- 🟡 The graph selects the universal-adapter branch. Exact adapter and index
  sequences remain unassigned; the published MiniSeq 75-cycle kit does not by
  itself establish the split of acquired reads and indexes.

## Construct logic

The nonextendable ddR both normalizes heterogeneous input ends and prevents RNA–RNA
cross-priming. Adapter order is chemical: the initial duplex can accept ddR-tailed input,
whereas the 3′C template is recruited only after the intermediate cDNA receives dG.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC8594491 JATS192–196; local PMC supplementary ZIP member inventory (article image files, no SI Appendix).

Withdrew active oligos transcribed solely from superseded bioRxiv; published SI TableS1 remains unverified. Added shrimp alkaline phosphatase treatment between ddNTP tailing and relay; specified Q5 4–8 cycles for universal adapters and kept PCR-free full-length branch separate.

Deferred: Published SI TableS1 retrieval/verification, exact template/primer modifications, RNA hydrolysis/recovery and acquired read mask; no absence claim accepted. This pass does not establish a complete source review or publication-wide absence.
