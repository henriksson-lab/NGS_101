# OTTR-seq

**Defining source:** Upton *et al.*, *PNAS* (2021), doi:10.1073/pnas.2107900118.
The PMC article and the authors' bioRxiv v1 preprint were read. The preprint appends
Table S1; the final PNAS SI Appendix could not be retrieved, so its agreement with the
preprint table remains unverified.

## Evidence used

- 🟢 The paper defines the Mn²⁺ ddA/ddG tail, phosphatase cleanup, +1Y first template
  jump, dG-assisted second jump onto the 3′C adapter template, and PCR versus full-length branches.
- 🟢 BioRxiv v1 Table S1 gives the primer-duplex and adapter-template sequences and blocks.
- 🟡 The final construct shows the universal-adapter branch; indexes are variable features.

## Oligos in the universal-adapter branch

The sequences below are printed in the authors' bioRxiv v1 Table S1 (PDF text lines
768–1200). The final PNAS article cites an SI Table S1 that is currently inaccessible;
whether it retained all the preprint variants is unresolved. Modifications and RNA
bases are retained as written; `N` positions are variable index bases.

| Role | Ordered oligo, 5′→3′ | Evidence |
|---|---|---|
| +1T DNA primer | `/5Cy5/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT` | 🟢 Table S1, universal duplex v2 |
| blocked duplex RNA | `rGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrU/3SpC3/` | 🟢 Table S1, RNA complement v1 |
| terminating adapter template | `/5Cy5/ACACTCTTTCCCTACACGACGCTCTTCCrGrArUrCrUrC` | 🟢 Table S1 |
| P7–i7 PCR primer | `CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTG` | 🟢 Table S1 |
| P5–i5 PCR primer | `AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGAC` | 🟢 Table S1 |

🟡 The P5 and P7 PCR oligos are intentionally partial: their 3′ ends land in
the adapter-derived cDNA and polymerase fills in the rest of the standard
Illumina read-primer regions. The final article states that the resulting libraries
were run on a MiniSeq 75-cycle kit. The inspected article and preprint do not name a separate *custom*
sequencing primer. The current role-only final-library drawing therefore
cannot establish an exact primer landing coordinate, even though the
adapter and PCR oligos are available.

## Construct logic

The nonextendable ddR both normalizes heterogeneous input ends and prevents RNA–RNA
cross-priming. Adapter order is chemical: the initial duplex can accept ddR-tailed input,
whereas the 3′C template is recruited only after the intermediate cDNA receives dG.
