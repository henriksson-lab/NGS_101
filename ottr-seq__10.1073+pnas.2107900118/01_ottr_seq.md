# OTTR-seq

**Defining source:** Upton *et al.*, *PNAS* (2021), doi:10.1073/pnas.2107900118.
The PMC article and its supplementary Table S1 were read; the linked preprint supplied a
text-extractable copy of the same oligo table.

## Evidence used

- 🟢 The paper defines the Mn²⁺ ddA/ddG tail, phosphatase cleanup, +1Y first template
  jump, dG-assisted second jump onto the 3′C adapter template, and PCR versus full-length branches.
- 🟢 Supplementary Table S1 gives the primer-duplex and adapter-template sequences and blocks.
- 🟡 The final construct shows the universal-adapter branch; indexes are variable features.

## Construct logic

The nonextendable ddR both normalizes heterogeneous input ends and prevents RNA–RNA
cross-priming. Adapter order is chemical: the initial duplex can accept ddR-tailed input,
whereas the 3′C template is recruited only after the intermediate cDNA receives dG.
