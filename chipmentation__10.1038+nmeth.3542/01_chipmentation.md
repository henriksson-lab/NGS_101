# ChIPmentation — research notes

## 1. What it is

ChIPmentation replaces conventional post-ChIP end repair and adapter ligation with one
Tn5 tagmentation reaction performed on bead-bound immunoprecipitated chromatin.

## 2. Sources read

- 🟢 Schmidl et al. 2015, DOI
  [10.1038/nmeth.3542](https://doi.org/10.1038/nmeth.3542), defining paper and Methods.

## 3. Molecular order

🟢 Crosslinked chromatin is sonicated and immunoprecipitated. While the complexes remain
on beads, Nextera Tn5 fragments and adapter-tags the selected DNA. Proteinase treatment
releases DNA and reverses crosslinks; indexed PCR completes the flow-cell arms.

## 4. Relationship to neighboring protocols

ChIPmentation selects target-associated chromatin before free Tn5 acts. CUT&Tag instead
tethers Tn5 to the antibody inside permeabilized cells. That difference in targeting
geometry is why these remain separate protocol pages despite sharing a Nextera endpoint.
