# SPRITE — research notes

## Sources read

- 🟢 Quinodoz et al. 2018, DOI [10.1016/j.cell.2018.05.024](https://doi.org/10.1016/j.cell.2018.05.024), PMC6548320, main text, Methods and tag-design section.

## Molecular order

🟢 DSG/formaldehyde-crosslinked nuclear complexes are fragmented, sparsely coupled to
NHS beads, washed, blunt-ended, phosphorylated and dA-tailed. Beads undergo five
split–ligate–quench–pool rounds: DPM, Odd, Even, Odd and Terminal. Crosslinks are then
reversed and the DPM/Terminal arms support Illumina library PCR. Reads sharing the entire
tag chain are assigned to one multiway SPRITE cluster.

## Identifier semantics

🟢 DPM and Terminal each contain a 9-nt well identifier; Odd and Even tags each contain
17 nt. The five pieces jointly identify a spatial complex. The model records each part as
a structured combinatorial identifier and also validates a composite cluster identifier
over the five segments. Genomic molecules in the cluster are never proximity-ligated to
one another.
