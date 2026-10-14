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

## Oligo source review

🟢 The linked author preprint (doi:10.1101/219683, tag-design section) specifies
phosphorylated DPM, Odd and Even adaptors, a 3′ spacer on the DPM top strand
that prevents ligation to genomic DNA, and a terminal adaptor whose top strand
alone is 5′ phosphorylated. Its Table S5 contains the individual barcode
oligos, but that attachment is not in the open preprint PDF or the fetched
2018 author manuscript. The later author protocol (doi:10.1038/s41596-021-00633-y)
publishes a 96-member plate for each of DPM, Odd, Even and Terminal strand
types, plus 96 indexed amplification primers. Its first DPM top strand
`/5Phos/AAACACCCAAGATCGGAAGAGCGTCGTGTA/3SpC3/` matches the exact bases
and 3′-spacer role in the 2018 preprint example; the corresponding bottom strand
is `/5Phos/TGACTTGTCATGTCTTCCGATCTTGGGTGTTTT` (written 5′→3′).
The terminal top example `/5Phos/AGTTGTCACCATAATAAGATCGGAAGA` also matches.
Only these cross-matched examples are established for the 2018 design; the full
2022 plate must not silently replace the unavailable original Table S5.

🟢 The preprint prints the SPRITE library amplification primers:
Read 1 `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
and Read 2
`CAAGCAGAAGACGGCATACGAGATGCCTAGCCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.
These are library PCR primers, distinct from the shorter platform sequencing primers.
