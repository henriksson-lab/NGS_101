# ChAR-seq — research notes

## Sources read

- 🟢 Bell et al. 2018, DOI `10.7554/eLife.27024`, eLife paper and appended detailed protocol.
- 🟢 Bell et al. 2020, PMCID `PMC7607860`, detailed Current Protocols implementation with the two bridge oligos printed with modifications.

## Construct path

🟢 Crosslinked chromatin-associated RNA is fragmented and ligated in situ to the RNA-compatible end of a biotinylated bridge. The bridge’s opposite end is joined to nearby genomic DNA. RNA is converted to cDNA, material is fragmented and bridge-bearing junctions are affinity enriched before single-end library construction. The internal bridge provides the orientation used to distinguish RNA-derived from DNA-derived sequence.

## Bridge geometry

🟢 The bridge top is `/5rApp/AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG/3SpC3/`. Its activated 5′ single-stranded end ligates specifically to fragmented RNA 3′-OH ends, `NNN` is the molecular identifier, and the 3′ C3 spacer prevents an undigested bridge from entering the DNA ligation.

🟢 The bottom strand is `/5Phos/GATCTGCGACTTAATTAAAGATCCTTGGACGCCGG/iBiodT/T`. The two strands contain a 31-base exact antiparallel core; their remaining ends encode the RNA-ligation tail, internal biotin and DNA-facing DpnII arrangement rather than a generic double-stranded linker.

🟢 DpnII digestion removes the blocked distal bridge piece and cuts nearby genomic DNA. T4 DNA ligase then joins compatible DpnII cohesive ends. RNase H and DNA polymerase I complete the cDNA-containing duplex, and sonication followed by streptavidin capture enriches only molecules retaining the bridge.

🟢 The defining paper used NEBNext end repair/dA tailing, NEBNext hairpin adapters, USER opening and TruSeq-compatible indexing, then single-end 152-base sequencing. 🟡 The named commercial workflow determines the sequencing architecture, but the exact proprietary hairpin oligo bases are not printed in the sources read and remain visibly inferred on the schematic.
