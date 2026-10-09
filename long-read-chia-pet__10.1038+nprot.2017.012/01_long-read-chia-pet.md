# Long-read ChIA-PET

## 1. What it is

Long-read ChIA-PET enriches chromatin contacts associated with a chosen protein, joins paired genomic fragments through a bridge linker, and converts the products into longer paired-end tags.

## 2. Authoritative source

🟢 Li et al., *Long-read ChIA-PET for base-pair-resolution mapping of haplotype-specific chromatin interactions*, Nature Protocols (2017), DOI `10.1038/nprot.2017.012`.

The PMC full text, including the oligo list and stepwise procedure, was read.

## 3. Construct-changing operations

1. 🟢 Dual-cross-linked, sonicated chromatin is immunoprecipitated for the target protein.
2. 🟢 Fragment ends are repaired and dA-tailed.
3. 🟢 Bridge linker-F is `/5Phos/CGCGATATC/iBIOdT/TATCTGACT`; bridge linker-R is
   `/5Phos/GTCAGATAAGATATCGCGT`. Their 18-bp core leaves a 3′ T at each end.
4. 🟢 One linker bridges two dA-tailed chromatin fragments in a single ligation reaction.
5. 🟢 Nextera Tn5 tagmentation runs for 5 min at 55 °C.
6. 🟢 Streptavidin captures tagmented fragments retaining the internal biotin linker;
   on-bead indexed PCR uses no more than 13 cycles.

## 4. Final library and reads

The final library has canonical Nextera sequencing ends around `genomic tag — bridge —
genomic tag`. Paired 150- or 250-base reads approach the bridge from the outer Tn5 ends.
The paper accepts either bridge orientation in a read, trims it, and maps the two flanking
genomic tags as the contact pair.
