# sciPlex-ATAC-seq

**Defining preprint:** Booth et al. (2023), doi:10.1101/2023.03.05.531201; published in *Genome Research* as “High-capacity sample multiplexing for single cell chromatin accessibility profiling.” Open manuscripts were read from PMCID PMC10028990 and PMC10696879.

## Molecular boundary

🟢 Permeabilized nuclei absorb unmodified, sample-specific single-stranded DNA hash oligos. An indexed primer-extension step makes polyadenylated hashes recoverable in the ATAC workflow.

🟢 Two-level combinatorial indexing gives each nucleus an ordered cell identity independent of its sample hash. Hash reads and ATAC reads are sequenced separately in the two-level implementation; the later three-level implementation handles them differently.

## Evidence boundary

🟢 The hash is modelled as a sample index, never as a cell barcode. The two combinatorial well indices jointly identify a nucleus.
