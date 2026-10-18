# NEBNext Single Cell / Low Input RNA

The page selects **Section 1, cells**, from the original **E6420 manual v6.0_01/24**, and its NEBNext adapter route requiring USER. Printed kit and revision identity were checked. The local manifest lacks a populated retrieval row; this reading does not upgrade that historical provenance mapping.

## Source oligos

🟢 RT primer: `AAGCAGTGGTATCAACGCAGAGTAC(T)30V`.

🟢 TSO: `GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG`; its terminal guanosines are RNA.

🟢 One cDNA PCR primer species: `AAGCAGTGGTATCAACGCAGAGT`.

## Actual molecular route

First-strand and switched cDNA are single strands beside independent source RNA and possible cellular DNA. The switched strand enters single-primer PCR; a completed duplex is not substituted for that precursor. The displayed copied TSO boundary and amplicon extent are explicitly inferred model references. RNA, original cDNA and off-target material retain unresolved physical fates.

The first cDNA recovery uses 0.6X beads and elutes 50 µl. Then 45 µl Bead Reconstitution Buffer is added to the eluted cDNA plus the same beads for a second binding/wash/recovery (§1.6.11); it is a separate physical operation.

Ultra II FS fragmentation, repair and 3′ dA addition precede native adapter ligation. The selected NEBNext adapter is processed with 3 µl USER at 37°C for 15 min. Its native structure/cleavage geometry is unassigned; USER changes the adapter state and produces unassigned molecular fragments, rather than an unchanged final indexed library. The manual's indexed-UMI adapter route does not require USER and is separate/not recommended for this workflow.

A 0.8X recovery precedes index PCR on actual USER-processed native fragments. Final PCR recovery uses 0.9X beads. The separately purchased oligo kit can supply single- or dual-index arrangements; no generic dual-eight-base shell is assigned.

## Sequencing and deferred work

The manual's typical paired 75-base read length is a recommendation. Actual run/index acquisition, chosen adapter/primer kit, USER-sensitive strand structure, first-copy reaction proof and class-specific RNA/DNA recovery remain unassigned. Source trimming motifs are partial read-through sequences and do not establish complete library ends.

## Bounded migration pass

Read one original, E6420 v6.0_01/24; mandatory bundle scraper hits checked against its oligo table. No newer oligo manual substitutes for the chosen historical kit. Section 2 purified-RNA handling and native unresolved details remain deferred; no absence claim or source-status upgrade is accepted.
