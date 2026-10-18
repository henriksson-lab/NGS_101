# SMART-seq — original Tn5 branch

Ramskold et al. (2012), [doi:10.1038/nbt.2282](https://doi.org/10.1038/nbt.2282).

The selected page models the original Tn5 branch. The paper also prepares PE libraries by Covaris fragmentation, end repair, dA tailing, adapter ligation and PCR; its figure data use that route unless a legend says otherwise. That alternative route is not represented by a generic Nextera endpoint.

## Source facts and molecular scope

🟢 Methods print the CDS primer as `AAGCAGTGGTATCAACGCAGAGTACT(30)VN`. They name SMARTer II A and kit IS PCR reagents. Their native sequences remain unassigned in the current model; later SMART-seq2 oligos are not evidence of original kit identity.

🟢 The original Tn5 route uses approximately 5 ng amplified cDNA, 55 C tagmentation for 5 min, PB stripping and SPRI XP purification at sample:beads 1:1.6, followed by nine cycles of standard Nextera PCR. It is the modified Epicentre Nextera preparation, rather than Nextera XT.

🟢 Supplementary Table 1's opening caption explicitly says sequencing was single-end. Its Tn5 rows list 50-nt GAIIx or 40-nt MiSeq reads. The previous model overlooked that caption. Native primer sites, index members and acquired index reads remain unassigned.

🟡 Actual first-strand and template-switched cDNA are represented as single strands before amplification. Native tagmented substrates enter index PCR. Alternative end classes have unresolved recovery, because nonamplification does not establish physical removal.

## Oligos

The source CDS primer is given above; unassigned kit reagents remain named roles.

## Migration source pass (2026-10-18)

Read: PMC3467340 Methods; published Supplementary Figure 1; published Supplementary Table 1 caption and Tn5 rows. The bundle scraper was run and relevant source pointers checked. The contemporaneous Clontech manual, original Epicentre guide and alternative PE-route details remain deferred. This bounded source pass does not establish publication-wide absence or a complete review.

## Critical implementation follow-up (2026-10-18)

Retained the independent RNA template beside first-strand and switched cDNA. The PCR output includes an RNA-derived molecular population with unresolved thermal fate and later cleanup recovery; this conceptual class distinction does not assert physical isolation. Copied TSO orientation was rechecked against complement_segments: it already reverses order and complements bases while preserving UMI feature identity.
