# Ribo-seq — original 2009 poly(A)-tail workflow

## Sources and version boundary

- 🟢 Ingolia et al., *Science* (2009), [doi:10.1126/science.1168978](https://doi.org/10.1126/science.1168978), Results and Figure 1: ribosome-protected RNA footprints are converted into a sequencing library by a strategy selected to avoid RNA-ligase bias.
- 🟢 Ingolia et al., *Nature Protocols* (2012), [doi:10.1038/nprot.2012.086](https://doi.org/10.1038/nprot.2012.086), “Linker ligation” (local JATS lines 98–104): the authors explicitly say their earlier method added a poly(A) tail to each RNA, then describe replacing that step with optimized ligation of a preadenylated linker. They also changed RT and PCR primer sequences so the later library could use standard Illumina primers and an optional six-nucleotide index.
- 🟢 The 2012 protocol prints `/5rApp/CTGTAGGCACCATCAAT/3ddC/`, a modified long RT primer, and indexed PCR primers (local JATS lines 188–208). These are **2012-version oligos**, not verified 2009 oligos.
- 🔴 The original 2009 supplementary methods, `NIHMS131482-supplement-MS_131482_Supplementary_Data.pdf`, are linked by the article but currently inaccessible through the local source route. Its exact oligos, junctions and run layout remain unverified. The source review checklist keeps this missing-information claim blocked.

## Molecular path represented here

The 2009 branch starts with nuclease protection, monosome/size selection, 3′ poly(A) tailing, reverse transcription and DNA library preparation. The graph shows the poly(A)-derived junction as a role, without inventing a tail length or oligo bases. Primer landing and final sequencing arms remain role-level pending the original supplement.

The 2012 article is a linked successor and a retrospective source for the tailing distinction. Its preadenylated linker, six-base sample index and standard-Illumina read sites are intentionally absent from the 2009 graph. This corrects an earlier repository model that placed the 2012 chemistry under the 2009 DOI; it is not a discrepancy within either publication.

## Bounded primary-source check (2026-10-16)

Inspected: Original PMC2746483 JATS24–44.

Original main article supports nuclease footprinting and avoidance of RNA ligase. Digested RNA is now an unresolved molecular product rather than presumed physical waste; tail length labelled unassigned. Preserved explicit 2009/2012 version boundary.

Deferred: Original 2009 SOM remains deferred/inaccessible; exact oligos, route details and run; no later linker route substituted. This is a capped check, not a full source audit or confirmation of missing information.
