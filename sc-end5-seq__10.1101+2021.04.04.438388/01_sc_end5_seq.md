# sc-end5-seq

**Defining source:** doi:10.1101/2021.04.04.438388.

## Molecular boundary

🟢 Oligo(dT)-primed and random-hexamer-primed variants use 10x-compatible template switching to preserve transcript 5′ ends with cell barcodes and UMIs.

## Evidence boundary

🟢 The defining source is a preprint. Both published priming variants are shown and neither is silently chosen as canonical.

## Primer evidence

🟢 The preprint's library-preparation methods (extracted PDF lines 1409–1457)
give the actual 5′-assay RT primers. `sc-end5-dT` uses
5′-AAGCAGTGGTATCAACGCAGAGTAC-T(30)-VN-3′; `sc-end5-rand` uses
5′-AAGCAGTGGTATCAACGCAGAGTACNNNNNN-3′. The run used 75 + 75 bases.
Both published sequences now appear as exact primer inputs in the reaction graph and
in the page's oligo list; the poly(dT) shorthand expands to 30 thymidines in the model.

🟢 The nearby template-switch oligo
5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′ is explicitly described for
the **sc-end3-dT control** made with the Chromium 3′ Reagent Kit v2. The two
sc-end5 variants use the initial Chromium 5′ Library and Gel Bead Kit; the
preprint does not give the exact oligo carried by those 5′ gel beads. The
3′-kit oligo must not be substituted for the 5′-kit bead oligo.

🟢 Official 10x CG000086 Rev K (2019) covers the same 5′ Library & Gel Bead Kit
product numbers PN-1000006/1000014 and prints the Read 1–16-base barcode–10-base
UMI–TSO bead oligo, RT product, Read 2 adapter and i7 PCR family in its appendix.
Two independent readings of that later original are recorded in the source checklist.
It is a same-kit architecture comparator, not proof of the exact guide revision
used for the sc-end5 experiments.

🟡 The preprint names a kit cDNA primer mix and i7 sample-index plate but does
not state the exact kit revision or which i7 well each library used. Those
run-specific details remain blocked in the source review. The preprint reports
75-base paired-end sequencing; a version-matched guide and run assignment would
be needed to place the original experiment's exact primer and index reads.

## Reaction order and run details (2026-10-17 bounded pass)

🟢 Source identity: bioRxiv v2 (posted 21 April 2021), Moody, Kouno *et al.*,
"Profiling of transcribed cis-regulatory elements in single cells". No published
version of this preprint was identified in this pass; the later tCRE-atlas paper by
the same group is a separate study and was not used.

🟢 Methods (PDF text lines 1380–1440): sc-end5 cells were loaded on Chip A with the
initial Chromium Single Cell 5′ Library and Gel Bead Kit; the random primer was used
at 17 µM; RT ran in droplets at 53 °C for 45 min; cDNA was amplified with the kit
cDNA primer mix for 12 (sc-end5-dT) or
14 (sc-end5-rand) cycles, "followed by the
standard steps according to manufacturer's instructions"; libraries took different
i7 sample-index plate indexes and were sequenced on a HiSeq 2500 as 75-base paired
reads (PBMC sc-end5-dT likewise, lines 1440–1460). The model now draws the cDNA PCR
as its own step and the unspecified kit steps as a warned unresolved boundary,
followed by the index PCR. Splitting the cDNA primer mix into one primer per cDNA end
is 🟡 inferred.

🟢 Methods "Removing strand invasion artifacts" (lines ~1790–1815) aligns genomic
sequence to the TS-oligo sequence `TTTCTTATATGGG`, and the
unencoded-G procedure (lines ~1740–1785) looks for its last 5 nt (ATGGG) before the
soft-clipped junction. This is the only sc-end5 bead-oligo sequence the preprint
prints; it is written as DNA read sequence, so the bead oligo's 5′ portion and
terminal ribonucleotides are not given there. It matches the 3′ end of the
CG000086 Rev K bead TSO comparator (…TTTCTTATATrGrGrG). Strand-invader removal is
a computational filter, not a physical separation.

🔴 Still open: exact kit guide revision, cDNA primer bases, fragmentation/adapter
steps, chosen i7 wells, and acquired index-read cycles.
