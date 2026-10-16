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
