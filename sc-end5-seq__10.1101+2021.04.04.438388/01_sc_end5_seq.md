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
These published sequences are absent from the current role-only molecular model.

🟢 The nearby template-switch oligo
5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′ is explicitly described for
the **sc-end3-dT control** made with the Chromium 3′ Reagent Kit v2. The two
sc-end5 variants use the initial Chromium 5′ Library and Gel Bead Kit; the
preprint does not give the exact oligo carried by those 5′ gel beads. The
3′-kit oligo must not be substituted for the 5′-kit bead oligo.

🟡 The source names a kit cDNA primer mix and i7 sample-index plate but
does not state which index or the full 5′-kit barcode/UMI oligo sequences.
Exact read-primer placement needs a version-matched 10x kit specification
and a chosen index, in addition to transcribing the source-published RT primers.
