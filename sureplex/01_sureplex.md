# SurePlex with VeriSeq PGS

## Defining source revisions

🟢 SurePlex Summary Protocol **15053626 v02, September 2020**; VeriSeq PGS Library Prep Reference Guide **15052877 v05, August 2021**; Guidance Note **1000000138837 v00**. Printed titles/revision identities were checked. The source bundle scraper returned zero oligos; that is not evidence of source-wide absence. Later adapter-sequence tables do not supply historical native strands or index widths here.

## Native amplification and tag preparation

The selected cell route performs source lysis/extraction, then the commercial **12-cycle preamplification** and separate **14-cycle amplification**. Their named buffers/enzymes and cycle programs are explicit. Buffers are declared as conditions; no physical primer pool is assigned to a buffer whose composition is unverified. Primer-pool species, first-copy architecture and native product ends remain unresolved, so these two stages use bounded commercial transitions. They do not assert an isothermal mechanism or a single symmetric PCR primer.

A parallel PBS no-cell control passes through the same two amplification stages as a separate reaction. Possible artifacts retain unresolved molecular identity/fate. The guidance note says they may yield unmapped sequencing reads which the system filters; that computational filter does not discard the material physically.

The library route selects **1 ng unpurified WGA product**, tagmentates with VeriSeq ATM/TD at 55°C for 5 min, neutralizes with NT at room temperature for 5 min, then performs the source initial **72°C 3-min extension** before 12-cycle PCR. Named S503/S504 and N701–N712 primer families remain native unknown strands; no eight-base physical/index-read default or completed Nextera shell is inserted.

Actual NT-treated/extended native fragment strands enter PCR. Other Tn5 end classes and original/off-target PCR molecules retain unresolved recovery. They are not discarded merely because they may fail a preferred amplification route.

## Physical normalization and sequencing

🟢 Source PCR cleanup uses 45 µl AMPure per 45 µl reaction, two 80% ethanol washes and RSB elution. LNA1/LNB1 normalization binding is followed by **two separate LNW1 washes**. NaOH elution yields source **single-stranded** library, then 25 µl supernatant is transferred into LNS1. Equal volumes of normalized libraries are pooled, diluted with HT1 and heat-treated before loading.

The final graph material stays single-stranded; no separate completed duplex panel is drawn. Normalization capture chemistry, selected strand identity, exact ends/index members/widths and class-specific carryover remain unassigned.

The Guidance Note describes a **36-cycle short-read recipe with single- or dual-index options**. This page selects the named dual-index preparation; acquired index-cycle counts and the exact native read mask are unassigned. No read-2 acquisition is implied.

## Bounded migration scope

This implementation uses the parent's three-original source handoff, with the same named source spans checked for identity, native amplification/NT/PCR/normalization and control/short-read qualifications. It adds no fourth original, recursive source hunt or missing-information acceptance. Single-index preparation, native primer chemistry, normalization-bead mechanism and run-specific details remain deferred; historical review status is preserved.
