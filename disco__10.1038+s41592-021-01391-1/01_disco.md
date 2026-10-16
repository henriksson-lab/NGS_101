# DisCo

**Defining source:** doi:10.1038/s41592-021-01391-1.

## Molecular boundary

🟢 Machine vision and multilayer valves position one low-input cell and one capture bead, co-eject them into a droplet and route valid pairs to the sample outlet.

## Evidence boundary

🟢 The library chemistry is Drop-seq-derived; deterministic co-encapsulation is the defining innovation.

## Author protocol: published library oligos

🟢 The authors' [DisCo step-by-step laboratory protocol](https://www.epfl.ch/labs/deplanckelab/wp-content/uploads/2021/09/DisCo_protocol.pdf)
(p. 2, “Required primers”; p. 15, “Sequencing Strategy”) publishes the
following ordering-level oligos for the DisCo implementation:

| Role | 5′→3′ sequence as published |
|---|---|
| Template switch oligo | AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG |
| SMART PCR primer | AAGCAGTGGTATCAACGCAGAGT |
| P5 SMART PCR | AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C |
| Custom Drop-seq SeqB sequencing primer | GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC |

🟢 Here `rGrGrG` and the asterisks on the P5 SMART PCR oligo are kept
as written; they must not be silently converted to plain DNA. The author
protocol directs Read 1 = 20 cycles and Read 2 = 50 cycles on NextSeq and
names a Buenrostro et al. i7 index-primer set without fixing a sample index.

🟡 The role-only page currently omits these published bases and therefore
cannot place its known custom SeqB primer. The selected i7 index and the
lot-specific Drop-seq bead barcode remain variable. The defining Nature
Methods article is subscription-gated through the source fetcher, but the
author-maintained step-by-step protocol resolves the primer identities above.

## Bounded curation round 3 (2026-10-16)

Inspected: Published MOESM1 PDF title and SupplementaryFigure1 caption; MOESM2 PDF extracted text; MOESM7 XLSX SupplementaryTable1 header/rows and Table2 header.

No correction to molecular claims: available publisher supplements inspected here concern live-cell gating and biological/sample results; they do not support completing the molecular workflow. Defining full text remains unavailable in the local manifest.

Deferred: Defining article and linked detailed protocol remain unavailable/unread; bead primer, RT/cleanup/library stages and run assignment not independently verified; no missing-information claim accepted. This capped pass does not establish a complete source review or publication-wide absence.
