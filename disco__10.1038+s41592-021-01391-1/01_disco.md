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

🟡 The page now carries these four oligos verbatim and computes the SeqB
landing on the printed P5 SMART PCR end (`chem._decompose()` checks that each
printed oligo decomposes into the canonical SMART handle, P5, the GCCTGTCCGCGG
spacer, GAAT and AC). The selected i7 index and the lot-specific bead barcode
remain variable.

## Library workflow in the author protocol (read 2026-10-17)

Source: `DisCo_protocol.pdf` (EPFL Deplancke lab, PDF created 14 Sep 2021,
15 pages; local `.txt` lines 94-122, 389-406, 478-694). It calls itself the
authors' "most up-to-date" protocol, so its equivalence with the published
experiments is not established.

1. 🟢 Beads (40-80 µL of stock, 40 µm strained) are resuspended in lysis buffer
   (Optiprep, Sarkosyl, EDTA, Tris, DTT); one cell and one bead are co-encapsulated
   on the DisCo chip and droplets collected under stuffer droplets.
2. 🟢 The collection tip is moved to a bead-capture (cp) chip; beads retained on
   pillars are washed twice with 6X SSC and once with RT buffer.
3. 🟢 RT on chip: Maxima H− RTase, Maxima RT buffer, Ficoll PM-400, dNTPs, RNase
   inhibitor and the Template Switch Oligo (10 µL of 50 µM per 200 µL mix),
   42 °C for 90 min; then TE-SDS once, TE-TW twice.
4. 🟢 10 mM Tris wash, Exonuclease I 37 °C 45 min on chip, TE-SDS, 2x TE-TW,
   water; beads eluted with water, pelleted and moved to a PCR tube.
5. 🟢 SMART PCR on beads with 2x KAPA HiFi HotStart ReadyMix: 4 cycles
   (65 °C anneal) then 15-19 cycles (67 °C); supernatant moved away from the beads;
   0.6x AMPure ("to remove fragments smaller than 200 bp"); Qubit.
6. 🟢 Tagmentation: 2 µL cDNA at 1 ng/µL, in-house Tn5 "loaded with p7|ME"
   (Picelli et al. 2014), 5x TAPS, 55 °C 6 min; stop with 0.2% SDS.
7. 🟢 Library PCR: i7 index primer (Buenrostro et al. 2013 set) + P5 SMART PCR
   with KAPA HiFi; 72 °C 3 min, 98 °C 30 s, 16 cycles (98 °C 10 s, 63 °C 30 s,
   72 °C 1 min), 72 °C 5 min.
8. 🟢 0.7x AMPure XP; Bioanalyzer/Fragment Analyzer.
9. 🟢 NextSeq: Read 1 = 20 bp with the custom Drop-seq SeqB primer, Read 2 = 50 bp.

🟡 Inferred, not printed: the non-templated C / rGrGrG pairing intermediate;
that Exonuclease I removes unextended single-stranded bead oligos (the protocol
gives no purpose); that only bead-end fragments (handle followed by AC) are
selected by the P5 SMART PCR and SeqB 3′ end.

🔴 Not in the three inspected originals (not confirmed absent): the bead oligo
sequence (barcode/UMI lengths, constant bases), Tn5 p7|ME adapter strands, i7
index-primer bases, Read 2 and index-read primers and index-read length.

## Bounded pass 2026-10-17

Retrieved: author protocol via the existing `extra_sources.tsv` row (now on
disk); defining Nature Methods full text still subscription-gated (nature.com
page shows only abstract/legends; not in PMC or Europe PMC). As labelled
fallback the superseded bioRxiv preprint 10.1101/2020.05.19.103812 v3 and its
media-1 supplement were fetched into `sources/doi_10.1101+2020.05.19.103812/`.

Preprint use: lead only, not adopted into the model (published version first;
its applicability to the published experiments cannot be verified). Its
supplement Materials list (media-1 `.txt` lines 320-360) names beads
"Chemgenes, lot 051917, Macosko-2011-10", prints the same TSO, SMART PCR and
P5 SMART PCR sequences, and cites "custom Nextera oligos" for indexing. It
differs from the later author protocol in cycle numbers and cleanup ratios
(cDNA PCR 18-23 cycles vs 4 + 15-19; library PCR 15 vs 16 cycles; library
cleanup CleanPCR 0.6x vs AMPure XP 0.7x). These are different document
versions, not a published discrepancy.

Model corrections: the generic shared RT-then-single-PCR route was replaced
by the protocol's own graph (cp-chip capture, RT + template switching, washes,
Exonuclease I with removal by washing, bead elution, SMART PCR, bead removal,
0.6x cleanup, p7|ME Tn5 tagmentation, SDS stop, P5 SMART / i7 PCR, 0.7x
cleanup); the false "no universal binding oligo" sequencing statement was
replaced by the printed SeqB Read 1 landing; the 8-nt shared sample-index role
and legacy role library were dropped.

## Bounded curation round 3 (2026-10-16)

Inspected: Published MOESM1 PDF title and SupplementaryFigure1 caption; MOESM2 PDF extracted text; MOESM7 XLSX SupplementaryTable1 header/rows and Table2 header.

No correction to molecular claims: available publisher supplements inspected here concern live-cell gating and biological/sample results; they do not support completing the molecular workflow. Defining full text remains unavailable in the local manifest.

Deferred: Defining article and linked detailed protocol remain unavailable/unread; bead primer, RT/cleanup/library stages and run assignment not independently verified; no missing-information claim accepted. This capped pass does not establish a complete source review or publication-wide absence.
