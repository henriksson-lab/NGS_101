# Microbe-seq — research notes

Defining source: Zheng et al., 2022, *High-throughput, single-microbe genomics with strain resolution, applied to a human gut microbiome*, Science, [doi:10.1126/science.abm1483](https://doi.org/10.1126/science.abm1483).

## Sources read

- 🔴 Published *Science* article and supplement: **not readable** in this audit. `get_sources.py` reports no open full text (manual row, publisher page). PMC idconv has no PMCID. OpenAlex and Unpaywall list one green open-access location, OSTI purl 1959222. `www.osti.gov` still times out from this host, but the Internet Archive snapshot of that purl (2025-02-02) returned the file on 2026-10-17: it is a 24-page PDF whose first page is the **bioRxiv v1 preprint** cover (`doi: 10.1101/2020.12.14.422699`, posted 14 December 2020), 3,131,102 bytes against the 3,131,106-byte copy already on disk. The indexed "accepted manuscript" is therefore the preprint, not the published text, and the OSTI route is closed rather than pending. A plain request for the publisher supplement (`science.org/doi/suppl/10.1126/science.abm1483/suppl_file/science.abm1483_sm.pdf`) answers HTTP 403, and Europe PMC reports no PMCID and no open full text for PMID 35653470. The published article and supplement must still be obtained by hand and compared with everything below.
- 🟡 Linked bioRxiv preprint 10.1101/2020.12.14.422699 v1 (superseded): supplementary Methods (`media-1`, text lines 40–113: MDA, tagmentation, bead synthesis, bead preparation, barcoding, pooling and library PCR, Illumina sequencing) and Table S6 (`media-3`, lines 217–700). The model uses these as **preprint-derived** facts because the published version is unavailable. They are never a published value.

## Molecular path (preprint)

🟡 Microbes are encapsulated with lysis reagents. Sample droplets are then merged with MDA-reagent droplets and incubated at 30 °C for 6–8 h, then 65 °C for 10 min. A second merger adds commercial Nextera tagmentation reagent (55 °C, 10 min). A third merger adds PCR reagents and one barcoding hydrogel bead. The barcoding program starts at 72 °C for 4 min, which covers gap filling of the tagmented ends. It continues at 98 °C for 30 s, then 10 cycles of 98 °C 7 s / 60 °C 30 s / 72 °C 40 s.

🟡 Droplets that merged during PCR are removed with a droplet size-filter device. The emulsion is broken with 20% PFO, and the aqueous phase is purified with 1.1× AMPure. A bulk Phusion PCR (5–10 cycles) uses `P5PE1` (IDT) and a Nextera i7 primer (Illumina) to add the sequencing adapters and **one** sample index. A 0.8× AMPure cleanup follows.

🟡 Beads are made by polymerising droplets that contain an acrydite-modified, photocleavable-linker hydrogel oligo. Two split-and-pool Bst 2.0 extension rounds follow, in 96 and then 384 wells, giving 96 × 384 = 36,864 combinations.

## Table S6 oligos (preprint, verbatim)

| Role | 5′→3′ oligo, spaces as printed between constant arms and barcode |
|---|---|
| Barcode-1 oligo-1 | `AAGGCGTCACAAGCAATCACTC AAACAAAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAG` |
| Barcode-2 oligo-1 | `TCTCTTATACACATCTGACGCTGCCGACGA AAACAAAC AAGGCGTCACAAGCAATCACTC` |

🟢 (in the preprint) Barcode 1 has 24 oligos each of length 8, 9, 10 and 11 nt; barcode 2 has 384 oligos of 8 nt. The self-test re-derives these counts from the table text.

🟡 Derived interlock (computed in `microbe_seq.py`): the barcode-1 3′ arm anneals to the hydrogel oligo, so the bead strand carries `{{= chemdraw.revcomp(microbe_seq.BC1_OLIGO_3P) }}`. The barcode-1 5′ arm is the 22-nt linker that barcode 2 anneals to. The barcode-2 5′ arm copies into `{{= microbe_seq.BEAD_TN5_ARM }}`, which is the Nextera S5 + mosaic-end adapter strand minus its last {{= len(microbe_seq.TN5_COPIED) }} bases. The bead primer therefore primes on the S5-adapter end of tagmented fragments.

## Final library and reads

{{= ' · '.join(s.name for s in microbe_seq.FINAL_LIBRARY) }}

🟡 Sequencing: a custom Read-1 primer (IDT) is required. A 100-bp run uses Read 1 = 45 bp of barcode only, Index 1 = 8 bp and Read 2 = the remainder (microbial). A 300-bp run uses Read 1 = 150 bp (first 45 barcode, last 75 microbial). Because the custom primer sequence is not printed, the page gives a source-declared read table instead of a computed primer landing.

🔴 Not printed in the inspected preprint sources: the initiating hydrogel oligo, `P5PE1`, the custom Read-1 primer, the commercial MDA primers, and any opposite (S7-side) primer in the barcoding droplet. These are open questions, not confirmed absences.

## Model corrections (2026-10-17 bounded pass)

- The earlier model drew a standard dual-index Nextera library: an i5 index, the standard Nextera Read 1 primer, a 16-nt barcode with no linker, and 150 + 150 reads. None of that is in the inspected source. It now follows the preprint: single Nextera i7 index, a bead barcode 1 (8–11 nt) + 22-nt linker + barcode 2 (8 nt) upstream of the Tn5 S5/ME arm, a custom unprinted Read-1 primer, and the 45/8/remaining read configuration.
- The model adds the droplet size filter, emulsion breaking and both AMPure cleanups. Preprint MDA, tagmentation and PCR programs are on the steps.

## Open questions

- Every preprint-derived item above, until it is compared with the published article and supplement. No open route remains: the one indexed green copy is the preprint, so this needs a hand-obtained *Science* article and supplement.
- Hydrogel oligo, P5PE1 and custom Read-1 primer sequences. Klein 2015 and Zilionis 2017 bead methods are cited for bead synthesis and were not followed in this capped pass.
