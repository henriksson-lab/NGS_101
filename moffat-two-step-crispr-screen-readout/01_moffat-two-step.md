# Moffat two-step pooled CRISPR screen readout

The Moffat Lab TKOv3 workflow deliberately separates locus enrichment from indexed-library construction. PCR1 enriches the integrated guide region; PCR2 uses nested primers to add complete dual-index TruSeq arms.

## Source read

- 🟢 Moffat Lab, *Pooled CRISPR Screen Protocol v2*, revision 2019-07-23, section 3.5.3 and Tables 1–2, hosted by Addgene. The page follows the LCV2::TKOv3 branch.
- 🟢 The protocol specifies an approximately 600-bp PCR1 product, approximately 200-bp PCR2 product, and a run of 21 dark cycles + 26 Read 1 cycles + 8 cycles each for i7 and i5.

## Reactions (section 3.5.3, PDF text lines 625–863; bounded pass 2026-10-17)

- 🟢 PCR1: 2x NEBNext Ultra II Q5 Master Mix, 2.5 µL each 10 µM LCV2 primer, 3.5 µg gDNA per 50 µL;
  100 µg gDNA per sample (~200-fold TKOv3 coverage) in identical reactions; 98 °C 30 s, then 25 cycles of
  98 °C 10 s / 66 °C 30 s / 72 °C 15 s, 72 °C 2 min. 2 µL checked on a 1% gel (600 bp); all reactions of a
  sample are pooled.
- 🟢 PCR2: one 50 µL reaction per sample, 5 µL pooled PCR1, same master mix, unique i5/i7 combination;
  98 °C 30 s, then 10 cycles of 98 °C 10 s / 55 °C 30 s / 65 °C 15 s, 65 °C 5 min.
- 🟢 The whole 50 µL is run on a 2% gel, the 200-bp band is excised on a blue-light transilluminator and
  gel-extracted, then quantified by NanoDrop and Qubit. The graph models this as a cleanup whose gel
  remainder is discarded.
- 🟢 PCR1 primers are regular desalted oligos; PCR2 primers are desalted IDT Ultramers. No terminal
  modifications are printed.
- 🟡 The PCR1 product drawn on the page carries the 22-nt PCR2 forward annealing part
  (`TTGTGGAAAGGACGAAACACCG`, ending in the transcription-start G) and the 24-nt reverse annealing part,
  because the nested primers must anneal there; the remaining ~600-bp LCV2 sequence is not printed and is
  drawn as placeholders.
- 🟢 Section 3.5.4 (lines 1037–1066): HiSeq 2500 SR50 or NextSeq 500 SR75 with dual index. Without dark
  cycles, a 0–5-nt stagger is placed 5′ of the i5 annealing part and Read 1 is lengthened; that variant is
  not drawn. The plan records 26 acquired Read 1 cycles separately; their start after 21 dark additions remains unplaced in the ordinary primer-start geometry.

🟢 [Table 1 of the Addgene-hosted protocol](https://media.addgene.org/cms/filer_public/44/9f/449f22c6-f818-4af6-8af4-2bab1d20f00a/moffat_lab_crispr_screen_protocol_v2.pdf)
prints the LCV2 PCR1 pair and the full i5/i7 PCR2 primer families. The schematic uses
one D501/D701 index combination as a representative library. Table 2 prints a distinct
pLCKO2 branch, whose different PCR1 reverse primer and i5 family must not be read as
interchangeable with this page's LCV2 construct.

## Boundary

The same document mentions a possible one-PCR alternative. That alternative is not drawn here: the page is specifically the two-step workflow, while the Joung protocol has its own one-PCR page.

## Native migration scope (2026-10-18)

Actual native genomic mixtures feed PCR1 and the actual pooled enriched mixture feeds PCR2. Possible original DNA, RNA and primers remain in the carried mixture; amplification does not discard them. The graph includes the analytical aliquot, pool, PCR2 aliquot and remainder, and physical final gel fraction. Native flank extents and class-specific gel recovery remain unassigned. The nested-site panel represents the target component rather than complete mixture geometry.
