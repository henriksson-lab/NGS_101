# UDA-seq — research notes

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published in the sources inspected. The construction model is `tools/chem.py`.

## Sources read

- 🟢 Wang et al. 2025, *Nature Methods*, DOI `10.1038/s41592-024-02586-y`, PMC12165859
  (PMC Cloud JATS `s41592-024-02586-y_PMC12165859.xml.txt`): Results "Theory and design of
  UDA-seq" (lines 149–155); Methods lines 321–399 (nuclei isolation, fixation, i7-only
  Tn5 transposome assembly, sc3′/sc5′ RNA UDA-seq, Multiome UDA-seq and its library
  construction); preprocessing paragraph on round-1/round-2 barcodes (lines 433, 445).
- 🟢 Supplementary Table 7 (`MOESM8`, sheets "primer summary", "Truseq-i5 index primer",
  "Nextera-i7 index primer").
- Not read (deferred): protocols.io refs 62 (sc RNA UDA-seq, doi:10.17504/protocols.io.n2bvjne15gk5/v1)
  and 63 (Multiome UDA-seq, doi:10.17504/protocols.io.5qpvokd4xl4o/v1), both `(manual)`
  in the source manifest; 10x Multiome (PN-1000283), 3′ (PN-1000268), 5′ (PN-1000263),
  GEM-X 5′ v3 and VDJ user guides; Supplementary Information PDF (`MOESM1`, image-only text).

## Protocol boundary

🟢 UDA-seq is an indexing wrapper rather than one analyte-specific chemistry. Fixed cells or
nuclei are intentionally overloaded into standard 10x droplets and receive a droplet
(round-1) barcode, are released intact when the emulsion is broken, are mixed and
dispensed into 96-well plates, and receive a well-specific (round-2) index by PCR. The
cell barcode is the 10x barcode combined with the post-index (Methods lines 433, 445). The
authors demonstrate RNA–ATAC (Multiome), RNA–VDJ and RNA–CRISPR configurations. The page
draws the Multiome route, which the Methods describe end to end.

## Oligos (Supplementary Table 7, primer summary) 🟢

```
Truseq-i5 index primer (96)  AATGATACGGCGACCACCGAGATCT-[8-bp sample index]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT
                             (each row ends ...CGCTCTTCCGATC-s-T: 3' phosphorothioate; first index GGTTTACT)
Nextare-i7 index primer (96) CAAGCAGAAGACGGCATACGAGAT-[8-bp sample index]-GTCTCGTGGGCTCGG   (first index CGCCGTTT)
P5 primer                    AATGATACGGCGACCACCGAGA
Partial TSO                  AAGCAGTGGTATCAACGCAGAG
Bio-Partial TSO              Bio-AAGCAGTGGTATCAACGCAGAG
P7 primer                    CAAGCAGAAGACGGCATACGAGAT
Tn5-top_ME                   [5Phos]CTGTCTCTTATACACATCT[3ddC]
Tn5-bottom_Read2N            GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
```

🟡 (computed with `lib/`): the i5 primer is `illumina.P5[:25]` + index + `illumina.TRUSEQ_READ1`;
it lacks the `ACAC` that follows `GATCT` in the 29-nt P5 of standard dual-index adapters.
The i7 primer is `illumina.P7` + index + `nextera.S7`. Tn5-top_ME is `nextera.ME_RC` and
Tn5-bottom_Read2N is `nextera.S7 + nextera.ME`, so the "i7-only" transposome inserts only
s7–ME ends. The P5 primer is `illumina.P5[:22]`; Partial TSO is `rt.SMART_HANDLE` minus
its 3′ T.

Published naming discrepancy (shown on the page): the Methods use "Partial P5",
"Partial TSO/IS", "Bio-Tso/ISPCR", "Tso/ISPCR" and "TN5_R2_index", whereas Table 7 prints
only "P5 primer", "Partial TSO", "Bio-Partial TSO" and "Tn5-bottom_Read2N". The model takes
the Table 7 oligo for each Methods name; whether "Partial P5" is a distinct unprinted oligo
remains open.

## Multiome UDA-seq, step by step 🟢 (Methods lines 325–395)

1. Fix nuclei: formaldehyde 1% (cell lines) or 1.5% (tissue), RT 10 min; wash twice.
2. Tagment **before** loading: 5 µL nuclei + 7 µL ATAC buffer B + 3 µL ATAC enzyme B
   (10x PN-1000283), 37 °C 60 min. (Results line 155 also says nuclei are tagmented before
   loading in the Multiome case.)
3. Overload the Chromium Controller (Multiome kit); in-GEM RT and ligation 37 °C 45 min,
   25 °C 30 min. Bead oligo, splint and TSO structures are 10x kit content 🔴.
4. Release: 125 µL recovery agent; remove the oil phase; add 120 µL PBS to the 80 µL aqueous
   phase; dispense 2 µL per well of a 96-well plate.
5. 1 µL proteinase K + 5 µL EB per well, 55 °C 10 min, 85 °C 5 min; Dynabeads SILANE
   purification, 10 µL EB.
6. Well-index PCR for both modalities: 15 µL NEBNext Ultra II Q5, 0.5 µM each Partial TSO/IS,
   Truseq-i5 index (well-unique), P5 primer and Nextera-i7 index (well-unique); 72 °C 5 min,
   98 °C 3 min, 7 × (98 °C 20 s, 63 °C 30 s, 72 °C 1 min), 72 °C 1 min. Half of the plate's
   products pooled, 1.6× XP, 160 µL EB.
7. ATAC: 40 µL + KAPA HiFi, Partial P5 + P7 primer; 98 °C 45 s, 7–9 × (98 °C 20 s,
   67 °C 30 s, 72 °C 20 s), 72 °C 1 min; 0.5–1.4× XP; MGISeq-T7 PE100.
8. cDNA: 40 µL + KAPA, Partial P5 + Bio-Tso/ISPCR, 4–5 cycles; MyOne Streptavidin C1 capture
   (B&W buffer, RT 60 min); on-bead PCR with Partial P5 + Tso/ISPCR, 4–5 cycles; 0.6× XP.
9. RNA library: i7-only Tn5 = cassette annealed 1:1 at 10 µM (75 °C 15 min, 60/50/40 °C 10 min
   each, 25 °C 30 min), 14 µL cassette + 8 µL TruePrep Tagment Enzyme + 78 µL coupling buffer,
   30 °C 1 h. Tagment 50 ng cDNA, 55 °C 10 min; 0.8× XP; NEBNext Q5 with Partial P5 +
   Nextera-i7 index, 72 °C 5 min, 98 °C 45 s, 8–9 × (98 °C 20 s, 60 °C 30 s, 72 °C 1 min),
   72 °C 5 min; 0.6–0.8× XP; MGISeq-T7 PE150.

## Final libraries and reads 🟡

Derived from the Table 7 primers and the step order; 10x regions inferred.

- RNA: P5[:25] · well index (8, round 2, i5 position) · TruSeq Read 1 · 10x barcode (16) ·
  UMI (12) · poly(dT) · cDNA · ME′ · s7′ · sample i7 (8) · P7′.
- ATAC: P5 · 10x barcode (16) · s5 · ME · gDNA · ME′ · s7′ · well index (8, round 2, i7
  position) · P7′.

So the round-2 well index sits in the i5 read for RNA and in the i7 read for ATAC. The
cell barcode is 10x barcode × post-index, and FASTQs were grouped by post-index (lines 433, 445).
Sequencing was on MGISeq-T7 (RNA PE150, ATAC PE100). 🔴 No Illumina-to-MGI conversion, read
primers or index-read lengths are given in the inspected text.

## Model history

- Before the 2026-10-17 bounded pass, the module came from shared helpers. It placed
  tagmentation *after* in-droplet RT and put both barcodes inside Read 1 of generic
  TruSeq/Nextera libraries. The pass replaced it with a local graph that follows the
  Methods order and the Table 7 primer positions.

## Open questions

- 🔴 Whether "Partial P5" is a distinct oligo; MGI conversion and read lengths of the index reads.
- 🔴 sc3′/sc5′ RNA, VDJ and CRISPR routes are not drawn; their step-by-step protocol is
  protocols.io ref 62 (not read). The sc RNA route uses the Truseq-i5 index + Partial TSO/IS
  well PCR (KAPA, 10–11 cycles) followed by the same i7-only Tn5 library (Methods 349–355).
- 🔴 ExoSAP-IT and direct-PCR purification alternatives for the RNA route (Methods 347).
