# MulTI-Tag

**Defining source:** Meers, Llagas, Janssens, Codomo & Henikoff, *Nature Biotechnology*
(2022), doi:10.1038/s41587-022-01522-9 (PMC10188359, CC BY). Published JATS and
Supplementary Table 1 (MOESM2) are the authority; the bioRxiv preprint
10.1101/2021.07.08.451691 is superseded and not used for active claims.

## Sources read (bounded pass, round 2026-10-17-breadth-first-13)

- Published JATS `s41587-022-01522-9_PMC10188359.xml.txt`: Fig. 1a legend (line 49),
  Results single-cell paragraph (line 53), Discussion (line 85), Methods lines 105–121
  ("CUT&Tag", "Conjugates for MulTI-Tag", "Bulk MulTI-Tag protocol", "Single-cell
  MulTI-Tag", "Sequencing and data pre-processing").
- Supplementary Table 1: published MOESM2 is byte-identical (MD5
  `5cd0f4d7b6de26836d8c04dac29ecb93`, source `MANIFEST.tsv`) to the local
  `2021.07.08.451691_media-2.xlsx`; sheets README, P5_i5 adapters, P7_i7 adapters,
  i5 primers, i7 primers were read in its `.txt` twin.
- Not read: the published Supplementary Information PDF (MOESM1) has no text layer;
  the cited CUT&Tag protocols.io (bcuhiwt6) is on disk but was not re-read; ICELL8
  reference 62 (single-cell denaturation/amplification) was not followed.

## Oligos (Supplementary Table 1, verbatim)

🟢 Tn5MErev: 5′-[phos]CTGTCTCTTATACACATCT-3′ (Methods "Conjugates"; table README).
P5_i5 and P7_i7 adapters are annealed to it and loaded into pA–Tn5; i5 and i7 primers
are for amplification.

🟢 P5_i5 adapters (target-coded; Methods: "5′ aminated, barcode-containing
oligonucleotides" — the table does not print the amine). Example `P5_i5_1` (H3K27me3):
5′-TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG-3′, barcode TATAGCCT. Layout:
TCGTCGGCAGCGTCTCCACGC – 8-nt target barcode – GCGATCGAGGACGGC – AGATGTGTATAAGAGACAG.
Rows 1, 2, 3, 5, 6 are listed; `P5_i5_2` is assigned to PolIIS5P in one experiment and
H3K4me1 in another.

🟢 P7_i7 adapters: rows 1–12 for bulk, rows 13–108 (96) for single-cell. Example
`P7_i7_13`: 5′-GTCTCGTGGGCTCGGCTGTCCCTGTCCGTCTACATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG-3′, barcode GTCTACAT. Layout:
GTCTCGTGGGCTCGGCTGTCCCTGTCC – 8-nt well barcode – ACCGTCTCCGCCTCA – AGATGTGTATAAGAGACAG.

🟢 72 i5 primers (`sci_Ad1.n`) and 72 i7 primers (`sci_Ad2.n`), all dispensed on the
ICELL8 chip for single-cell work. Examples: `sci_Ad1.1_TAGATCGC`
5′-AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCTCCACGC-3′; `sci_Ad2.1_TAAGGCGA` 5′-CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGCTGTCCCTGTCC-3′.

🟡 Derived (checked by computation): both adapters end in the reverse complement of
Tn5MErev (the mosaic end); the i5 primer 3′ end equals the first 21 nt of every P5_i5
adapter, and the i7 primer 3′ end equals the first 27 nt of every P7_i7 adapter.

## Step by step (single-cell branch, Methods lines 109–119)

1. 🟢 Anneal 200 µM adapter to equimolar Tn5MErev in 1× PBS; conjugate to 100 µg
   antibody (Abcam Oligonucleotide Conjugation Kit, random amine conjugation, ~10:1 adapter
   excess).
2. 🟢 Per target: 0.5 µg conjugate + 1 µl ~5 µM pA–Tn5 + 16 pmol unconjugated annealed
   i5 adapter of the same sequence, 30–60 min RT (conjugate i5 transposome). Separately
   1 µl pA–Tn5 + 32 pmol i7 adapter (i7 transposome). Used within 24 h.
3. 🟢 Nuclei on beads (bulk: 50,000 on ConA; single-cell: 1.8 million on biotin-WGA /
   streptavidin T1 Dynabeads). First conjugate transposome in Wash-300 + 2 mM EDTA,
   1 h RT or overnight 4 °C; three Wash-300 washes; tagment in Wash-300 + 10 mM MgCl2,
   1 h 37 °C; 5 min Wash-300 + 5 mM EDTA. Repeat for each conjugate; the published orders
   always end with H3K27me3.
4. 🟢 Secondary antibodies (species of the conjugates), 1 h RT, two Wash Buffer washes;
   i7 transposome in Wash-300 1 h RT (single-cell: nuclei split across 96 uniquely
   barcoded i7 transposomes), three washes, final tagmentation with 10 mM MgCl2 1 h 37 °C.
5. 🟢 Single-cell: re-aggregate, TAPS washes, RQ1 DNase 10 min on a 0 °C block, stop
   buffer, two 20-µm strainers, DAPI + ICELL8 Second Diluent, 857 nuclei/µl, 640 µl
   dispensed into an ICELL8 chip; SDS denaturation, Triton X-100 neutralization and
   amplification in microwells "as described previously" (ref. 62); re-aggregate, two
   1.3× HighPrep cleanups.
6. 🟢 Bulk: TAPS, TAPS-SDS 58 °C 1 h, Triton X-100, primers + 2× NEBNext Master Mix:
   58 °C 5 min, 72 °C 5 min, 98 °C 30 s, then 98 °C 10 s / 60 °C 10 s with "repeat steps
   4–5 14 times", 72 °C 2 min; 1.1× HighPrep.

## Final library and reads

🟡 Role order derived from the oligo layouts: P5 – i5 primer barcode – target barcode –
ME – insert – ME – i7 transposome-well barcode – i7 primer barcode – P7. Cell identity
is the combination of i7 well (96), ICELL8 i5 primer and i7 primer.

🟢 Illumina HiSeq, paired-end 25 × 25 (Methods "Sequencing and data pre-processing").

## Open questions

- 🔴/unread: which sequencing primers read Read 1 and the index reads, and the index-read
  lengths. The target barcode lies between TCGTCGGCAGCGTCTCCACGC and the mosaic end, so
  a standard Nextera Read 1 primer is not implied. Not stated in the inspected Methods;
  MOESM1 and ref. 62 are unread. This is not a confirmed absence.
- Single-cell PCR mix and cycling (cited to ref. 62).
- Whether "repeat steps 4–5 14 times" means 14 or 15 cycles in total.
- The protocols.io DOI 10.17504/protocols.io.6qpvr4bezgmk/v2 that an earlier note
  version named was not found in any inspected source text and is no longer asserted.
- Bulk branch (12 bulk P7_i7 rows, tube PCR with specific primer pairs) is described
  here but not drawn.
