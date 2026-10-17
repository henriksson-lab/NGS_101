# smol-seq

**Defining paper:** Tan JH, Fraser AG. *Quantifying metabolites using structure-switching aptamers coupled to DNA sequencing.* *Nature Biotechnology* 44 (2026), doi:10.1038/s41587-025-02554-7. Read from the PMC Cloud open-access JATS (PMC12807852, `s41587-025-02554-7_PMC12807852.xml.txt`) and the published Supplementary Data 1 (`…MOESM2_ESM.xlsx.txt`, sheet `SuppData1 Oligonucleotides`). The MOESM1 PDF is the three-page Nature Portfolio reporting summary.

## Sections read (2026-10-17 bounded pass)

- JATS line 32 (Main, SSA design), line 38 (Fig. 1a legend), line 51 (multiplexing), Methods lines 61 (reagents), 65 (ligand-binding assay), 77 (SRO readout), 81–89 (multiplexing and sequencing of bSSAs).
- Supplementary Data 1 lines 56–75 (mefloquine, ATP and piperaquine bSSAs; round-1 and round-2 PCR primers).

## Molecular design

🟢 Each structure-switching aptamer (SSA) is two oligos: a longer ligand-binding oligo (LBO), base-paired to a short release oligo (SRO) through a short stem. Cognate ligand drives the LBO into a stable stem-loop, which displaces and releases the SRO. The LBO is immobilized on a solid matrix, so ligand binding releases the SRO from the matrix (line 32). Each barcode SSA (bSSA) has an SRO with a barcode unique to its target.

🟢 Mefloquine LBO (Supplementary Data 1): `/5Biosg/ATACCGCTCTCGGGACGACGGCAGTCTATACCCCGTATCGCCGAAGGTTGTCGTCCC`. Mefloquine barcode 1 SRO: `CGTCGTCCCGAGAGTAGC AGTCTGAC N5-7 AGATCGGAAGAGCGTCGTGTAG`. Barcodes 2 and 3 are `CTGATGCA` and `CGATCATG`; ATP and piperaquine bSSAs follow the same SRO pattern with their own stems and three barcodes each.

🟢 SRO design (Methods line 81–83): `5′ stem-specific sequence – 8-nt barcode – N5–7 – AGATCGGAAGAGCGTCGTGTAG 3′`; the 3′ end is a PCR handle "complementary to the Illumina P5 adaptor sequence" and the N5–7 random bases help sequence these low-diversity amplicons. They are not described as molecular identifiers.

🟡 Computed with `chemdraw.Scene` in the model: the MQ SRO's 5′ 14 nt pair with LBO bases 7–20 (`chem.LBO_PAIRED`, checked by `sensor_scene()`); the SRO's `TAGC`, barcode, N5–7 and handle are single-stranded. The LBO 3′ end (`…GTCGTCCC`) can also pair with that LBO stem region, consistent with the described ligand-stabilized stem-loop displacing the SRO; the fold itself is not drawn.

## Step by step (Methods lines 65, 85–87)

1. 🟢 25 pmol LBO + 125 pmol SRO heated to 95 °C for 5 min, slowly cooled to 25 °C.
2. 🟢 Bound to 5 µL Dynabeads MyOne Streptavidin C1 for 30 min; washed twice, incubated 45 min, washed once.
3. 🟢 In the sequenced experiment, ATP, PQ and MQ bSSAs (three barcodes each) were mixed and incubated 45 min with 125 µM ATP, 30 µM MQ, 25 µM PQ or PBS; the supernatant is the "released" sample.
4. 🟢 The remaining bound SROs were boiled off (ligand-binding assay: Mg-free strand separation buffer, 95 °C 4 min) to give the "boiled" sample used to compute percentage release. The retained fraction is therefore sequenced, not discarded.
5. 🟢 Replicates (marked by a short identifier sequence) were pooled, spike-in controls (ten sequences, twofold titration) added, and the single-stranded DNA purified with the Monarch kit oligonucleotide protocol.
6. 🟢 Round-1 PCR adds partial P5/P7 handles (14 cycles released, 11 cycles boiled), primers `ACACTCTTTCCCTACACGACGCTCTTCCGATC` (Reverse P5) and `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` + stem (e.g. Forward MQ `…CGTCGTCCCGAGAGTAGC`). Monarch PCR cleanup, 20 µL elution, 1:5 dilution.
7. 🟢 Round-2 PCR (5 µL template, 100 µL, 8 cycles) with `AATGATACGGCGACCACCGAGATCTACAC [i5 index] ACACTCTTTCCCTACACGACG` and `CAAGCAGAAGACGGCATACGAGAT [i7 index] GTGACTGGAGTTCAGACGTG`; Monarch cleanup; NucleoMag NGS size selection at 1.2× bead:sample.
8. 🟢 MiSeq, ~4 million 50-bp single-end reads.

The model checks that both round-2 primer 3′ ends land on the round-1 product.

## Final library and read

🟡 Derived: `P5 – i5 – ACACTCTTTCCCTACACGACGCTCTTCCGATCT – N5–7′ – barcode′ – stem′ – AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC – i7′ – P7′`. The Read 1 site is assembled from the Reverse P5 primer plus the first templated T; Read 1 therefore reads the diversity bases, then the 8-nt analyte barcode (reverse complement), then the stem. The paper does not name the sequencing primers or index-read cycles: the TruSeq Read 1 primer is inferred, and the 8-nt i5/i7 windows are drawn as inferred placeholders.

## Model corrections (2026-10-17)

The earlier generic schematic described a quenched capture strand, a separate "reporter-release condition", a single PCR and discarded unbound sensor. It now draws the published LBO/SRO pair on streptavidin beads, ligand-driven release, the sequenced boiled-off remainder, both PCR rounds with source primers, cleanups and the 50-bp single-end read.

## Open questions

- PCR polymerase and spike-in sequences are not given in the inspected Methods/Supplementary Data 1.
- Index lengths/values and index-read cycles are not stated.
- Predecessor SSA sources (Yang et al., Coonahan et al., etc., line 61) were not followed in this capped pass.
