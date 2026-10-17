# Illumina DNA Prep

## Sources read

🟢 Illumina, *Illumina DNA Prep Reference Guide*, document 1000000025416 v12, pp. 1–5 and 21–31 (workflow, input recommendations, product contents, Tagment Genomic DNA, Post Tagmentation Clean Up, Amplify Tagmented DNA, Clean Up Libraries, Pool Libraries, Check Library Quality). It defines bead-linked tagmentation, the TSB stop and TWB washes, input-dependent limited-cycle PCR, double-sided bead cleanup and pooling.

🟢 Illumina, *Illumina Adapter Sequences*: the guide (p. 5) refers index-adapter bases to this reference. Its umbrella landing page (bulletin URL in `ref/MANIFEST.md`) returned HTTP 404 on 2026-10-17; the child page *Sequences for Illumina Prep and Nextera Kits* and the *IDT for Illumina UD Indexes* page were read from the shotgun-metagenomics source bundle. They print the Read 1/Read 2 transposase adapters, the index-PCR primer arms `CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG` and `AATGATACGGCGACCACCGAGATCTACAC[i5]TCGTCGGCAGCGTC`, and state that UD index adapters are 10 bases. Canonical sequences are imported through `lib/nextera.py` and `lib/illumina.py`.

## Molecular workflow

1. 🟢 BLT and TB1 (11 µl each per sample, 20 µl dispensed) are mixed with 1–500 ng DNA in 30 µl; TAG program 55 °C 15 min.
2. 🟢 10 µl TSB is added and the PTC program runs 37 °C 15 min. On the magnet the supernatant is discarded, two 100 µl TWB washes follow, and a third TWB volume is left on the beads until PCR.
3. 🟢 TWB is removed; 40 µl EPM/water master mix and 10 µl pre-paired i7/i5 index adapters are added to the beads. BLT PCR: 68 °C 3 min, 98 °C 3 min, X cycles of 98 °C 45 s / 62 °C 30 s / 68 °C 2 min, 68 °C 1 min. The guide states that this PCR adds the i7 and i5 index adapters and the cluster-generation sequences.
4. 🟢 Double-sided IPB cleanup (standard gDNA): 45 µl supernatant + 40 µl water + 45 µl IPB; 125 µl supernatant moved to 15 µl IPB; first plate discarded; two 80 % ethanol washes; elution in 32 µl RSB, 30 µl recovered.
5. 🟢 Typical libraries span 250–1000 bp, with average fragment size about 600 bp on the guide's example traces.

🟡 The purpose of the initial 68 °C 3 min step is not stated by the guide; it is drawn only as a program step, not as an established gap fill.

## Input-dependent behavior

🟢 The guide accepts 1–500 ng or more. It recommends 100–500 ng for human and other complex genomes, while small genomes can use 1 ng.

🟢 PCR cycles by total DNA input: 1–9 ng → 12, 10–24 ng → 8, 25–49 ng → 6, 50–99 ng → 5, 100–500 ng → 5, blood/saliva → 5 (transcribed in `tools/chem.py` as `PCR_CYCLES`).

🟢 At 100–500 ng, the finite bead-linked transposome population provides normalized yield suitable for equal-volume pooling (< 20 % CV). Below 100 ng, libraries are not normalized and each library is quantified individually. The saturation limits yield; unbound excess DNA is physically removed with the post-tagmentation supernatant and washes, which is how the graph now records its fate.

## Index adapters and reads

🟢 Compatible index plates (p. 5): IDT for Illumina DNA/RNA UD Indexes, IDT for Illumina Nextera DNA UD Indexes, or Nextera DNA CD Indexes. The guide (p. 26) says UD codes are 10 bp, unlike the 8 bp Nextera DNA CD codes, and that the sequencing system must be configured for 10 bp index reads.

🟡 The schematic draws the UD variant with 10-base symbolic i5/i7 windows; a CD-indexed library has 8-base windows. The guide does not fix paired read lengths, so none are declared.

## Open questions

🔴 Illumina does not disclose the reagent-level bead/transposome composition. The schematic therefore treats BLT as reaction geometry while using the public Nextera endpoint and computed standard primer sites.

- The *Illumina Adapter Sequences* umbrella document (1000000002694) at its bulletin URL is unreachable (HTTP 404); only its child pages were read.
- Blood, dried-blood-spot, saliva and small-amplicon (< 500 bp, 81 µl IPB) branches are not drawn.
