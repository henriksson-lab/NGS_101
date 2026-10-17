# Spear-ATAC — source intake

Defining paper: Pierce et al., *Nature Communications* (2021), DOI
[10.1038/s41467-021-23213-w](https://doi.org/10.1038/s41467-021-23213-w), read as
PMC8137922 (JATS `.txt`). Originals inspected in this pass: the article Methods
(extracted lines 100–124: plasmid cloning, 10x modifications, sgRNA amplification,
sgRNA alignment), Supplementary Note 1 (MOESM1 lines 2323–2502, cloning) and
Supplementary Note 3 (MOESM1 lines 2590–2721, targeted sgRNA amplification), and
Supplementary Data 8 — the primer table (MOESM11 xlsx, rows 7–11).

🟢 The protocol produces chromatin-accessibility and guide-identity libraries that retain
the same droplet cell barcode.

## Exact oligos (Supplementary Data 8, confirmed in the Methods and Note 3)

🟢 `oSP2053`, 5′-biotinylated:
`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTaagtatcccttggagaaccaccttg` — TruSeq Read 2
adapter plus a 25-nt sgRNA-specific 3′ half (printed lowercase in the source).
🟢 `oSP1594`: `AATGATACGGCGACCACCGAGA` (P5).
🟢 Indexed P7 primer:
`CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTG` — note that it stops after
the first 22 bases of TruSeq Read 2, which the model now reproduces (an earlier model
version used the full 34-base Read 2 sequence).
🟢 `oSP1735`: `gctacattttacatgataggcttgg` — 1.2 µL of 50 µM is spiked into the 10x GEM
master mix, and GEM PCR1 is run for 15 instead of 12 cycles. These are the only two
modifications to the 10x scATAC v1 kit (PN-1000110).
🟢 Custom guide Read 1 primer `oMCB1672`:
`GCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC`.
Sequencing from the Nextera Read 1 adapter instead is possible but needs ~150 cycles.

## Reaction order (Supplementary Note 3)

🟢 The sgRNA library is built **from the finished scATAC library** (~150 nM in 20 µL):
2.5 µL per sample is split over eight 100 µL 2× NEBNext reactions with biotinylated
oSP2053 and linearly amplified for 25 cycles (98 °C 10 s / 63 °C 15 s / 72 °C 20 s);
pooled and MinElute-purified into 40.5 µL; captured on 25 µL Streptavidin MyOne C1
beads in B&W buffer (2 M NaCl, 10 mM Tris pH 7.5, 1 mM EDTA), washed 2× B&W and 2×
water; PCR2 off the beads in eight 100 µL reactions with oSP1594 and the indexed P7
primer for 15 cycles (69 °C anneal); the beads are magnetised away and the
**supernatant** is MinEluted into 20.5 µL EB and KAPA-quantified.

🟢 The note warns that exponential amplification before the streptavidin step, or no
streptavidin step at all, contaminates the sgRNA library with the scATAC library.

## What is not source-exact

🔴 The pSP618 cassette between the 10-bp sgRNA barcode and the spacer is described only
in words (mouse U6 promoter carrying a 34 bp Nextera Read 2 adapter in place of the loxP
site, the Perturb-seq `cr1` constant region, a constant region and a 34 bp Nextera
Read 1 adapter). No inspected original prints those bases, so the model draws a
placeholder there. 🟡 The sequence downstream of the spacer is derived as the complement
of the sgRNA-specific 3′ half of oSP2053, which must be present for that primer to
anneal. The Addgene plasmid pSP618 / GenBank MW852482 would settle the interior and was
not read in this pass.

🔴 The 10x scATAC v1 gel-bead oligo (P5 – 16-nt cell barcode – S5) and the indexed N7xx
ATAC primer are vendor material: named in the Methods, not printed in the inspected
originals. They stay declared by role with vendor-side bases from the canonical
Illumina/Nextera modules.

🔴 No instrument, read length or cycle count is given for either library. The 16-nt
droplet cell barcode is acquired as an index read and is the identifier that links a
guide to its accessibility profile; the i7 index is used for sample demultiplexing.

## Open questions

- Addgene pSP618 / GenBank MW852482 for the exact spacer context and barcode position.
- The 10x Chromium scATAC v1 (PN-1000110) user guide revision used, for the bead oligo
  and the N7xx index primers.
- Whether the 10-bp sgRNA-barcode read-out route (direct from the scATAC library) was
  used for any library reported here, or only the targeted amplification route.
