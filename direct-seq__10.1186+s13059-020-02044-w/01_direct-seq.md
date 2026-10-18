# Direct-seq — source notes

Defining paper: Song et al., *Genome Biology* (2020), [doi:10.1186/s13059-020-02044-w](https://doi.org/10.1186/s13059-020-02044-w).

- 🟢 PMC full text, figures and supplementary oligo tables were read.
- 🟢 A 30-nt mixed A/G “8A8G” tract is inserted into the sgRNA scaffold so ordinary poly(dT) primers capture the Pol III guide transcript.
- 🟢 In the 10x 3′ v3 workflow, guide and endogenous RNA acquire the same cell barcode and UMI during RT.
- 🟢 Pre-amplified cDNA is split: one branch makes the normal mRNA library and nested guide-specific PCR enriches a guide-index library before P5/P7 completion.
- 🟢 The exact capture tract and terminal PCR primers are printed in the supplementary tables.
- 🟢 The 10x Chromium 3′ v3 CRISPR Feature Barcode guide CG000184 Rev A,
  oligonucleotide appendix, prints the gel-bead Read 1 handle
  `CTACACGACGCTCTTCCGATCT`, followed by a 16-base cell barcode, 12-base
  UMI and T30. The capture drawing and graph now use this vendor-published
  handle rather than an unknown-length placeholder. The paper cites CG000184
  but does not identify individual barcode members.

## Programmed transcript and capture

🟢 Supplementary Tables S2–S3 give the complete Tail-8A8G scaffold. The 30-nt `AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA` tract sits between the programmed scaffold and Pol III terminator. In the 10x 3′ implementation a published 72-nt human Gln tRNA precedes the guide spacer and supplies the nested-PCR priming site.

🟡 The 10x poly(dT) bead primer can anneal across the capture tract, but the three encoded G residues are intentional mismatches. The schematic therefore places the exact oligos without claiming a perfect 30-bp duplex. Reverse transcription still gives the guide transcript the normal 16-nt cell barcode and 12-nt UMI.

## Nested guide-library PCR

🟢 The first enrichment PCR uses `P5_Read1` and `tRNA_Read2`. The 3′ 23 nt of `tRNA_Read2` exactly match bases 25–47 of the human Gln tRNA; its 5′ 34 nt are the TruSeq Read 2 primer site. A third PCR retains `P5_Read1` and replaces the reverse primer with one of two printed P7/index/Read-2 primers.

🟡 The final construct therefore makes its read structure computable: Read 1 begins with the ordinary 16-nt cell barcode plus 12-nt UMI. Read 2 positions 1–23 are the tRNA primer site, 24–48 the remaining tRNA sequence, 49–68 the guide spacer and 69 onward the programmed scaffold. The guide spacer—not the 8A8G tract—is the perturbation identifier.

🟡 Sequencing instrument, exact native primer assignments and acquired run lengths remain
unassigned by this bounded model. The model previously declared acquired Read 1 = 28 and Read 2 = 90 cycles;
those were unsourced 10x-default values and were removed in the 2026-10-17 bounded pass.
Positions above are derived spans, not acquired-cycle claims.

## Library-construction conditions (Methods, 10x 3′ platform)

🟢 PMC7278172 JATS lines 132–144 ("scRNA-seq on the 10 × 3′ RNA-seq platform"):
Chromium Single-Cell 3′ v3 (PN-1000075), B Chip (PN-1000153) and i7 Multiplex Kit
(PN-120262); CG000184 followed through cDNA pre-amplification (the 1st PCR). The cDNA is
divided into two 40 µl aliquots. Aliquot A is size-selected following the 10x protocol:
the 0.6× left-sided products make the mRNA library and the 0.6–1.2× double-sided products
(eluted in 25 µl) make index gRNA library I. Aliquot B is purified with 1.2× AMPure (25 µl)
for index gRNA library II. The 2nd PCR is eight 50 µl reactions per sample (3 µl template,
NEBNext Ultra II Q5 Master Mix M0544S, 0.5 µM each tRNA_Read2 and P5_read1; 98 °C 30 s,
14 × [98 °C 10 s, 60 °C 10 s, 72 °C 10 s], 72 °C 2 min), pooled and size-selected
0.7–1.0× (80 µl). The 3rd PCR is five 50 µl reactions (10 µl 2nd-PCR product,
P7_Read2_Index1 or Index2 plus P5_read1; 98 °C 30 s, 5 × [98 °C 10 s, 54 °C 15 s,
65 °C 20 s], 72 °C 2 min), then 0.7–1.0× AMPure XP. The page draws library I; library II
is a declared unresolved terminal receiving the same PCRs. The detailed Protocol Exchange
protocol (doi:10.21203/rs.3.pex-952/v1) and the Fluidigm C1/TruePrep branch remain unread.


## Endpoint source correction (2026-10-16)

The gene-expression and guide endpoints are now assembled locally. The earlier
shared endpoint recipe inserted a Nextera Read 1 arm even though the
paper-cited 10x bead handle and Table S3d `P5_Read1` both identify the TruSeq
arm. This was a repository model error. Construction now requires the
independently printed Table S3d PCR pair to have inward sites on the final
guide endpoint; substituting Nextera fails that check.

The representative guide cDNA now retains the primer-derived poly(dT)30
tract followed by upstream scaffold cDNA. The encoded 8A8G capture tract is
the site to which the bead primer anneals, not an additional full tract to
append after that primer. Its mismatched G:T capture geometry remains
explicit in the capture view; individual priming registers are not measured
by the schematic.

## Native graph migration, 2026-10-18

Read three published originals: article Methods (PMC7278172 JATS 132–144, 10x platform/aliquots and both guide PCRs), Table S2 (MOESM3, Tail-8A8G and human Gln tRNA rows), and Table S3d (MOESM4, P5_Read1/tRNA_Read2/P7_Read2_Index1/Index2 rows 45–48). Ran the bundle oligo scraper and verified relevant pointers in these originals. The separately cached CG000184 Rev A bead scaffold retains its prior independently reviewed provenance; the vendor guide was not re-read as a fourth original, nor was the detailed Protocol Exchange protocol followed in this capped pass.

The cell's guide RNA and endogenous mRNA are independent Scene populations, as are their single-stranded RT copies and the preamplified populations. Both source 40 µL aliquots retain the mixed pool; neither aliquoting nor size selection is treated as sequence-specific purification. Empty/unresolved GEMs have unresolved physical recovery, while removed bead supernatants remain physical waste. The mRNA library is a separately sequenced manufacturer-route boundary, replacing the generic terminal PCR and unsupported exact transcript endpoint.

The old first guide-PCR template inserted a 20-nt P5_Read1 landing region without reconstructing the native outer preamplification end. The source bead-derived handle instead shares only a 13-nt terminal overlap with that primer. The new template retains the known bead-derived sequence, scaffold, guide identity and tRNA sequence, with explicit unknown outer ends; it does not extend that handle to manufacture a longer PCR site. Actual first-enrichment PCR geometry is therefore unavailable on the incoming size-selected mixture. The local tRNA primer-binding view establishes its named local relation only. The resulting candidate amplicon is assembled from the printed primer ends and copied guide/tRNA interval, preserving barcode/UMI/guide features. The indexed PCR retains exact geometry on this actual first-PCR product, not on the expected final library. Both printed eight-base index-primer alternatives remain in the oligo inventory; Index1 is the represented indexing example.

Source modification, index and read semantics remain local. No acquired cycle counts or absent-information conclusion is introduced. Versioned vendor cDNA/mRNA-library chemistry, native first-PCR end geometry, the second guide-library copy and the separately scoped C1 route remain open targets.

Independent implementation review corrected the first guide-amplicon assembly: P5_Read1 shares the known 13-nt terminal overlap with the reviewed bead handle, so copying retains the remaining nine handle bases before the cell barcode. The suffix is derived computationally from those two sequences. A constructor checks this prefix in both first and indexed guide products; a negative regression rejects deletion of the retained suffix. This establishes the local source-sequence relation only. First-PCR geometry remains unavailable on the actual mixed input, and second-PCR exact geometry checks its actual corrected first amplicon. No new original was read; the migration source scope above and prior CG000184 bead-scaffold review apply.
