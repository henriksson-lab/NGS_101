# ScISOr-seq2 — research notes

## Sources read

- 🟢 Joglekar et al. 2024, *Nature Neuroscience*, DOI `10.1038/s41593-024-01616-4` (PMC11156538 author-manuscript JATS, local `s41593-024-01616-4_PMC11156538.xml.txt`): Methods "Linear/asymmetric PCR to remove non-barcoded cDNA" (line 212), "Exome capture to enrich for spliced cDNA" (line 216), "10× 3' library preparation and sequencing" (line 220), "Library preparation for PacBio" (line 224) and "Library preparation for ONT" (line 228).
- 🟢 Springer MOESM1 (supplementary figures/notes) and MOESM2 (reporting summary): scraper and keyword pass found no library oligo or reaction detail beyond the Methods.
- Not read: the PacBio SMRTbell Express Template Prep Kit 2.0 manual (100-938-900), the ONT SQK-LSK110 "Amplicons by Ligation" protocol ACDE_9110_v110_revC_10Nov2020, the Agilent SureSelectXT HSQ manual and the 10x CG000183 Rev C user guide.

## Construct path

🟢 Input: 30 ng of 10x Chromium Single Cell 3′ v3 cDNA (12 cycles of 10x cDNA amplification), so the cell barcode and UMI come from the original 10x RT product (lines 212, 220).

🟢 LAP: 12 linear/asymmetric cycles with Partial Read1 only (95 °C 3 min; 98 °C 20 s, 64 °C 30 s, 72 °C 60 s), 0.8× SPRIselect; then 6 exponential cycles with Partial Read1 + Partial TSO under the same conditions, 0.6× SPRIselect, eluted in 30 µl EB. KAPA HiFi HotStart ReadyMix is the polymerase for all these PCRs (line 212).

🟢 CAP: Partial Read1 and Partial TSO (`CTACACGACGCTCTTCCGATCT`, `AAGCAGTGGTATCAACGCAGAGTACAT`) are mixed as block oligos and bound to the cDNA ends (95 °C 5 min, 65 °C hold). Agilent SureSelect Human All Exon V8 or Mouse All Exon probes with SureSelectXT HSQ reagents hybridize for 24 h at 65 °C. The hybrids are captured on M-270 streptavidin Dynabeads and washed in WB1 (room temperature) and WB2 (65 °C, three changes). Bead-bound cDNA is then amplified for 12 cycles with Partial Read1 + Partial TSO (98 °C 20 s, 64 °C 60 s, 72 °C 3 min); the supernatant is kept and purified with 0.8× SPRIselect (line 216). The same two primers act as blockers and as PCR primers; they are not additional undisclosed oligos.

🟢 PacBio: ~500 ng LAP-CAP cDNA, SMRTbell Express 2.0: DNA damage repair 37 °C 30 min, end-repair/A-tailing 20 °C 30 min then 65 °C 30 min, adapter ligation 20 °C 60 min, 0.6× SPRIselect (line 224).

🟢 ONT: ~75 fmol LAP-CAP cDNA, SQK-LSK110 following protocol ACDE_9110_v110_revC_10Nov2020, PromethION flow cell FLO-PRO002 for 72 h, and Guppy basecalls with quality >7 (line 228). The page draws ONT end-prep plus adapter ligation as one vendor step with unavailable adapter geometry.

🟡 The PacBio hairpin T-overhang/dA geometry is the standard SMRTbell Express mechanism and is drawn as inferred; the kit manual was not read in this pass.

🔴 The paper names the PacBio and ONT kits but does not print their adapter or sequencing-primer bases. Those terminal commercial regions remain explicit placeholders.
