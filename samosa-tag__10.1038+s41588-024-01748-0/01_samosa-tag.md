# SAMOSA-Tag — research notes

## Active native migration, 2026-10-18

The graph selects in-situ cell-line SAMOSA-Tag. A separate no-EcoGII control allocation is reserved before marking; it retains any native DNA modifications and is not assumed chemically unmethylated. EcoGII250U/37°C30min and source SAM replenishment mark accessible adenines on original genomic strands without assigning a fixed methylation pattern. The source pellet/resuspension/filter and a real10,000–50,000-nucleus tagmentation allocation are explicit; unallocated nuclei are not discarded merely by sampling.

SupplementaryTable5's representative SMRT-A_bc001 is a5′phosphorylated single DNA hairpin oligo. Its complete sequence is source-transcription tested; its27bp reverse-complementary stem is checked by HairpinTopology, and its logical8base sample barcode survives both hairpin arms and the native final component inventory. The adapter substrate is distinct from a completed SMRTbell template. Source hairpin-loaded Tn5 transposition leaves9nt gaps; the graph does not claim transposition alone closes a genomic circuit.

Source RNaseA37°C15min and ProteinaseK/SDS/EDTA60°C1–2h precede a2XSPRI recovery and a separate0.6XSPRI enrichment for>500bp. Coupled Phusion-HF/normal dNTP gap filling and TaqDNA-ligase nick sealing at37°C1h are separately explicit chemical transformations within the same source reaction. Short normal-DNA repair patches are distinct from both original methylated genomic strands; no whole-genome PCR replacement is introduced. Incomplete closure classes retain unresolved actual recovery. A separate2XSPRI recovery precedes ExoIII37°C1h; ExoIII hydrolysates/possible surviving open templates remain molecular classes until a final separate2XSPRI cleanup. Failure to close or amplify never itself physically discards a molecule.

The final scene inventories components of the source-known candidate closed native template, preserving both original strands, integrated hairpin features and short repair-DNA class. Exact oriented genomic/hairpin junction sequence, per-strand methylation positions and sequencing-primer binding polarity remain unassigned. SourceSequelII8MSMRTcells,30h movies/2h pre-extension/4h immobilization and size-dependent polymerase2.1/2.2 are distinguished from barcode width and short-read cycles. No new claim that sequencing-primer sequence is proprietary/absent follows this capped pass. Ex-situ/PDX, optional high/low-size selections and all other adapter members are separate experiments.

This capped pass read 3 originals within the selected round and reran the mandatory bundle primer scraper. All historical review statuses and the independent two-reader missing-information gate are unchanged.

## Historical observations

The earlier material below retains its original scope; current source-owned declarations above replace family-selected geometry and later-version execution defaults.


## Sources read

- 🟢 Dennis et al. 2024, DOI [10.1038/s41588-024-01748-0](https://doi.org/10.1038/s41588-024-01748-0), full text, Methods and Supplementary Note 1.

## Sections read (2026-10-17 bounded pass)

- PMC11176058 JATS (`s41588-024-01748-0_PMC11176058.xml.txt`) Methods lines 177–185
  (hairpin annealing and Tn5 loading), 195 (SMRT-Tag tagmentation, Repair Mix, ExoDigest
  Mix), 207 (optional AMPure PB size selection), 213–225 (nuclei isolation, SAMOSA
  footprinting, in situ and ex situ SAMOSA-Tag), 249 (QC and PacBio sequencing).
- MOESM7 Supplementary Table 5 rows 1–8 (SMRT-A_bc001 matches the page oligo base for base).
- MOESM1 Supplementary Fig. 1 legend lines 366–374 (Phusion + Taq ligase chosen for gap repair).

## Molecular order

🟢 Permeabilized nuclei are treated with the nonspecific adenine methyltransferase EcoGII
and SAM. 🟢 Hairpin-adapter-loaded Tn5 then fragments the marked chromatin in situ while
transferring PacBio-compatible hairpin/mosaic-end structures. 🟢 Purified DNA is gap
repaired and sequenced without PCR. The closed template preserves endogenous CpG
methylation and exogenous m6dA accessibility marks on the same molecule.

## Interpretation boundary

🟢 Supplementary Table 5 publishes 26 complete 5′-phosphorylated SMRT-A barcode
hairpins. The page uses SMRT-A_bc001 as a representative exact molecule: reverse-
complement mosaic end, reverse-complement 8-nt barcode, 45-nt PacBio hairpin loop,
barcode, and mosaic end. The shared hairpin constructor verifies that the 27-nt
barcode-plus-mosaic-end arms self-pair before rendering the nucleotide-resolved fold.

🔴 The separate PacBio sequencing-primer oligo remains proprietary. This does not make
the published hairpin sequence unknown. The page now distinguishes those two facts.
Tn5 transfer still leaves 9-nt gaps; polymerase/ligase gap repair is a separate reaction
before exonuclease selection of closed templates.

## Reaction conditions modelled (in situ SAMOSA-Tag)

🟢 Nuclei in Buffer M with 1 mM SAM are treated with 250 U EcoGII for 30 min at 37 °C (SAM
replenished to 1.16 mM at 15 min); half is reserved as an unmethylated control. Hairpin
oligos are annealed (95 °C 5 min, 25 °C 30 min) and loaded onto triple-mutant Tn5
(R27S,E54K,L372P) for 55 min at 23 °C. 10,000–50,000 nuclei are tagmented with 9.4–18.8 pmol
transposome in Omni-ATAC buffer at 55 °C for 45–60 min, stopped with RNase A then proteinase
K/SDS/EDTA (60 °C, 1–2 h), purified with 2× SPRI and enriched for >500 bp with 0.6× SPRI.
Repair Mix (2 U Phusion-HF, 80 U Taq DNA ligase, 0.8 mM dNTPs; 37 °C 1 h) fills the gaps
and seals the nicks; after 2× SPRI, ExoDigest Mix (exonuclease III, 100 U per 160 ng,
NEBuffer 2; 37 °C 1 h) digests unrepaired molecules, followed by a final 2× SPRI. Libraries
were sequenced on Sequel II 8M SMRT Cells (30 h movies) with polymerase 2.1 or 2.2.

🟡 Not drawn: the unmethylated control branch (same chemistry without EcoGII), the ex situ
variant (purified EcoGII-footprinted DNA tagmented with 0.046 pmol Tn5 and processed as for
SMRT-Tag), the PDX FACS front end and the optional AMPure PB size selection.
