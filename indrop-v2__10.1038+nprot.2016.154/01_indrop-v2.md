# inDrop v2 — research-note provenance

The bead synthesis, source reading, v1/v2 comparison, and read-layout audit are retained
in the combined [inDrop research note](../indrop-v1__10.1016+j.cell.2015.04.044/01_indrop.md).
This pointer avoids duplicating that evidence trail; the schematic is tied to Zilionis
*et al.* 2017 and its supplementary oligo tables.


**Bounded pass 2026-10-17.** Read Supplementary Methods 1–2 (MOESM454), the oligo table
(MOESM455 rows 1–36) and the MOESM460 pipeline README; the main Nature Protocols text
(doi:10.1038/nprot.2016.154) and its ICCB mirror remain unavailable. The schematic now
shows the Exonuclease I bead clean-up with the published 102 vs 101 nt band discrepancy,
treats empty/multiplet droplets as computationally excluded rather than discarded, and
no longer states acquired cycle counts (none appear in the inspected supplements). Details
and deferrals are in §10 of the combined note.

**Source follow-up 2026-10-17.** No legitimate open copy of the main text was found:
Unpaywall and OpenAlex report it closed with no repository full text; Europe PMC has no
PMCID; PubMed links no GEO series; Harvard DASH has no deposit; the Klein lab site links
only the DOI; the ICCB mirror still returns HTTP 403. The schematic is unchanged.

**Bounded pass 2026-10-17 (deferred-followup-03): main text read.** The publisher PDF of
Zilionis *et al.* 2017 (30 pages, printed doi:10.1038/nprot.2016.154; md5 in the source
manifest) is now in the source bundle and was checked against its printed DOI, title and
authors before use. Sections read: Fig. 3 and legend (library preparation), Materials
reagent oligos, Reagent Setup (PE2-N6, PE1/PE2 PCR primers, custom sequencing primers,
BA19), Procedure Steps 30–79 (bead extension, HinfI, BA19/ExoI clean-up, BA19 removal),
Fig. 6 legend, Steps 127–166 and Table 1.

- 🟢 Every oligo in the model matches the main-text Reagent Setup and MOESM455 bases,
  checked by computation (PE2-N6, PE2 P5 primer, PE1 index primer, custom Read 1 / Index /
  Read 2, acrydite primer, W1 probe). BA19 (`BAAAAAAAAAAAAAAAAAAA`) is now a graph input.
- 🟢 Bead clean-up: HinfI 37 °C 2 h (Step 66); BA19 hybridization (Steps 70–72); ExoI
  0.27 U/µl, 2 h at room temperature "to remove primers missing poly-T tails" (Step 74);
  BA19 removed by alkaline denaturation (Step 76). Fig. 6 calls the lower bands synthesis
  intermediates hydrolysed by ExoI. The purpose of the bead HinfI digestion is not stated.
- 🟢 Library steps: in-drop SuperScript III RT 50 °C 2 h, 70 °C 15 min (Step 128); ExoI +
  HinfI digestion 37 °C 30 min, beads removed on the filter, 1.2× AMPure (Steps 134–138;
  the target survives as a cDNA:mRNA hybrid with ssDNA T7 promoter, Fig. 3 legend);
  NEBNext second strand 16 °C 2.5 h (Step 140); HiScribe T7 IVT 37 °C 13–15 h (Step 142);
  zinc-ion fragmentation 70 °C 3 min (Step 148); PE2-N6 PrimeScript RT (Steps 151–154);
  Kapa HiFi PCR with 2 + qPCR-determined (typically 8–12) cycles (Step 161); 0.7× AMPure,
  400–600 bp average (Steps 163–164).
- 🟢 Reads: Read 1 is the transcript, Read 2 the cell barcode/UMI read of at least 50 bp
  (Step 166). Table 1 gives instrument-specific examples (for example NextSeq 37/6/51); the
  page lists all four rather than choosing one.
- Published naming slip shown on the page: Reagent Setup prints the indexed primer as "PE1
  primer" but refers to "24 PE2 primer variants"; Supplementary Table 2 names them PE1.
- The 102 vs 101 nt bead-QC band conflict stays within Supplementary Method 2. The main text
  gives only the full-length size (134–136 nt).

Still open: the bead HinfI purpose; whether the AMPure clean-ups remove the digestion
products; MOESM459 Hamilton methods; a second independent reader.


## Critical re-review, 2026-10-17

ExoI hydrolysis is distinguished from physical discard. Short-primer hydrolysis products remain unresolved because soluble versus bead-attached products and later wash partition are not reconstructed.

Inspected scope: nprot.2016.154.pdf.txt Fig. 6 and Procedure Steps 74–76; existing source-qualified bead-primer graph. This bounded correction does not complete the linked-source or two-reader absence audit.
