# Multi-CUT&Tag — research notes

## Sources read

- 🟢 Gopalan et al. 2021, DOI `10.1016/j.molcel.2021.09.019`, PMC full text and STAR Methods.
- 🟢 Supplementary Table 1 in the bioRxiv supplement (`2021.04.27.441642_media-1.pdf`): complete target-coded A/B oligos, Tn5MErev and custom sequencing primers.
- 🟡 Provenance of the sequences: the published Table S1 (PMC `NIHMS1743381-supplement-2.pdf`) sits behind PMC's
  download gate and Europe PMC reports the manuscript as not open access (checked 2026-10-17). The preprint
  Table S1 is used only because the published one could not be read; the published STAR Methods name the same
  `P5_i5_N_Universal_Connector_A` / `P7_i7_N_Universal_Connector_B` series and custom read/index primers. Any
  difference between the two Table S1 versions is unchecked.

## Reaction order (published STAR Methods, PMC8604773 XML text lines 161–181; bounded pass 2026-10-17)

- 🟢 Adapters A and B are each annealed at 45 µM (10 mM HEPES pH 7.2, 50 mM NaCl, 1 mM EDTA; 95 °C 5 min, 0.1 °C/min
  to 25 °C) and loaded onto 6-His pA–Tn5 with glycerol for 60 min at room temperature.
- 🟢 Loaded pA–Tn5 is coupled to its antibody (~2:1 pA–Tn5:antibody, 4 °C 4 h or overnight), then purified with
  Dynabeads His-Tag pulldown, eluted with 300 mM imidazole and buffer-exchanged on Amicon filters; the Results
  state this removes uncomplexed antibody and free adapters. The graph now draws this separation.
- 🟢 Bulk: 100,000 cells on concanavalin A beads, 2 µg of each complex in Dig-wash + 2 mM EDTA overnight at 4 °C,
  three Dig-med washes, 10 mM MgCl2 at 37 °C for 1 h, stop with EDTA/SDS/proteinase K (55 °C 1 h),
  phenol:chloroform extraction and ethanol precipitation.
- 🟢 Bulk PCR in NEBNext HiFi 2x master mix: 72 °C 5 min gap fill, 98 °C 30 s, 17 cycles of 98 °C 10 s / 63 °C 30 s,
  72 °C 1 min; 1.1x AMPure XP cleanup.
- 🔴/🟡 The barcoded i5 and i7 amplification primers are named but their sequences are not in the preprint
  Table S1, and the published Table S1 could not be read. The schematic's P5–i5 and i7–P7 arms are therefore
  inferred; their positions are constrained only by the custom Index2/Index1 primers, which read toward the
  library ends from the Read 1/Read 2 sites. (Earlier versions of this repository wrongly said Table S1 prints
  these primer families; that was a model error, not a publication error.)
- 🟢 Single-cell: no concanavalin A, Dig-wash washes with 600 g pelleting, equal volume of 20 mM MgCl2, then 10x
  scATAC GEMs with 72 °C 5 min and a single 98/59/72 °C cycle (vendor recommends more, to avoid
  incorporating uninserted adapters), SILANE/AMPure cleanup, 14-cycle sample-index PCR and double AMPure size
  selection.

## Exact adapter design

🟢 Four A-family and four B-family top oligos are printed in Table S1. Each contains a custom sequencing-primer site, an 8-nt antibody barcode, a conserved connector and the 19-nt Tn5 mosaic end. The common bottom strand is `5′-/5Phos/CTGTCTCTTATACACATCT-3′`; it pairs only to the mosaic end, leaving the barcode-containing part single-stranded in the loaded transposome.

🟢 A-family structure: `TCGTCGGCAGCGTCTCCACGC` + 8-nt antibody barcode + `GCGATCGAGGACGGC` + `AGATGTGTATAAGAGACAG`.

🟢 B-family structure: `GTCTCGTGGGCTCGGCTGTCCCTGTCC` + 8-nt antibody barcode + `CACCGTCTCCGCCTC` + `AGATGTGTATAAGAGACAG`.

🟢 Each antibody is pre-coupled to pA–Tn5 loaded with its matched barcoded adapters. Magnesium activation transfers those adapters to nearby chromatin. A fragment may carry the same antibody code at both ends (homogeneous) or different codes at its ends (mixed, interpreted as co-localization).

## Completion and sequencing

🟢 Bulk libraries receive a 72 °C, 5-minute gap fill, then 17 PCR cycles with indexed i5 and i7 primers. The source run is Read 1 72 cycles, Read 2 72 cycles, Index 1 8 cycles and Index 2 8 cycles. Both genomic reads begin with 8-nt antibody barcode + 15-nt connector + 19-nt mosaic end; the paper trims these first 42 bases before alignment.

🟢 The custom sequencing primers are Read1 `TCGTCGGCAGCGTCTCCACGC`, Read2 `GTCTCGTGGGCTCGGCTGTCCCTGTCC`, Index1 `GGACAGGGACAGCCGAGCCCACGAGAC`, and Index2 `GCGTGGAGACGCTGCCGACGA`.

🟢 For single-cell Multi-CUT&Tag, the target-coded fragments enter the 10x Chromium Single Cell ATAC workflow. The paper uses one GEM amplification cycle and 16 Index 2 cycles. The named [10x CG000209 v1.1 guide](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf) shows the gel-bead oligo's P5–16-base barcode–Read 1N architecture (Appendix A), and the vendor [Cell Ranger ATAC software](https://github.com/10XGenomics/cellranger-atac) distributes its barcode whitelist (`lib/python/barcodes/737K-cratac-v1.txt`). 🔴 The paper does not assign a particular whitelist barcode to a particular cell. The public schematic leaves that variable molecule abstract.
