# APEX-seq

Fazal et al.(2019), *Atlas of Subcellular RNA Localization Revealed by APEX-seq*. DOI https://doi.org/10.1016/j.cell.2019.05.027 .

🟢 Associated accepted manuscript PMC6786773, provenance71; Results93–99 direct RNA labeling versus protein crosslink capture, Methods223/229/245/249 inspected. Main published VOR applicability unverified; linked NIHMS chemistry supplements are gated. Automatically fetched preprint454470 is superseded and not read or used. Mandatory full acquired-bundle scraper run; relevant sequencing preparation uses named historical kit, while extracted5S/qPCR sequences belong to validation, not sequencing-library adapters. No absence claim.

The source-qualified route uses APEX2-expressing cells, BP500µM30min37°C, then H2O2final1mM1min. Quench5mMTrolox/10mMascorbate/10mMazide, wash3times with Trolox/ascorbate DPBS. Harvest/pellet cells, RNeasyPlusMini RNA extraction with gDNA-eliminator column and RWT wash substitution. Neither microscopy fixation nor the separate dot-blotDNase/UV procedure is inserted into the sequencing route.

Direct RNA labeling is distinct from APEX-RIP protein labeling and formaldehyde capture. Guanine/phenol adduct characterization and selected RT stops are validation evidence, not a guarantee of RNA site specificity or an obligatory sequence mismatch. The exact RNA-adduct structure and copying distribution remain unassembled.

Source RNA is captured on Pierce streptavidin magnetic beads2h4°C and washed3times inB&W. ProteinaseK/sarcosine42°C1h then55°C1h liberates RNA-containing preparation; physical bead/protein remnant removal is not inferred from proteolysis alone. RNA Clean&Concentrator-5 purification follows. TruSeqstrandedmRNA preparation includes polyA selection; MiSeqQC precedes HiSeq4000PE75. Total-RNA MITO experiments omit polyA selection as a separate branch. Historical kit internals and exact oligos/strand/read geometry remain unreviewed and opaque. No whole-recipe import from APEX-RIP.
