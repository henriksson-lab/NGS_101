# SHAPE-Seq

**Defining source:** Lucks *et al.*, *PNAS* (2011),
doi:10.1073/pnas.1106501108. The open article and experimental overview were read.

## Evidence used

- 🟢 Each designed RNA carries a DNA barcode at its 3′ end before pooling.
- 🟢 SHAPE chemistry acylates flexible 2′ hydroxyls; primer-extension stops encode the
  modification sites relative to a minus-reagent control.
- 🟢 A second adapter is ligated to cDNA and paired-end sequencing connects each stop
  with the barcode identifying the RNA template.
- 🟡 Exact oligo sequences have not yet been transcribed into the shared model.

## Construct logic

The original protocol is a multiplexed in-vitro accessibility assay for defined RNA
pools. It is distinct from SHAPE-MaP, where modifications are read as mutations rather
than termination sites.
