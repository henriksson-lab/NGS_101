# ICICLE-seq — source notes

Defining paper: Swanson et al., *eLife* (2021), [doi:10.7554/eLife.63632](https://doi.org/10.7554/eLife.63632).

- 🟢 The manuscript Methods and Supplementary Table 5 were read.
- 🟢 Supplementary Table 5 supplies the ordered `Poly-A Top-L` with 30 T bases,
  degenerate `VN` and the 19-nt mosaic end; the `/5phos/` Tn5 bottom strand;
  the S7 top strand; SI-P5-22, F BC, ADT-Rev-AMP and SI-PCR primers;
  the TruSeq Read 1 sequencing primer; and the ADT-i7 primer with an eight-base
  variable index and terminal phosphorothioate bond (`...CC*A`). These are
  ordering oligos, separate from the commercial antibody-tag barcode family.
- 🟢 Custom Tn5 complexes carry either a poly(dT) or Nextera R2N overhang; 10x 3′ GEMs barcode ATAC fragments and polyadenylated TotalSeq-A tags.
- 🟢 Separate ATAC and ADT libraries use 28-bp Read 1, 8-bp i7 and 100-bp Read 2; a custom TruSeq Read 1 primer was spiked in.
- 🟡 The custom capture branch is drawn structurally rather than mislabeled as a biological barcode.
- 🟡 The page now lists representative exact custom oligos, while the 10x bead
  oligo and TotalSeq-A antibody-tag members remain vendor-controlled inputs.
