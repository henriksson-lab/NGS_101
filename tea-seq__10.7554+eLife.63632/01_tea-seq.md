# TEA-seq — source notes

Defining paper: Swanson et al., *eLife* (2021), [doi:10.7554/eLife.63632](https://doi.org/10.7554/eLife.63632).

- 🟢 The manuscript Methods and Supplementary Table 5 were read. The TEA-seq
  Methods explicitly use that table's `SI-PCR-Oligo` and `ADT-i7` primers
  for the ADT library; the latter has a variable eight-base index and
  terminal phosphorothioate bond (`...CC*A`).
- 🟢 TEA-seq jointly measures transcripts, antibody-derived tags and chromatin accessibility with the 10x Single Cell Multiome platform.
- 🟢 RNA, ADT and ATAC outputs are amplified separately while retaining a shared cell barcode.
- 🟡 The Multiome gel-bead, TSO, adapter and TotalSeq-A antibody-barcode
  families are vendor reagents whose full ordered bases are not printed in
  this paper. They remain typed placeholders, distinct from the published
  ADT PCR primers listed on the page.
