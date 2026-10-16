# nano-CUT&Tag — source notes

Defining paper: Bartosovic et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01535-4](https://doi.org/10.1038/s41587-022-01535-4).

- 🟢 PMC full text, Methods and Supplementary Table 1 were read.
- 🟢 Secondary nanobody–Tn5 fusions bind mouse or rabbit primary antibodies and carry distinct MeA barcode oligos.
- 🟢 Targeted MeA tagmentation is followed by 10x Chromium linear amplification, a second standard-Tn5 tagmentation with MeB/Me-Rev, and PCR.
- 🟢 Supplementary Table 1 prints the exact representative BcdA top oligo as
  S5–`TCCACGC`–`TATAGCCT`–`GCGATCGAGGACGGC`–ME, plus MeB and
  5′-phosphorylated Tn5_Rev. Only the terminal 19-nt ME duplex is loaded on Tn5.
- 🟢 **Published sequence correction:** The author protocols.io v2 PDF materials
  table, extracted lines 123–138, gives BcdA/B/C as S5–`T`–barcode–connector–ME,
  omitting `CCACGC`. Its correction at lines 271–293 explicitly calls the earlier
  strings wrong and restores S5–`CCACGC`–barcode–connector–ME; published
  Supplementary Table 1 rows 3–14 agrees. The model uses the corrected BcdA
  oligo. 🔴 Whether any earlier experiment used the erroneous strings is unknown.
- 🟢 The Methods state the molecular order: deterministic MeA tagmentation, 10x
  droplet linear amplification and cell barcoding, MeB-only second tagmentation,
  then 10x scATAC PCR.
- 🟢 The published run is 36-8-48-36 (R1-I1-I2-R2), with exact custom R1 and I2
  primers and standard R2/I1.
- 🟢 The page lists the source's full BcdA/B/C × 0N–3N adapter family. The
  public reaction graph draws the representative BcdA_0N branch; the other
  members differ only in the modality barcode and variable sequencing stagger.
- 🟢 The linked 10x Chromium Next GEM Single Cell ATAC v1.1 guide CG000209 Rev G, Appendix A1 (lines 3729–3829 in the extracted copy), prints the gel-bead P5–N16–partial Read1N oligo and SI-PCR primer scaffolds. The cell barcode and sample index are variable windows; strand-resolved PCR binding geometry remains to reconstruct.
  Their molecular position is retained, but their bases and exact PCR landing geometry
  are not reconstructed.
- 🟢 The author protocol's optional-barcode citation to Amini et al. 2014
  (doi:10.1038/ng.3119) leads to its original Supplementary Table 4. Rows 9–11
  print the same P5 adapter scaffold and first three eight-base barcodes as the
  nano-CT BcdA/B/C rows. The workbook's correction notice says its first five
  rows were missing from the initial 2014 online supplement and were restored
  on 20 May 2015. This predecessor supplies optional barcode candidates; the
  nano-CT paper's own table remains the authority for its named BcdA/B/C oligos.
  Two independent readers checked the Amini original and corrected table in
  the source-review checklist.
