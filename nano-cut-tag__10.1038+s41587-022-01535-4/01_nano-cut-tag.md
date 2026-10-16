# nano-CUT&Tag — source notes

Defining paper: Bartosovic et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01535-4](https://doi.org/10.1038/s41587-022-01535-4).

- 🟢 PMC full text, Methods and Supplementary Table 1 were read.
- 🟢 Secondary nanobody–Tn5 fusions bind mouse or rabbit primary antibodies and carry distinct MeA barcode oligos.
- 🟢 Targeted MeA tagmentation is followed by 10x Chromium linear amplification, a second standard-Tn5 tagmentation with MeB/Me-Rev, and PCR.
- 🟢 Supplementary Table 1 prints the exact representative BcdB top oligo as
  S5–`TCCACGC`–`ATAGAGGC`–`GCGATCGAGGACGGC`–ME, plus MeB and
  5′-phosphorylated Tn5_Rev. Only the terminal 19-nt ME duplex is loaded on Tn5.
- 🟢 **Published sequence correction:** The author protocols.io v2 PDF materials
  table, extracted lines 123–138, gives BcdA/B/C as S5–`T`–barcode–connector–ME,
  omitting `CCACGC`. Its correction at lines 271–293 explicitly calls the earlier
  strings wrong and restores S5–`CCACGC`–barcode–connector–ME; published
  Supplementary Table 1 rows 3–14 agrees. The model uses the corrected BcdB
  oligo. 🔴 Whether any earlier experiment used the erroneous strings is unknown.
- 🟢 The Methods state the molecular order: deterministic MeA tagmentation, 10x
  droplet linear amplification and cell barcoding, MeB-only second tagmentation,
  then 10x scATAC PCR.
- 🟢 The published run is 36-8-48-36 (R1-I1-I2-R2), with exact custom R1 and I2
  primers and standard R2/I1.
- 🟢 The page lists the source's full BcdA/B/C × 0N–3N adapter family. The
  public reaction graph draws the representative BcdB_0N/H3K27ac branch; the other
  members differ only in the modality barcode and variable sequencing stagger.
- 🟢 The linked 10x Chromium Next GEM Single Cell ATAC v1.1 guide CG000209 Rev G, Appendix A1 (lines 3729–3829 in the extracted copy), prints the gel-bead P5–N16–partial Read1N oligo and SI-PCR primer scaffolds. The cell barcode and sample index are variable windows. The model now locates both SI-PCR sites: the forward primer binds the P5 end, and the reverse primer binds the 15-base partial Read2N arm while carrying P7 and the sample index as a 5′ tail. The index placeholder does not assign an experimental plate member.
- 🟢 The author protocol's optional-barcode citation to Amini et al. 2014
  (doi:10.1038/ng.3119) leads to its original Supplementary Table 4. Rows 9–11
  print the same P5 adapter scaffold and first three eight-base barcodes as the
  nano-CT BcdA/B/C rows. The workbook's correction notice says its first five
  rows were missing from the initial 2014 online supplement and were restored
  on 20 May 2015. This predecessor supplies optional barcode candidates; the
  nano-CT paper's own table remains the authority for its named BcdA/B/C oligos.
  Two independent readers checked the Amini original and corrected table in
  the source-review checklist.

- 🟢 **Conflicting fusion-deposit numbers:** Published Supplementary Table 4
  lists anti-rabbit-Tn5 as Addgene #183637 and anti-mouse-Tn5 as #183638.
  The original public [#183637 deposit](https://www.addgene.org/183637/) names
  psfaMs-Tn5 and explicitly says anti-mouse;
  [#183638](https://www.addgene.org/183638/) names psfaRb-Tn5 and says anti-rabbit.
  Both cite this paper. The graph retains the named species affinities rather
  than assigning the conflicting deposit IDs to historical samples. Current
  deposit identities do not independently prove the physical reagent used.
- 🟢 The author's linked `debarcode.py` detects the eight-base modality barcode
  immediately before the reverse-complement outer-connector/S5 pattern and
  extracts 16 cell-barcode bases after it. Its barcode choices are detected
  from input reads or supplied by the caller; its default of three candidates
  is not an experimental barcode-to-antibody assignment. The molecular model
  uses BcdB/H3K27ac from the author study configuration as its representative branch.

- 🟢 **Nanobody provenance checked in the original:** Pleiner et al. 2018
  ([doi:10.1083/jcb.201709115](https://doi.org/10.1083/jcb.201709115)),
  Supplementary Table S1, identifies TP1170 as an anti-mouse κ-chain/Fab
  binder and TP897 as an anti-rabbit Fc binder. The table also records
  cross-reactivity, so these names are reported target affinities rather
  than claims of absolute species exclusivity. Nano-CT's own published
  supplement supplies the full fusion sequences. This predecessor is
  antibody provenance; its dye-labeling and expression procedures are not
  transferred into the molecular library construction.
- 🟢 **Optional accessibility branch:** Nano-CT Methods “nano-CT with
  ATAC-seq” cites Corces et al. 2017
  ([doi:10.1038/nmeth.4396](https://doi.org/10.1038/nmeth.4396)). The original
  Omni-ATAC Methods and Supplementary Protocol 1 describe the accessibility
  preparation and transposition front end. Nano-CT explicitly changes the
  loaded adapter to barcoded MeA/Me-Rev and subsequently follows its own
  Chromium/MeB library construction. Omni-ATAC's ordinary Ad1/Ad2 PCR
  continuation therefore does not define the Nano-CT endpoint. The current
  representative graph covers antibody-targeted molecules; the optional
  accessibility branch remains a declared modeling limit.

- 🟢 **Study barcode assignments recovered from the author configuration:**
  The article-linked repository README points to
  [config/config.yaml](https://github.com/mardzix/bcd_nano_CUTnTag/blob/7657a3a9f15e798b0a8d2ba01ab7f3006a238d21/config/config.yaml).
  In its named study samples, BcdA `TATAGCCT` labels ATAC, BcdB `ATAGAGGC`
  labels H3K27ac, and BcdC `CCTATCCT` labels H3K27me3. The linked
  [nanoscope instructions](https://github.com/bartosovic-lab/nanoscope)
  corroborate this mapping. These assignments belong to those configurations;
  barcode sequences do not have inherent antibody affinities. The earlier
  repository drawing used BcdA as a generic antibody-targeted example. It
  now draws BcdB/H3K27ac to agree with the source-defined study branch. This
  was a repository modeling correction, not a publication error.
