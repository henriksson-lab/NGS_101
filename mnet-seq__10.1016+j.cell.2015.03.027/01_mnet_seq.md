# mNET-seq

[Nojima et al., Cell (2015)](https://doi.org/10.1016/j.cell.2015.03.027) profiles native Pol II-associated RNA, including cotranscriptional processing intermediates. It is distinct from Mayer's human NET-seq library and does not inherit its linker, UMI, RT primer or circularization route.

## Selected 2015 chemistry

🟢 The defining main article (JATS lines 83, 199–211) and Extended Experimental Procedures (Document S1, current text lines 2497–2510) specify:

1. Wash isolated chromatin (from ~1.6 × 10^8 cells) in 1× MNase buffer, digest with MNase (40 U/µl, 37 °C, 90 s, 1,400 rpm) and stop immediately with 25 mM EGTA.
2. Centrifuge (13,000 rpm, 5 min) and collect the soluble chromatin supernatant; dilute in 9 ml NET-2 and immunoprecipitate the selected Pol II CTD state with antibody-conjugated beads (40 µg antibody, 4 °C, 1 h); wash six times with NET-2 and once with PNKT.
3. Treat the washed immunoprecipitation beads with **T4 PNK, 3′ phosphatase minus**, and ATP. This phosphorylates RNA 5′ ends before extraction while preserving the informative Pol II-protected 3′ hydroxyl.
4. Wash once with NET-2, extract RNA with Trizol, resolve on an 8% acrylamide/7 M urea gel, elute the 35–100-nt interval (1 M NaOAc, 1 mM EDTA) and ethanol-precipitate.
5. Prepare libraries using the named TruSeq Small RNA kit and sequence paired reads on HiSeq 2000/2500.

The former repository model extracted RNA before a generic kinase reaction. It also supplied an eight-base index through a shared helper. These were repository modelling errors, not publication errors. The new graph owns its route and keeps unassigned molecular details as explicit roles.

## Positive kit-family information and version boundaries

🟢 The recovered manufacturer guide is **15004197 Rev. C, March 2011**, mirrored by Boston University. It describes sequential 3′/5′ adapter ligation, truncated T4 RNA ligase 2, kit T4 RNA ligase with ATP, SuperScript II, RT/PCR and a six-base RPI index family. Its title identifies its own edition, not the exact guide or kit lot used in 2015.

🟢 The current official manufacturer sequence table, document **1000000002694 v22**, prints these family references:

- RA5 (RNA): `GUUCAGAGUUCUACAGUCCGACGAUC`.
- RA3 (DNA): `TGGAATTCTCGGGTGCCAAGG`.
- RTP: `GCCTTGGCACCCGAGAATTCCA`.
- RP1: `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`.
- RPI: `CAAGCAGAAGACGGCATACGAGAT[6 bases]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`.

These are source-positive manufacturer-family data, displayed separately from actual 2015 reagent assignment. They do not identify a selected sample index, exact historical modifications or the physical sequencing primers for that experiment. The declaration leaves the ordered 2015 endpoint unassembled; it does not call these sequences unpublished.

🟢 The separately recovered [Nojima et al. 2016 detailed protocol](https://doi.org/10.1038/nprot.2016.012), PMC5152764, is a successor. Its author PDF steps 69–95 names deletion-mutant ligase for RA3, wild-type ligase for RA5, SuperScript III, RP1/RPI, a partial RT transfer with a retained remainder and final gel selection of 150–230-bp product excluding 120–125-bp background. Those enzyme, transfer and gel details are not automatically assigned to the 2015 implementation. The difference from the 2011 manufacturer guide is retained as source-version context, not labelled a publication error.

## Published source difference

The defining article (JATS line 201) and original Document S1 (text line 2502) specify **0.1% Triton X-100** in the PNKT bead wash. The repeated mNET-seq extraction paragraph in the deposited GSE60358 metadata instead specifies **0.05%**. Both original values are retained in a public schematic notice. The graph preserves the physical wash and phosphorylation order without selecting this concentration; the actual experimental value remains unresolved. No typographical error is assumed.

The same sources also differ in the MNase-buffer chromatin wash: **100 µl** in the article (JATS line 201) and Document S1 (text line 2498), **1 ml** in the GEO paragraph repeated for 18 samples. This is also shown on the schematic; the graph selects no wash volume. The GEO paragraph renders several other microlitre/microgram amounts as ml/mg (500 ml PNKT, 100 ml PNK mix, 40 mg antibody, 40 u/mL MNase), so those unit renderings are not interpreted as dose changes. The PNK amount printed in the article (0.05 U/ml) is transcribed but not assigned as a model concentration.

## Readout and molecule fates

🟢 The defining processing section states that **Read 2's first base** marks the RNA 3′ end; Read 1 supplies strand direction. Exact acquisition cycles and index reads remain unassigned in this model. Figure S1's explicit 51-base statement describes the separate ChrRNA-seq panel and is not silently transferred to mNET-seq.

The insoluble chromatin pellet, unbound immunoprecipitation material and excluded RNA gel intervals are physical separations. Adapter trimming, exclusion of short/unpaired reads, mapping and keeping the informative nucleotide are computational operations. They do not imply molecular discard. Signals from processing intermediates are not all equated to polymerase positions.

## Supplement and linked-source scope

The separate original Table S1 is a six-factor CPA CLIP peak-position result table, not an oligo table; both readers inspected its complete body. Document S1 contains the original article, extended methods and Figures S1–S7. Its siRNAs, validation RT-PCR primers and ChIP primers have distinct experimental roles and are not assigned as sequencing-library oligos.

Original GSE60358 metadata contains 28 samples, including 18 mNET-seq extraction descriptions and ten comparator/control descriptions; its instrument fields report HiSeq 2000 and HiSeq 2500. RAW archives and processed bigWig files were not downloaded. Identified linked predecessors, remaining figure contents, experiment-level metadata and software/vendor scope remain explicitly pending in the source review. No source absence is accepted from an unassigned study revision or sample member.
