# TruePrime Single Cell WGA

**Source.** Picher et al., *Nature Communications* 7, 13296 (2016), [doi:10.1038/ncomms13296](https://doi.org/10.1038/ncomms13296), PMC5141293, Methods “Single-cell WGA” in `ncomms13296_PMC5141293.xml.txt`. A SYGNIS *TruePrime Single Cell WGA Kit version 2.0 Handbook* is available from a distributor [mirror](https://www.pacificscience.co.th/userfiles2/fckeditor/121912/file/Sygnis-TruePrime-Single-Cell-WGA-Handbook(1).pdf), pp. 10–12.

- 🟢 The paper puts one visually inspected HEK293 cell in 2.5 µl PBS, adds 2.5 µl Buffer L2, incubates 10 min on ice, adds 2.5 µl Buffer N, then adds 42.5 µl reaction mix. The mix is 26.8 µl water, 5 µl reaction buffer, 5 µl dNTPs, 5 µl TthPrimPol (Enzyme 1), and 0.7 µl phi29 polymerase (Enzyme 2). Incubation is 3 or 6 h at 30 °C followed by 10 min at 65 °C.
- 🟢 TthPrimPol synthesizes the short DNA primers de novo, then phi29 extends them by strand displacement. No synthetic random primer is added in this reaction.
- 🟢 The `GTCC`-bearing synthetic template in Figure 3a and the mycoplasma PCR primers in the Methods are assay and cell-culture reagents, respectively; neither is a primer added to the single-cell WGA reaction.
- 🟢 The paper's downstream sequencing uses Covaris fragmentation, NEBNext library preparation, and paired-end HiSeq 2500.
- 🟢 The paper-linked [ENA study SRP085855](https://www.ebi.ac.uk/ena/browser/api/xml/SRP085855) identifies a TruePrime single-cell WGA run and comparator runs; its experiment metadata do not specify the downstream NEBNext kit or adapter.
- 🟡 The handbook version 2.0 differs at lysis: p. 11 says 5 min at room temperature then 1 min at 95 °C, followed by ice; p. 12 gives 3 h at 30 °C with 6 h optional. The paper gives 10 min on ice and does not identify the handbook revision used. The model follows the paper. This is a source-version difference, not evidence that either source is erroneous.
- 🟡 [SYGNIS announced version 2 on 11 July 2016](https://www.eqs-news.com/news/corporate/sygnis-ag-announces-global-launch-of-trueprimetm-single-cell-wga-kit-v2/b6944fd3-d5aa-471d-a00c-5203ce0835df_en), before the paper's publication. Publication order alone does not date the experiments or identify their kit revision.
- 🔴 The paper does not identify the NEBNext kit generation, chosen adapter or index oligos, or sequencing-primer landing sites. The unresolved vendor revision also limits extrapolation from the version 2.0 manual.

The separately modelled TruePrime purified-DNA WGA branch starts from extracted genomic DNA and Buffer D alkaline denaturation; it must not be substituted for intact-cell lysis.
