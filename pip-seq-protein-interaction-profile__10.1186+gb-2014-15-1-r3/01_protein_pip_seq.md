# PIP-seq (protein interaction profile sequencing)

**Defining source:** Silverman *et al.*, "RNase-mediated protein footprint sequencing reveals
protein-binding sites throughout the human transcriptome", *Genome Biology* 15:R3 (2014),
doi:10.1186/gb-2014-15-1-r3. Read: the open PMC4053792 JATS text, specifically the Results
overview (lines 61–71) with the Figure 1A legend (line 65), and the Methods sections Cell lines,
Cross-linking experiments and PIP-seq library preparation (lines 155–171). Also read:
Additional file 1 (`S1.xlsx`, Supplemental Table 1) and, for the cited predecessor, Zheng
*et al.* 2010 Text S1 ("dsRNA-seq library preparation", lines 58–59). Last checked 2026-10-17.

## Evidence used

- 🟢 **Cross-linking.** HeLa cells at 90% confluence are treated in one of three ways:
  1% formaldehyde for 10 min at room temperature, quenched with 125 mM glycine for 5 min;
  UV-C at 254 nm, 400 mJ/cm² (Stratalinker 2400); or not cross-linked.
- 🟢 **Lysis.** Cells are lysed in RIP buffer (25 mM Tris pH 7.4, 150 mM KCl, 5 mM EDTA,
  0.5% NP40, 10 µM DTT, protease inhibitors; 850 µL per 10 million cells), ground by hand,
  and treated with RNase-free DNase (Qiagen).
- 🟢 **Footprint samples.** The lysate is digested either with RNaseONE (Promega;
  100 U/mL, 200 µg/mL BSA, 1× RNaseONE buffer, 1 h at room temperature) or with RNaseV1
  (Ambion; 2.5 U/mL in 1× RNA structure buffer, 1 h at 37 °C). Proteins are then removed
  with 1% SDS and 0.1 mg/mL proteinase K for 15 min. Cross-links are reversed at 65 °C for
  2 h (formaldehyde) or by a second treatment with 1% SDS and 1 mg/mL proteinase K for
  30 min (UV).
- 🟢 **RNase digestion controls.** The lysate is first treated with 1% SDS and 0.1 mg/mL
  proteinase K and ethanol-precipitated. The same RNase digestions follow, and otherwise
  "identical treatments". Each replicate therefore yields four libraries:
  ss/ds RNase × footprint/control.
- 🟢 RNA is isolated with the Qiagen miRNeasy kit. A strand-specific library is then prepared
  "as previously described [29,30]", that is, by Li *et al.* 2012 (*Cell Reports*,
  doi:10.1016/j.celrep.2011.10.002) and Zheng *et al.* 2010 (*PLoS Genetics*, PMC2947979).
  The one change is the addition of Illumina DSN normalization: 100 ng of the final library
  is denatured at 95 °C, re-annealed for 5 h at 68 °C, and treated with 2 µL of DSN
  (1 U/µL) to remove abundant duplexes that re-anneal quickly.
- 🟢 Libraries were sequenced on an Illumina HiSeq2000 with the standard 50-bp single-read
  protocol. The paper says nothing about indexing or multiplexing, so the model has no sample
  index. Earlier versions of this page declared one without a source.
- 🟢 The cited Zheng predecessor (Text S1) describes the following route: RNase One, then
  phenol:chloroform extraction, Applied Biosystems fragmentation, T4 PNK end treatment, and
  the Illumina Small RNA Sample Prep v1.5 kit. Whether PIP-seq applied exactly this route,
  or the Li 2012 variant, is not stated. The model therefore marks library construction as a
  warned "magic happens" boundary rather than reconstructing ligation, RT or PCR.
- 🔴 The defining article and the predecessor text that was checked do not print the complete
  library adapter, RT, PCR or sequencing-primer oligos. The oligo scraper's hits in the article
  are the T7 and SNP probe oligos used for the UV-cross-linking validation assay; its hits in
  the Zheng supplement are qPCR and genotyping primers. Neither set belongs to the library.

## Published discrepancy

- The Results (line 71) report "approximately 31 to 60 million raw reads per library", but the
  raw-read column of Additional file 1 runs from 24,604,010 reads (UV dsRNase footprint,
  replicate 2) to 103,702,805 reads (formaldehyde dsRNase footprint, replicate 2). The page
  shows this discrepancy and does not model read depth.

## Open questions

- Li *et al.* 2012 (*Cell Reports*) is not in PMC. The manual URL is
  https://doi.org/10.1016/j.celrep.2011.10.002. Its strand-specific library chemistry is the
  most direct remaining source for the adapters.
- The historical Illumina Small RNA Sample Prep v1.5 guide (adapter and primer sequences).
- The detailed protocol, Current Protocols in Molecular Biology doi:10.1002/cpmb.21, which is
  listed in `ref/MANIFEST.md`. It was not read in this pass.
- The PMC-gated Additional files 2–9 (figure supplements).

## Naming

This is unrelated to the particle-templated single-cell RNA-seq method also called
PIP-seq. The catalogue and URL qualify this page as “protein interaction profile.”
