# CUT&RUN

## 1. What it is

CUT&RUN maps a chromatin protein by antibody-tethering Protein A–MNase in intact permeabilized nuclei, activating cleavage with calcium, and sequencing the target-proximal fragments released into solution.

## 2. Sources and evidence

- 🟢 Skene PJ, Henikoff S. *An efficient targeted nuclease strategy for high-resolution mapping of DNA binding sites.* eLife 6:e21856 (2017). [doi:10.7554/eLife.21856](https://doi.org/10.7554/eLife.21856).
- 🟢 The defining protocol immobilizes nuclei on concanavalin-A beads, binds antibody and pA–MNase, activates with 2 mM calcium at 0 °C, stops with EDTA/EGTA, and recovers soluble two-cut fragments.
- 🟢 Extracted fragments entered the KAPA DNA polymerase library-preparation workflow without size selection and were sequenced for 25 cycles paired-end on HiSeq 2500.
- 🟢 The contemporary KAPA kit supplied end repair, dA-tailing, ligation and amplification chemistry but required adapters to be supplied separately.
- 🟢 The CUT&RUN Methods say that the sequencing libraries were prepared “as described” by Kasinathan et al. (2014) and Henikoff et al. (2011), and elsewhere state that the same Illumina library-preparation procedure was used for the direct comparison with Kasinathan et al.
- 🟢 Kasinathan et al., *Nature Methods* 11:203–209 (2014), [doi:10.1038/nmeth.2766](https://doi.org/10.1038/nmeth.2766), Supplementary Figure 1c explicitly labels the ligation reagent “Illumina Tru-seq adapters” and shows end repair, 3′-A addition, adapter ligation, PCR amplification and paired-end HiSeq 2000 sequencing.
- 🟢 The public Illumina adapter-sequence guide supplies the canonical TruSeq adapter and sequencing-primer sequences used by the shared repository model.
- 🔴 Neither the CUT&RUN paper nor the cited Kasinathan supplement identifies the particular eight-base TruSeq sample index selected for an individual library. The model therefore preserves the exact single-index position as a placeholder rather than assigning a barcode.

## 3. Molecular path

CUT&RUN first makes unadapted, antibody-selected DNA fragments. Only after fragment release and extraction are ends repaired, dA-tailed, ligated to the cited TruSeq adapters and amplified. This separation is the defining contrast with CUT&Tag.
