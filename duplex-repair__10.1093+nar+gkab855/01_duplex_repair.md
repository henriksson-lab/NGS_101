# Duplex-Repair

**Source.** Xiong et al., *Nucleic Acids Research* (2022), [doi:10.1093/nar/gkab855](https://doi.org/10.1093/nar/gkab855), PMCID PMC8755016; Methods and Supplementary Figures S2–S10.

- 🟢 Step 1 excises damaged bases and removes single-stranded overhangs with EndoIV, Fpg, UDG, T4 PDG, EndoVIII and ExoVII (37°C, 30 min).
- 🟢 Step 2 adds T4 PNK and T4 DNA polymerase with ATP and all four dNTPs to blunt ends and fill only remaining short gaps and overhangs (37°C, 30 min).
- 🟢 Step 3 seals nicks with HiFi Taq DNA ligase while the cycler ramps from 35°C to 65°C over 45 min.
- 🟢 Step 4 uses Klenow fragment (exo−) and Taq DNA polymerase with 0.2 mM dATP as the only nucleotide (room temperature 30 min, then 65°C 30 min). Restricting the nucleotide supply is what prevents templated extension during dA-tailing.
- 🟢 T4 DNA ligase, 5′-deadenylase and 10% PEG 8000 attach custom dual-index duplex-UMI adapters; the products are PCR-amplified and were sequenced paired-end on HiSeq 2500 or HiSeq X.
- 🟢 The paper describes the adapters by role but does not publish their base sequences. The schematic therefore stops at labelled adapter roles and does not invent sequencing-primer sites.
- 🟢 Supplementary PDF, “Quantification of library conversion efficiency by ddPCR,” prints adapter-flank **assay** primers `CACTCTTTCCCTACACGACG` and `AGTTCAGACGTGTGCTCTTC`. These amplify doubly ligated products for QC; they do not disclose the complete custom dual-index duplex-UMI adapter molecules or the library-amplification primers. Supplementary Table S1's modified oligos are synthetic damage/repair substrates for capillary electrophoresis, not the ligated library adapters.

The central distinction from ordinary combined end-repair/dA-tailing is reaction order: damaged ends are cleaned, limited fill-in is completed, and internal nicks are sealed before polymerases see the dATP-only tailing reaction.
