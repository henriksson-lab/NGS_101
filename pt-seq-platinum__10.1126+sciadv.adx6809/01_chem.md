# Pt-seq: platinum-DNA adduct mapping

Shao et al. (2025), **Pt-seq unveils the genomic binding pattern of platinum-based drugs**, ScienceAdvances11:eadx6809. [Published original](https://doi.org/10.1126/sciadv.adx6809). This is platinum Pt-seq, not the phosphorothioate PT-seq homonym.

🟢 Inspected three originals: published main PMC12571083 Results111–113/Methods193–205; published supplement FigS1, originalpage2 rendered; TableS1 workbook, originalcellsA5:C11 inspected directly. The cap is closed. Other workbook sheets/model/qPCR/target-validation oligos remain unreviewed for adoption. Mandatory primer scraper yielded38 candidates in2 sequence-positive files; relevant TableS1 oligos were checked in the original, not adopted from the scrape.

🟢 The library suppresses background using a restriction-site adapter and protected ends. After Pt antibody enrichment, primer extension can stall at a Pt adduct; non-stalled background can instead form the duplex restriction site. SbfI digestion and exonuclease processing are genuine physical preparation steps beyond a generic Damage-seq or cisplatin-seq analysis. The graph does not equate predicted stall specificity, antibody specificity or phosphorothioate resistance with perfect enrichment.

Source ordering oligos, verbatim sequence and separately printed modifications:

- 1st-Y-Adaptor-U: `T*T*T*T*T*T*TTTCCTGCAGGGCTCTTCCGATC*T` — Phosphorothioate.
- 1st-Y-Adaptor-B: `GATCGGAAGAGCACACGTCTGAACTC*C*A*G*T*C*A` — 3`NH2 C7/5`P/Phosphorothioate.
- primer_extension: `GACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.
- 2nd-Y-Adaptor-U: `ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNN` — 3`NH2 C7.
- 2nd-Y-Adaptor-B: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` — 3`NH2 C7/5`P.
- Universal Primer: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC`.
- Index Primer: `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.

The inventory retains internal phosphorothioate bonds, nucleotides after the final star, five-prime phosphate and three-prime NH2 C7 blocks as distinct properties. Only the first-adapter12bp core is paired by Scene; unpaired arms and terminal modifications remain in the separate full oligo inventories. This is not a drawing of complete Y-adapter ligation geometry. The second-adapter6N window has no adopted UMI function. The Index Primer6N sample window has structured sample-index meaning, with the actual index choice unassigned. The main paper's UMI processing belongs to separate HGPRT validation sequencing and is not imposed on genomic Pt-seq.

🟢 Selected source preparation: extracted500–1000ng genomic DNA →CovarisS220 median300bp fragmentation →physical addition of1ng independently prepared model sequences →NEBE7546 endrepair/dA-tail component →NEBE7595 first-adapterligation →denature/ice →save nonimmunoprecipitatedinput at−20°C separately →antiPtDNAantibody1µg overnight4°C →BSA/PBS-blocked DynabeadsProteinG2h4°C →three IP washes →proteinaseK56°C2h elution and source purification.

The model-sequence stock is an explicit input boundary: species, Pt-adduct composition and native versus model fraction remain unassigned. Actual model preparation/orthogonal MALDI/ICP-MS experiments are not silently inserted into production. The selected genomic extraction boundary does not assign a drug, cell type or sequencing platform to a particular specimen without evidence.

🟢 RecJf plus exonucleaseI digestion →Q5 primer extension →exonucleaseI primer digestion →denature/ice →second adapter/T4DNA ligase/PEG8000 at16°Covernight →Q5two-cycle preamplification →1.2×AMPure →FIRSTSbfI37°C1h →1.2×AMPure →Q5six-to-eight-cycle amplification →SECONDSbfI →final1.2×AMPure. The second digest's incubation conditions are not restated in Methods and remain unassigned. No cleanup is invented between exonucleases, extension, second ligation or the first two cycles; digest fragments and surviving substrates remain until source physical separation. No PCR is added after the second SbfI digest.

🟢 Separate Qubit4/Equalbit and Agilent4150TapeStation measurements precede NovaSeq6000/MGISEQ2000/MGISEQT7 acquisition, PE150. Specimen-to-run assignments and historical sequencing/ligation/amplicon geometry remain unassembled. Source post-acquisition SbfI terminal-read filtering does not physically remove DNA during preparation.

🟡 Private graph and ordering inventory only. No perfect cleavage/protection, exact adduct-specific stop position, whole library sequence, final read/index landing or accepted absence claim is asserted. Historical E7546/E7595 documents and cited IP/library-dedication predecessor sources remain to inspect. Source Results113 credits Paired-seq2019 (reference17, DOI10.1038/s41594-019-0323-x) and PairedTag2021 (reference18, DOI10.1038/s41592-021-01060-3) for the restriction-site library-dedication inspiration. These exact reference identities were verified in the main JATS bibliography. This supports component design-inspiration links, without establishing equal molecular inputs, ligation geometry or whole-library inheritance; those predecessor originals remain unreviewed in this pass.
