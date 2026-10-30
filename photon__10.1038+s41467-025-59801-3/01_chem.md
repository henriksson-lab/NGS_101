# PHOTON: selected cell preparation

Rajachandran et al. (2025), **Subcellular level spatial transcriptomics with PHOTON**, DOI: 10.1038/s41467-025-59801-3.

## Scope and sources

🟢 Published article PMC12078482, Fig.1 and Methods192–222; published Supplementary Information MOESM1, Table1, original PDF p6. Source text and the whole original table were inspected. Two originals; source investigation closed. The automatically fetched 2024 preprint is superseded and unused. Mandatory scraper ran, and its pointers were checked against these sources. This is a selected capped pass, not a full linked-source or missing-information audit.

## Distinct preparation

🟢 Fixed/permeabilized cells receive both oligo-dT and random-hexamer RT primers. Each carries a fluorescent photocleavable 5′ linker, four N-POM-caged thymidines and an internal biotin-dT. ROI light treatment uncages bases and exposes the primer-derived cDNA 5′ phosphate. After proteinase-K lysis and column purification, a preannealed splint/P7 handle is ligated, and biotinylated cDNA is captured on streptavidin beads. Continued cDNA extension/template switching occurs after capture. This differs from installing photogated transposase adapters on genomic DNA.

## Oligos and local geometry

🟢 Table1 NPOM-RT common ordering, preserving modifications:
`/5Alex546N//iSpPC/G/iNPOM-dT/C/iNPOM-dT/CG/iNPOM-dT/GGGC/iNPOM-dT/CGGAGATGTGTATAAGAGACAG/iBiodT/`
followed by the literal tail `TT TTT TTT TTT TTT TTT TTT TTV N`. NPOM-RT-Random shares the body and ends in `NN NNN N`. Random priming degeneracy is not a source-declared UMI.

🟢 Splint: `CCGAGCCCACGAGACTCGGACGATCATGGG`.
P7 handle: `CAAGCAGAAGACGGCATACGAGATCCCATGATCGTCCGA`.
TSO: `AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG`.
PCR_TSO: `AAGCAGTGGTATCAACGCAGAGT`.
PCR_P7: `CAAGCAGAAGACGGCATACGAGAT`.
The table's unblocked control RT oligos are not substituted for experimental caged primers.

🟡 The local Scene computes pairing of the two 15-base splint arms to the RT-derived 5′ handle and P7-handle 3′ tail, placing an adjacent nick. This models source-known ligation-ready ends; complete RNA/cDNA, sealed product, uncaging efficiency and recovery remain unassembled.

## Order and physical fates

🟢 Main192–222: .4%PFA10min → PBS washes → .5%Triton15min → PBS → RT30minRT/overnight37°C → ROI405nm → PK55°C30min → NucleoSpin → splinted ligation → streptavidin capture/washes → continued RT/TS → washes → five-cycle Q5 amplification → 2.5µl analytical qPCR allocation → additional production cycles → collect supernatant/.8×AMPure → NexteraXT → .7×SPRI → analytical QC → selected DASH route → PK → .8×SPRI → sequencing.

🟢 Adapter annealing makes 50µl stock; only 10µl enters ligation. The model keeps the 40µl remainder. The 2.5µl qPCR assay is not returned to production; its one-third-saturation signal determines additional cycles. Capture/workup washes remove physical fractions without asserting perfect chemical selection. The second source RT-labelled stage is represented as continued extension/template switching with TSO as acceptor, not a newly supplied RT primer. Continued PCR consumes the existing production reaction and adds no second Phusion/Q5 stock.

🟢 The selected HeLa SG route uses DASH. The mitochondrial RNA experiment instead uses SEQuoia RiboDepletion; that alternative is not merged. Cas9 digestion products remain in the mixture until actual workup. Read filters are computational rather than physical DNA depletion.

## Endpoint and limits

🟢 The source names NexteraXT, NextSeq2000 or NovaSeqXPlus; specific index primers, guide panel, complete kit internals, product/primer/read geometry and cycle masks remain unassigned locally. No predecessor kit constructor is imported to fill them. No sequence/reagent absence is claimed. Linked cited-component/vendor originals remain unreviewed within this capped pass.
