# CAGEscan

**Defining paper:** Plessy et al., *Nature Methods* (2010), doi:10.1038/nmeth.1470. The open article and supplement were read from PMCID PMC2900199.

## Molecular boundary

🟢 A random-hexamer-bearing RT primer creates a downstream cDNA end. Cap-dependent template switching then installs a tag at the transcript 5′ end, linking a transcription start site to a downstream transcript segment without restriction cleavage.

🟢 The original supplementary primer table publishes eight CAGEscan oligos: the template-switch primer with a 3- or 6-base DNA barcode and terminal RNA `rGrGrG`, the `N6` RT primer, two semi-suppressive PCR primers, two library PCR primers and two custom paired-end sequencing primers. All eight are transcribed in the page's oligo list. Source: `nmeth.1470_41592_2010_BFnmeth1470_MOESM315_ESM.pdf.txt`, lines 5135–5172 (the NIHMS supplement carries the same table).

## Evidence boundary

🟡 The molecular workflow still renders some primer arms by role rather than annealing each published oligo to a base-resolved construct. The historical PCR and custom sequencing oligos establish the cap-end and downstream-end roles; their exact landing sites require a complete historical library reconstruction.
