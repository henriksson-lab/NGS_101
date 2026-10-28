# FAIRE-seq

Giresi and Lieb (2009), Isolation of active regulatory elements from eukaryotic chromatin using FAIRE (Formaldehyde Assisted Isolation of Regulatory Elements). DOI: https://doi.org/10.1016/j.ymeth.2009.03.003

This entry models the paper's Illumina GAII sequencing branch of FAIRE, with the human cultured-cell fixation branch. Microarray LM-PCR, qPCR, and tissue-specific alternatives remain separate source branches.

## Source read and identity

🟢 Associated accepted author manuscript PMC2710428, Methods: FAIRE procedure, Crosslinking, Cell lysis, Sonication, Phenol/Chloroform extraction and DNA precipitation (text lines 41–71); Detection by DNA microarray (87–91, branch distinction); High-throughput sequencing (99–103). Manifest identifies an author manuscript, not a publisher version-of-record comparison.

## Preparation

1. 🟢 Human cultured cells: 1% formaldehyde, 5 min at room temperature; quench with 125 mM glycine, 5 min, and wash. Yeast's 30 min fixation is not silently used for human cells.
2. 🟢 Lyse with the stated detergent/glass-bead protocol and sonicate. The source checks 100–1000 bp bulk sonication fragments; the recovered FAIRE fraction is typically 75–200 bp. These are different points in the workflow, not contradictory insert measurements.
3. 🟢 Clear cellular debris by centrifugation. Perform phenol/chloroform extraction, recover the aqueous phase, repeat extraction, then chloroform extraction. Crosslinked protein-associated DNA preferentially partitions away from the aqueous FAIRE fraction. This is an enrichment, not an exact binary classifier of nucleosomes or sequence-specific occupancy.
4. 🟢 Sodium acetate/glycogen/ethanol precipitation and 70% ethanol wash. RNase A acts for 1 h at 37°C after resuspension. RNA digestion leaves debris in the same preparation until the subsequent spin-column or extraction/precipitation cleanup; it is not physically discarded at digestion.
5. 🟢 The source explicitly says the old overnight 65°C crosslink-reversal step can be omitted without detected change in downstream enzymatic efficiency. It is not added to this selected branch. The reversed-crosslink input control in the microarray branch is also not imposed here.
6. 🟢 Use 100–200 ng recovered FAIRE DNA for GAII library preparation: blunt with Epicentre ER0720, add A overhang with KL06041K, ligate double-stranded T-overhang Illumina adapters #1000521. Select 125–275 bp ligation products on 2% agarose; PCR with PfuUltraII and Illumina primers #1000537/#1000538.
7. 🟢 GAII, 36 bp genomic end reads. The paper's computational extension to 200 bp estimates signal density; it does not create 200 sequenced bases and is not a molecular reaction.

## Oligo and model boundaries

The scraper's sequence family belongs to the explicitly separate microarray LM-PCR linkers (5′-GCGGTGACCCGGGAGATCTGAATTC-3′ and 5′-GAATTCAGATC-3′). Those are not assigned to the Illumina sequencing branch. Historical Illumina product identities are source-positive, but their exact original polymer sequences, complete duplex geometry and run primers are unassembled in this pass. Source-known A/T terminal identities can be drawn without inventing adapter bodies.

🟡 The graph is a source-qualified partial molecular model and keeps protein-associated fractions and RNA debris separate from DNA strands. It does not demonstrate quantitative enrichment, absence of crosslinks, or a complete sequencing library sequence. No sample index or UMI is invented; unassigned run architecture is not an absence claim.

The article mentions one-page cell/tissue supplementary protocols, but their original attachment identity has not been resolved by the current JATS bundle. This is a real linked-source target; a source URL was not fabricated. Full publisher article, those original supplements and historical vendor library documentation remain deferred within the capped pass.

## Lineage

The article explicitly cites the original FAIRE report by Giresi et al. (2007), DOI 10.1101/gr.5533506. Component scope is formaldehyde/phenol fractionation for nucleosome-depleted DNA; the present 2009 Illumina preparation is not inherited merely from the assay name. LM-PCR microarray references establish a different branch and do not establish NGS adapter inheritance.
