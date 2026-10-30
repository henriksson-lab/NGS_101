# RoCK and ROI: selected published eGFP preparation

Moro et al., *RoCK and ROI: single-cell transcriptomics with multiplexed enrichment of selected transcripts and region-specific sequencing*, DOI10.1038/s41467-025-66248-z(2025). This is one method entry with complementary capture/priming components, not separate target-panel or assay-count entries.

🟢 Main288–302 changes the TSO portion of BD Rhapsody beads using T4DNApolymerase and5P splints. SelectedTable11eGFPsplint is `/5Phos/CATGGTCCTGCTGGAGTTCGTGACCCATACCTACTACGCATA`;5P protectiveA18 preserves dT oligos; lambdaexonuclease removes complementary strands. Local nativeTSO:splint pairing is computed from source-specified suffix complementarity, with the whole upstream bead construct unassembled. This differs from simply inheriting DART-seq ligation chemistry.

🟢 Full vial is processed as four500ul bead vials. Each gets source preincubation WITHOUT T4polymerase5min37°C, actual nucleotide/splint/protection additions,6ul5U/ulT4polymerase10minRT→75°C10min→washes→lambda37°C30min→75°C10min→washes/storage. The source's stated buffer/oligo concentration conventions are preserved, not recomputed into assumed corrections. Experiment allocation from these stocks is unassigned; the model neither pools all four factory tubes nor consumes the full stock silently. FluorescentQC20ul/scaled reactions remain a separate unassembled analytical branch.

🟢 ROIeGFP primer is `TCAGACGTGTGCTCTTCCGATCTGGCAACTACAAG`; source adds4ul from10ul dilution of1ul100uM primer, leaving6ul stock. This is second-strand priming alongside standard randomers after source RNA capture/cDNA synthesis. Tprimer is `ACAGGAAACTCATGGTGCGT`. Whole custom indexing and sequencing primers are model declarations derived from canonicalP5 where applicable and checked against originalTable11, including all fixed tails. Additional5P TSOprotective oligo is `CATACCTACTACGCATA`; this inventory does not assert a second experimental addition at100%modification.

🟢 WTA and TSO indexingPCRs use separate aliquots of the same primary library, same8bp sample index and different forward primers. Selected library aliquots form the actual1:1concentrationpool; stock remainders remain unresolved. CustomRead1primer and physical3%PhiXDNA stock are added for sequencing. Selected runR1:60/R2:62; leukemia/HeLaR2:180 branches are qualified alternatives. No whole primer landing or sealed-product proof is claimed.

🟡 Exact cited protocols.io v2 and BDdocs210966/23-21711-00 remain ordinary unreviewed, despite automatic source retrieval. The model does not replace them with modern kits or a newer revision. Whole bead/cDNA/amplicon/read structures and factory allocations are unassembled, not source-absent.

## Selected sources and review scope

- s41467-025-66248-z_PMC12695968.xml: Published main288–352 selected bead modification, ROI primer/library generation and sequencing; main100–118 preparative rationale.
- s41467-025-66248-z_supp_41467_2025_66248_MOESM1_ESM.pdf: PublishedSI1Table11 originalp41 rendered; selected whole eGFPsplint/ROI/T/index/customR1/protective oligos and5Pcolumn verified. Other SIfigures/experiments not fully inspected.

This capped pass is not a complete source, linked-information, chemistry or absence audit.
