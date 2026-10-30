# CAT-ATAC

[Published defining article](https://doi.org/10.1016/j.crmeth.2025.101222), *Simultaneous capture of single cell RNA-seq, ATAC-seq, and CRISPR perturbation enables multiomic screens to identify gene regulatory relationships* (2025).

🟢 Two published originals: PMC12859477main122/STARoligotable379–425/Methods486; DocumentS1MethodS1actualp16–20. Originalp16/p19 tables were rendered and read; mandatory whole-bundle scraper pointers verified. Otherguide-panel workbooks and combinedarticlePDF are not added as new chemistry sources. Prepared perturbed cells/nuclei, rather than lentiviral construction or editing treatment, are the graph boundary.

## Capture and GEM programme

🟢 Source nuclei extraction/tagmentation delegates to10xMultiome. Preanneal sourceRTprimer andLNAsplint bylinkedIDTannealing protocol. RTprimer: `/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNNNTTGCTAGGACCGGCCTTAAAGC`. ItsN12 is explicitly aUMI and the finalCS1-complementary part primes gRNA. Splint: `TGACGCTGCC+G+A+C+G+ACAGACGCG/3Phos/`; fiveLNA residues and3′P are preserved, not represented as ordinary unmodifiedDNA chemistry. Exact local15bp RT/splint complementarity is computable; occupied bead/gRNA strands, complete junctions and LNA thermodynamics remain unassembled.

🟢 Source64µlGEMmix:49.5barcodingreagentmix,1.1templateswitcholigo,1.9reducingagent,1µlRTsplintduplex describedas0.5µM,7.5barcodingenzymemix,3µlHi-T4ligase. The0.5µM descriptor is retained without silently deciding stock versus final convention.25°C30min precedes53°C45min, then25°C30min/4°Chold. Main122 states that priorATAC-fragment ligation prevents loss of capture at the higherRTtemperature; Hi-T4remains available after53°C. The source's hypothesis about RT filling ends is not asserted as a complete sideproduct prediction. Co-retained original primers/ligase are not freshly added at each thermal step.

🟢 PostGEMcleanup replaces1.8×beads with3.2×SPRI(160µl)+1.8×isopropanol(90µl), intended to retain single-strandedgRNAcDNA. Actual molecular recovery/classes remain unassigned. SourcestandardATAC/GEX library production is a separate opaque historical10xroute, not imported currentkitconstructors.

## CRISPR library branch

🟢 Take40µlpreamp plus60µlEB; first60µlSPRI then40µlSPRI remove largecDNA/retain selectedsmaller material; elute30µl. Otherpreamp/ATAC/GEX allocations remain source-qualified/opaque; the removed large fraction is not automatically the GEXproduction branch. FeaturePCR uses5µltemplate in100µl:50KapaHiFiHotStartReadyMix,4µl10µMpartialP5(oligo8),4µl10µMR2TSOprimer(oligo5),37water.10cycles98°C20s/58°C5s/72°C5s after98°C45s; final72°C1min. Actual2µlTapeStationQC precedes0.8×–0.9×SPRIcleanup/elute30µl; secondQC amount is not stated. Sourceprimer8 `AATGATACGGCGACCACCGAGA`; primer5 `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGCAGTGGTATCAACGCAGAG` is a PCRprimer, not automatically theGEMtemplateswitcholigo.

🟢 IndexPCR:5µltemplate/100µl, same50/4/4/37 stock layout withpartialP5 and one of sixP7-i7-R2 choices.9cycles98°C20s/60°C30s/72°C20s after98°C45s; final72°C1min. Actual1µlQC,0.7×–0.8×double-sidedSPRI/elute30µl,1µlpostcleanupQC. SourcewholeP7primerfamily consists of canonicalP7 +one8ntindex+TruSeqRead2. Six sourcei7values: TTCGCAGT,CGAGACTA,ACAGCTCA,AAGTGTCG,AGTAAACC,CCGTTTAG. The published next-page sample-sheet values are reverse complements, retained separately; actualindexchoice per study sample remains unassigned here. Full guide/cellbarcode/readcycle mapping remains unassembled.

🟢 Main486: NovaSeqX GEX28×10×10×90; actualATAC:CRISPRpool5:1 (ratio basis unassigned), run151×8×24×151. SI19 recommends100×8×24×100 and at least75bpPE for protospacer; these are retained as sourceactual/recommended run descriptions, not silently normalized. The schematic labels a16ntBC in thei5position; this is not assumed to be a standalone sampleindex or equal to the24cyclephysicalindexread. ATAC/CRISPR andGEX are actual separate readout pools.

## Deferred

LinkedIDTannealing and historical10xcell-line/Multiome manuals are directly linked but ordinary unread, beyond the selected two-original pass. Full bead arms/cellbarcodes/ligatedsubstrates/templateswitch products/PCRproducts and chemistry-specificLNA thermodynamics remain unassembled. No complete reaction, corrected-source or accepted missing-information audit is claimed. Biologicalpanel/drug/timepoint comparisons remain variants of this preparation.
