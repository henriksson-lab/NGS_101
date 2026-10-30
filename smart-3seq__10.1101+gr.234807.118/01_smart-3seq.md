# SMART-3SEQ

[Published defining article](https://doi.org/10.1101/gr.234807.118).

🟢 Two originals: published main and its SupplementalFile2, printed revisionJanuary8,2019. Selected standard nondegradedRNA / pre-SPRI pooling option2. FFPE/LCM variants and their distinctlysis/cyclechoices are not substituted. Main referenceRNA runs used this branch. Upstream cell collection/extraction is outside this selected purifiedRNA preparation.

🟢 Start1µltotalRNA. Add4µlfragmentationmix:2µl5×SMARTScribebuffer,1µl10mMeachdNTP,0.5µl20µM1S,0.5µl80mMMgCl2. RNA fragmentation80°C5min **precedes RT**. At42°C add5µlTS-RTmix:2µl5MTMG,1µl20mMDTT,0.5µl20×RNaseinhibitor,0.5µl20µM2S,1µl100U/µlSMARTScribe;42°C30min/70°C10min. Whole fragment positions, RNA/cDNA occupancy and template-switch products remain unassembled.

🟢 Actualwholeoligos verified originalPDFp10: 1S `/5Biosg/GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC`+T30+V. 2S `/5Biosg/CTACACGACGCTCTTCCGATCT`+N5+`rGrGrG`. N5 is the source explicitly assigned randomUMI introduced beforePCR; terminalthreeG are not UMI. 5′biotins retained. SingleindexP5=`AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T`; P7=`CAAGCAGAAGACGGCATACGAGAT`+[6nti7]+`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T`. Actual firstprintedindex `CGTGAT` retained as representative; study-specificindex allocation unassigned. PS bonds are not extra bases.

🟢 Add2.5µl2µMeachPCRprimerpair and12.5µl2×KAPAHiFiHotStartReadyMix,25µlfinal.98°C45s,cycles98°C15s/60°C30s/72°C10s,final72°C60s. Cyclecount depends on totalRNAinput and pooledsamplecount (sourceformula), not hardcoded asLCM25cycles. Model both same-source prepared template/PCR-stock components without freshduplicateenzyme/primeradditions.

🟢 Option2 physically combines indexedPCRproducts beforeSPRI. Add15×Nµlbeads/5minRT, discard actualsupernatant; resuspendin96µlRe-SPRImix(60water+36beadsorbindingbuffer),5minRT, discardactualsupernatant,2×1ml80%ethanol,airdry2min/elute10µl. Other indexedlibrary contributions are explicit source preparations; no assertion of equimolar firstpool. Main reports subsequentTapeStation/qPCR and equimolarsequencingpooling; exact normalization relationship to already-prepooledbatches remains unassembled. Actualselectedloadingaliquot andstockremainder are retained. Addactual1%PhiXDNA controlstock; NextSeq500R1=86/I7=6 forreferenceRNA (LCMR1=76 andMiSeq169 comparator notmixed). Sourcefirst5readbasesUMI/next3G, completewholeproduct/primerlanding unassembled.

🟡 Deferred: exact studysampleindex/cycleallocations, complete heterogeneouspolyA/fragmentedRNA/TSjunctions; dualindexalternative andLCM/FFPE recipes outside selectedgraph. Two originals closed with reviewer slot reserved; no missing-information claim.
