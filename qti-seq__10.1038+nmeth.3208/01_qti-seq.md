# QTI-seq

[Defining publication](https://doi.org/10.1038/nmeth.3208).

🟢 Source scope: associated author manuscriptNIHMS642756/PMC4344187 Methods110–142 with originalPDFp10 wholemodifiedoligos visually verified; publishedNatureSI MOESM592 Figure1/2captions; publishedNaturemain is a subscriptionpreview, not a finalMethodsbody. Three originals closed. The manuscript recipe is source-qualified rather than declared final-equivalent.

## Selected cultured-cell initiation footprint branch

🟢 Two10cmcell dishes into400µlice-coldpolysomebuffer (10mMHEPES7.4/100mMKCl/5mMMgCl2) containing5µMLTM; MatrixD6×20s with40siceintervals, clarify12000×g4°C10min and retainactualsupernatant. CHX100µgml is a separate comparator choice, not added together to selectedLTMlysate. Incubate lysate15min with25µMpuromycin,16mMHEPES7.4,10mMcreatinephosphate,0.1mMspermidine,40µgmlcreatinephosphokinase,0.8mMATP. The physical LTM/PMY treatment is intended to retain initiatingribosome populations while releasing elongating populations; exactselectivity/occupancy are not assumed. It is not a computational peak filter.

🟢15–45%sucrosegradient in10mMHEPES7.4/100mMKCl/5mMMgCl2/CHX100µgml/5mMDTT/20UmlSUPERaseIn,SW41,38000rpm100min4°C; fractionate0.375mlmin,1minfractions. The source names pooledpolysomes; exactselectedfractionnumbers/relativeweights are unassembled, not silently replaced by an isolatedmonosomefraction. E.coliRNaseI750Uper100A260,4°C1h;SUPERaseIn50Uper100URNaseI stopsdigestion,TRIzolRNAextract. LiverRiboTagantiHAcapture is a separate branch, not inherited into culturedcells.

🟢 Source15µlPNKdephosphorylation:20UT4PNK/10USUPERaseIn/1×PNKbuffer37°C1h/65°C20min. Do not silentlyaddATP or invent5Pformation here; source3OHsubstrateclassforPAP remainsqualified.15%denaturingRNA PAGE, actualSYBRGoldvisualization, excise~28ntband, crush/overnight300mMNaOAc5.5/1mMEDTA/inhibitor,SpinX/ethanol.65°C30sdenature,8µlPAPreaction1mMATP/3UE.coliPAP/.75UµlSUPERaseIn37°C45min.

🟢 Four whole5′Pbarcoded RTorderings from originalp10 are preserved literally, including source `Ø` between DNAblocks. `Ø` is not converted into an ordinarybase, assigned a particular commercialspacer, or silently dropped. Commonpre-Ø DNAarm derivesfromreversecomplementRA5 and sourceprefixes. WholeMCA02: `pCAGATCGTCGGACTGTAGAACTCTØCAAGCAGAAGACGGCATACGATTTTTTTTTTTTTTTTTTTTVN`. LGT03/YAG04/HTC05 replaceleadingthreebases byGTG/AGG/TCG. The source says barcodes; exactbarcodewindow boundary/studyassignment is unassembled (sourcevariableprefixes are not molecularUMIs). DownstreamDNAblock iscanonicalP7minusitslast3bases +T20+VN. SourceprintedRTprimerconcentration **2.5mM** retained without µMrepair;0.5mMdNTP,65°C5min/ice5min thenSSIII200U/RNaseOUT40U/TrisKClMgDTT, RTprogramme delegatesvendor.

🟢10%denaturingPAGE~100ntextendedcDNA selection/300mMNaCl1mMEDTA.20µlCircLigaseII100U/1×buffer/2.5mMMnCl2/1Mbetaine60°C1h/80°C10min. APE1 7.5U37°C1h re-linearizes sourcecircularssDNA in20mMTrisacetate/50mMKacetate/10mMMgacetate/1mMDTT; exactØchemicalidentity/cleavagebond/endproducts remain unassembled.10%denaturingPAGE~100ntlinearDNA selection. WholeqNTI200=`CAAGCAGAAGACGGCATA`; wholeqNTI201=`AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGACG`; wholecustomsequencingprimer=`CGACAGGTTCAGAGTTCTACAGTCCGACGATC`.

🟢 Phusion0.5U/0.2mMdNTP/0.5µMeachprimer/1×HFbuffer,12cycles98°C10s/60°C20s/72°C10s after98°C30s.8%nativePAGE~120bpband, BioanalyzerQC, actualequal-amountbarcodedlibrarypool withsourceotherlibrarycontributions/selectedaliquots anduntransferredstocks.3–5pMHiSeqloading. Processing25–35ntreads/12ntPsiteoffset is computational; it does not establish a complete physicalrunmask or exact initiationselectivity.

🟡 Deferred: publishedfullMethodsbody, directlylinkedvendor/predecessorrecipes, exactgradientfractions, Øchemistry/cleavageends and fulloccupiedcDNA/circle/PCR/readgeometry. No sourceabsence or completechemistryaudit.
