# IPOD-HR

[Published defining article](https://doi.org/10.1371/journal.pbio.3001306), *Dynamic landscape of protein occupancy across the Escherichia coli chromosome* (2021).

🟢 Three originals inspected: publishedPMC8282354main258–310/326–328; published2022correctionPMC9045836; publishedS5TableoriginalXLSXheaders/firstWTinput,IPOD,RNAPread/kitrows (ZIP/XMLverified). Mandatory whole-bundle scraper run; hits in promoter/reportersequences were not adopted as library oligos. Appended publication-history/referee material is not chemistry authority. No vendor-kit or later detailed-protocol inheritance.

## Physical DNA fractions

🟢 Source cells at targetOD6000.2: remove200µlculture into1mlDNA/RNAShield for a separateRNAcomparator. Remaining culture receives150µgmlrifampin10min, then1%formaldehyde/10mMNaPO4crosslink5minRT and0.333Mglycine5minRT, two10mlcoldPBSwashes, dry/snapfreeze. This pretreatment is distinct from later computationalRNApolymerase subtraction.

🟢 Perpellet600µllysis10mMTrispH8/50mMNaCl/proteaseinhibitors/52.5kUmlReadyLyse,30°C15min, thenBranson25%power3×10s/10sicepauses. CalibratedDNaseI digestion:60µgRNaseA,6µlDNaseI,5.4µl100mMMnCl2,4.5µl100mMCaCl2,30minonice for the study; stop50µl500mMEDTA, usually50–200bpDNA. Stationary-phase10×dilution/one-tenth use is a separate branch, not silently applied to exponential samples.

🟢 Originalmain274 reports pre-interface clarification **169,000×g**. Publishedcorrection[10.1371/journal.pbio.3001557](https://doi.org/10.1371/journal.pbio.3001557) explicitly replaces it with **15,900×g,10min,4°C**. The model uses the published correction, while retaining both original/corrected locations and values. This is a documented publication correction, not an inferred fix. Subsequentphase-separation21,130×g2min is a different operation and is not replaced by15,900×g.

🟢 After clarification, remove50µlinput and dilute1:9 inelutionbuffer50mMTrispH8/10mMEDTA/1%SDS. Actualremaining lysate plus1volume100mMTrisbase+2volumes25:24:1phenol/chloroform/isoamylalcohol, vortex,RT10min,21,130×g2min. Retain the **white aqueous–organic interface disc**. Aqueous/organic liquid is removed; the interface is washed with350µlTE+350µl100mMTrisbase+700µl24:1chloroform/isoamylalcohol, then700µlTE+700µl24:1chloroform/isoamylalcohol. Eachvortex/21,130×g2min/liquidremoval retains the interface, then resuspend in500µlelutionbuffer. This is not recovery of the aqueousDNA fraction as inFAIRE.

🟢 A separateRNAPChIP preparation shares the sourcepreparation through clarification. Whenpossible, lysates from pellets of the same culture were pooled/mixed before taking a singleinput. ExactIPOD/RNAPbranch volumes and provenance for a particular sample pair are unassembled. RNAPbranch mixeslysate1:1 with2×IPbuffer,10µlantiE.coliRNAPNeoCloneWP023overnight4°C, then50µlProteinGdynabeads2h4°C. SourcewashA/B/C/TE and500µlelution65°C30min follow. It is a real separate antibodycapture preparation; RNApol is not chemically depleted out of the IPODinterface preparation by the computationalcorrection.

## Common recovery and sequencing

🟢 Input,interface andRNAPfractions: reversecrosslinks65°C6–16h,100µgRNaseA37°C2h,200µgProteinaseK50°C2h, phenol/chloroform/ethanol recovery with actualGlycoBluecoprecipitant/NaClsalt and95%ethanolwash. FluorescentDNAquantitation andoptional2%gel use analyticalaliquots, not wholeproductionstock. StandardhistoricalNEBNextUltraDNAE7370 (single/dualNEBindexes) replaces pre-PCRAMPurecleanups byZymoClean&Concentrator5/OligoClean&Concentrator toretainshortDNA; finalAMPure/Axygencleanup mayrepeat toreduceadapterdimers. These internal kitreaction/adapterends/primers remainunassembled. TruseqNano/Truseq variants are not merged into this selecteddefaultkitroute.

🟢 S5Tablelists paired layouts including76×76,39×38,44×44 and102×102 and differingkits. It does not byitself assign one of those runmasks to every graphbranch. SourceHiSeq/NextSeqmodels andpairedreadmode retained; no wholeprimerlanding or exactper-sampleindex layoutclaimed.

🟢 Main326–328 derivesIPOD-HRtracks from separatelysequencedinput/IPOD/RNAPsignals: inputabundancenormalization and a modifiedlinearRNAPsignal subtraction. That is computation on DNAsequencing results, not physicalDNAdepletion, RNApolymerase digestion or a new capturedDNAclass. No full reaction/source/absenceaudit; linkedvendor andpredecessor originals remainordinaryunread.
