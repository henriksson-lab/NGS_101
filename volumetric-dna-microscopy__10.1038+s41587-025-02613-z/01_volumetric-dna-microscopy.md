# Volumetric DNA microscopy — bounded source notes

Defining DOI: https://doi.org/10.1038/s41587-025-02613-z

PUBLISHED CC-BY-NC-ND PMC12339006.2 (not associatedAM); publishedSI1 andactualSupplementaryTable3MOESM5XLSX. Three originals. Superseded2023preprint and later companionprotocol are not recipe authorities. Modality is RNA through DNA/RNA intermediates.

🟢 This method sequences RNA-derived cDNA and proximity-event records. It differs from original DNA microscopy by cDNA anchoring, pre-circularized identifier templates and localized RCA before IVT/unanchoreddiffusion. Source main40–76 andSI1S1–S4 distinguish molecule UMIs, UEIs marking joiningevents, and6ntsamplebarcodes(XXXXXX inS3caption). RandomN8RTpriming andN6adapterbinding windows are not automatically UMIs.

🟢 Selected24hpfembryo source190–240: pronase1mgml28°C5min/PFA4%4°Cnight/100%methanol15–30minRT/−80°C≥2h, rehydrate75/50/25%methanol5mineach/4xPBST. ThermolabilePK5e−5UulRT12min/55°C15min/4xPBST. PreRT4°C1h then65°C10min (20%formamide/.5UulSuperase/4.4mMDTT/.5ugulalbumin/PBS),cool4°C/waterrinse. RT actual400uMdNTP/**32uMaminoallyldUTP**/.5ugulalbumin/1UulSuperase/10UulSSIII/1uM21.068C-8N/1xFSbuffer;4°C1h/60°C3min/37°Cnight/50°C1h. Source optimizationSummary refers to32nM in **initialcultured-cell optimization**; retain both contextual values and choose detailed32uM for end-to-endembryo. No asserted publicationerror or silent1000x repair. ExoI1.43Uul4°C1h/37°C1h removesRTprimer/displacedcDNA; washes retain preparation, not perfectyield.

🟢 Table3 transposome:22tn5.003A has firsttwo phosphorothioatebonds `A*C*...`;22tn5.MOS-3p is `CTGTCTCTTATACACATCT/3Phos/`, **3P not substituted with familiar5P**. Oligos100uM anneal1:195°C2min/cool1°Cmin70min, then25uMduplex/1ugultagmentase23°C30min,50%glycerolstock−20°C. Stock3:10dilution then1:50reaction addition in5mMMgCl2/10mMTris7.5/10%DMF/9%PEG8000/850uMATP,4°C1h/55°C1h. BS(PEG)5 5mMPBSRT1h crosslinks aminoallyl-bearing preparation;1MTrisquench30min;4xSSC4°C10min/70°C15min;10%formamide2xSSC4°C10min/50°C10min/washes. No identical anchoring efficiency inferred for every molecule.

🟢 Actual3′adapters22tn5.005/.006 have5P, randomN6binding window, source two terminalphosphorothioatebonds, and3P. Source500nMeach/SplintR1.25Uul4°C1h/23°Cnight; washes, **QuickCIP.5Uul removes3P afterligation** at4°C1h/37°C1h, thenwashes. Known local5P/N6RNA-splinted adapter-end marks are typed, preserving3P/PS; complete occupied RNA-splint/DNA-junction geometry remains unassembled for heterogeneous binding; no guessedtargetRNAsequence.

🟢 Main namescirc6G1/circ7G1 andTable3 explicitly maps linear5P21.006G1/21.007G1 withsplint6F5/splint7F5 to those preparations. SameSI1S1/S2 labels24.006G1/24.007G1 andsplint6/splint7; names remain as printed and wholeversionequivalence is not assumed. Local closure check uses **actualTable3** endsequences andshortsplints to verify terminalcomplementarity; heterogeneousinternaldegenerateUMIs remain sourceclasses. T4ligase/splints→ZymoOligoConcentrator/TBEureaQC. CircDNA100nMeach/2xSSC10%formamide.1%Tween40°Cnight/washes.

🟢 Select enhancedpreRCA25ngulT4g32/1mMdNTP/.04UulT4DNApol/NEBuffer2.1; normal100uMdNTPalternative notmerged. RCA25ngulg32/.5ugulalbumin/**250uMd(AUGC)TP**/.2Uulphi29;4°C1h/30°Cnight. Source dU substitutesdT so USER can later digestnanoballs. Optional20uMfluorescein12dUTP imaging experiment is a separate analytical condition, not automatically added to unlabelled end-to-endproduction.

🟢 UEIanneal100nM21.004G1/2-BCmix/100nM21.073pt/100nM21.074B in2xSSC5%formamide.1%Tween4°C1h/50°C2h/washes;extendligate1mMdNTP/.15UulT4DNApol/20UulT4ligase+ATP4°C1h/RT40min/washes. Actual21.074B has5P/3P plusPS;21.073pt hasPS;21.075T7promoteroligo has3P. IVT30ul includes100nM21.075/7.5mMrNTP/100ngulg32/.05UulUSER/10%T7enzymeMix/buffer/6.4ugulTHIOCUREETTMP/73.6ugul4armPEG20Kvinylsulfone. FormgelRT2h then37°C20h. **USER releases products from dU-containing nanoballs**, T7produces aRNA, and unanchoredaRNA can generate longer-rangeIVT-UEIs. RCA-UEIs andIVT-UEIs are chemicalproductclasses; the eventual spatialreconstruction is computed. Do not simulate exactdiffusion/interaction fidelity as molecular proof.

🟢 Add12uldenaturation457.5mMKOH/100mMEDTA/42.5mMDTT4°C2h;12ulstop600mMTris7.5/.4NHCl→54ul; add30ulPKmix(.28%Tween/.09UulPK/8.6mMTris7.5)→84ul/50°C1h. 1.2xRNAClean/elute;DNaseI.1Uul/Superase.8Uul37°C30min/1.2xRNAClean. ExternalRTprimers21.077/21.085500nMeach plus500uMdNTP first65°C5min/ice; thenSSIII10Uul/FSbuffer/5mMDTT/Superase1Uul50°C1h/70°C15min. ExoI3.3Uul37°C30min/80°C20min.

🟢 Split actual postRTstock allocations into **fourcDNA PCRs** (template1:80) and **twoUEI PCRs** (1:40); remainder retained/unassigned, not six copies of entirestock. cDNA primers21.046G1-BC/21.081b300nMeach/HiFibuffer/200uMdNTP/2mMMgSO4/.02UulPlatinumTaqHiFi;95°C2min/5cycles95°C30s56°C30s68°C2min/20cycles95°C30s68°C2min/68°C5min. UEI primers21.077-G1/21.076BB300nMeach plus3.3uMeach **3Pblocked**4E4.interf1/4E.interf2 and5%DMSO;95°C2min/1cycle95°C30s66°C30s68°C2min/18cycles95°C30s68°C2min/68°C5min. Interferenceoligos are not extendablePCRprimers; no dephosphorylation inserted aftertheiraddition. .75xAMPure,quantify,sourceNextSeq500150cyclekitwith112R1/44R2 andcustomspikeSBS3B. Preserve reportedkit/readnumbers without guessingmask/indexcycles or silently shorteningreads. Actualreplicatepool/loading schedule unassembled.

🟡 Private graph retains sourcewholemodifiedoligoinventory, sourceknown endmods/features, computedlocalcircle/splint relation,actualsourceanchoring→RCA→UEI→USER-IVT transition andpairedcDNA/UEIlibraryclasses. Wholeheterogeneousbindings/ligation/amplicons, actualbarcodevaluechoices, completephysicaldiffusion andfinalprimerlanding remainunassembled. No primers/kit chemistry imported from originalDNA microscopy or companionprotocol; laterpolicy/absence/sourcegate notcomplete.

Source investigation closed at selected1/1/3originals beforeimplementation. Mandatorybundle scraper run; sourcefacts followed toselectedoriginals, inaccessibleoligos notcalledabsent.

Deferred:

- PublishedSupplementaryTable2reagentlist andreportingsummary notread; mainrecipesonly, no uninspectedvendorformulation imported.
- WholeheterogeneousRNA/adapter/context andactual UMI/UEI/readlanding products unassembled; publication sourcealias24-versus21 bodies notassertedequivalent.
- OriginalDNA microscopy/splint circularization predecessor/companionprotocols ordinaryunread, no recipeinheritance.


🟡 Final narrow peer refinements retain known local5P/N6RNA-splinted adapter-only interfaces without asserting whole occupied products, and expose the actual source1mgmlPronase catalytic input.37focusedchecks include23independent originalTable3wholemodifiedordering fixtures. Three default-guard unrepresented inputs contain explicit nonidentifier N8RT/N6adapter priming windows; no false UMI feature or exemption masks this limitation.
