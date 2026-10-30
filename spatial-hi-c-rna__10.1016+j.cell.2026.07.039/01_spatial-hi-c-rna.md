# Spatial Hi-C-RNA — bounded source notes

Defining DOI: https://doi.org/10.1016/j.cell.2026.07.039

Associated author manuscript NIHMS2201910/PMC13509204, not final/published-body equivalence. One original; corrected-proof/discovery identity does not authorize merging another version.

🟢 Source main238–286: DNA oligos are explicitly delegated to TableS7. Same-originalJATS SD8 identifies `NIHMS2201910-supplement-8.xlsx` as TableS7(primers/linkers/barcodes); otherS1–S6 captions describe contact statistics. No tablebody read/confirmed absence or claimed whole-primer constructor. Source80/200channel Ai/Bj uses 1:1annealedbarcode/linker and2xannealbuffer(100mMNaCl/20mMTris8/2mMEDTA); printed100μm not silently normalized into a derived stock model. Hi-C transposome uses Tn5MErev/ME-B **only**; ATAC's A+B preparation is a separate branch, not inherited.

🟢 Actual tissue10umsection: thaw37°C1min;1%FA/.05UulRNaseinhibitor5min;1.25Mglycine5min;3mMDSG/.05Uulinhibitor40min;glycine5min/waterwash/imaging. NP40lysis(.2%IGEPAL/10mMTris8/10mMNaCl/proteaseinhibitor/RNaseinhibitor)15min,DPBSstop,.5%SDS65°C10min,1%Triton5min,.5xDPBS5min. Adjacent H&Estain is a different tissue section, not productionRNA returned to same vessel.

🟢 RT mixture12.5ul5xbuffer/6.2ulMaximaHMinus/20.3ul.5xDPBS/9.6ulwater/.8ulinhibitor/10ul100uMRTprimer:RT30min then42°C90min. Whole primer/biotin placement/nucleotide inventory unassembled; source does not justify adding a guessed RToligo or repairing this printed recipe. Tissue stays onslide during subsequentcontactchemistry.

🟢 Hi-C digest200ul1xNEBr3.1,1Uul eachDpnII/DdeI/HinfI,.05Uulinhibitor/.1%Triton37°Covernight. WashNEBr2.1;blunt fill1xNEBr2.1/100uMdNTPmix/.2UulDNApolI RT2h, washT4buffer. Proximityligation1xT4buffer/16UulT4ligase/.1%Triton/inhibitorRTnight; **no biotin-dATP invented for this contactfill**. Then.1%SDS55°C1h,1%Triton/inhibitor5min;2xwash(10mMTris7.4/10mMNaCl/3mMMgCl2/1%BSA/.1%Tween/inhibitor).

🟢 Transposition50ul:1.87ulwater/5ulMErev-ME-Btransposome/16.5ulDPBS/25ul2xbuffer/.5ul1%digitonin/.5ul10%Tween/inhibitor,55°C1h;40mMEDTA5min/2xDPBSwash. Printed component summation/stock composition is retained without fabricated volume closure. A thenBorthogonal PDMSchannels applybarcodes37°C30min each, removing/washing chips between. Source master preparation1.5ulannealedbarcode+4ulmix gives5.5ul, while loading instruction says4.5ulperinlet: distinguish prepared stock and selected loading aliquot; leftover quantity/sourcehandling unassembled. DAPI300nMRT10min/2xwaterwash/imaging preserves spatialcoordinates; these are not extra molecularUMIs or physicalpixel aliquots.

🟢 WholeROI100ullysis(.4mgmlPK/1mMEDTA/50mMTris8/200mMNaCl/1%SDS)58°C2h; collectlysate60°Cnight. ZymoDNApurify/elute100ulwater. Actual40ulC1beads in100ul2xB&W capture100ulDNA RT1h; magnet separates **boundcDNA beads** from actual **gDNA supernatant**. No whole-lysate duplicate and no complete capture/yield claim. Source selects cDNA through streptavidin, but exact biotin-bearing oligo/nucleotide position remains unassigned without TableS7.

🟢 Supernatant→Zymopurify20ulDNA,96°C5min/ice2min. Prepare12uMsplintP5mixture30.6ulwater/6ul100uMSL_P5_rc/8.4ul100uMS5MEA_N10/5ul10xRNA-ligasebuffer;95°C1min/cool.1°Cs to10°C/hold5min; dilute.75uM anduse10ul withDNA. Ligationmaster80ul recipe includes40ul50%PEG8000/12.5ulSCR(666mMTris8/132mMMgCl2)/10ul100mMDTT/10mMATP/1.25ul10UulPNK/6.25ul400UulT4ligase;37°C45min/65°C20min. Printed80ulrecipe sums/ATPvolume unassembled, not silently repaired. Source200ulPCRmixture150ulHighFidelity2xMasterMix/5ul10uMN70X/5ul10uMN50X/40ulwater;95°C5min/10–18cycles95°C30s56°C30s72°C60s/72°C5min; .7xSPRI/elute20ul. WholeP5splintsubstrates/primerroles/N10meaning remain unassembled; suffixN10 alone is not proof of a UMI.

🟢 BoundcDNAbeads:2x400ul1xB&W(.05%Tween/inhibitor),10mMTris8(.1%Tween/inhibitor),400ulwater. SourceTSOmixture described200ul but listedcomponents sum220ul(22uldNTP/44ul5xMaximabuffer/44ul20%Ficoll/88ulwater/5.5ul100uMTSO/11ulMaxima/5.5ulinhibitor); RT30min/42°C90min, washes. Model retains source mixture description/components without silently choosing a final concentration. Source220ulPCRmixture110ul2xKAPAHiFiHotStart/8.8ul10uMprimer1andprimer2/92.4ulwater is split into50ultubes; ambiguity whether8.8ul is each primer not silently resolved. Explicit50ul aliquots and remainder are preserved.95°C3min,5cycles98°C20s65°C45s72°C3min; **product separated from beads** before1xSYBR/qPCR. qPCR95°C3min;up to15cycles98°C20s65°C20s72°C3min, stopatplateau;72°C5min;.6xSPRI20ul.

🟢 Actual2ngcDNAstock aliquot diluted5ul,10ulTagmentDNA buffer/5ulAmpliconTagment mix55°C5min;5ulNTbufferRT5min. Add25ulPCRmastermix(15ulmix/1ulN50X/1ulN70X/8ulwater); source12cycles95°C10s55°C30s72°C30s with initial95°C30s/final72°C5min. .7xSPRI finalRNA library. Both libraryclasses have separate D5000QC/loading allocations, NovaSeqXPlusPE150. Computational sharedspatial assignment does not physically merge gDNA andcDNA branches.

🟡 Private partial source graph retains source-known onephysicalmastermix components/capturefractions/aliquots, explicit unavailable P5splintgeometry and exact source volume descriptions. Detailed catalytic/end/primer products, tableS7 oligos, read landing, actual sample/branch tube pooling and historicalvendor recipes remain unassembled. No complete chemistry/absence or published-error claim.

Source investigation closed at selected1/1/3originals beforeimplementation. Mandatorybundle scraper run; sourcefacts followed toselectedoriginals, inaccessibleoligos notcalledabsent.

Deferred:

- Actual TableS7 SD8 oligo original remains gated; SI1PDF unknown chemistry relevance/body unread.
- SourcecitedinsituHiC andDiagenodevendorloading originals ordinaryunread; no inherited chemistry.
- Detailedwholeprimer/biotin/ligation/end geometry and per-tube pool/loading schedule unassembled, not absent.


🟡 One-stock inventory refinement exposes16known-but-unassembled oligo components: RT/TSO/threePCRprimerpairs/P5splintpair/loadedB-onlyTn5pair/channelbarcode-linkerpairs. They remain components of actual same stocks, not fresh free additions or channel pooling. BSA is a molecular component of the actual1%washbuffer; carryover unassigned. The default external-input source guard does not audit these derived oligos; no whole bodies guessed.
