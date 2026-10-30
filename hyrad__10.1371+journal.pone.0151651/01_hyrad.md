# hyRAD — bounded source notes

Defining DOI: https://doi.org/10.1371/journal.pone.0151651

🟢 One PUBLISHED CC-BY originalPMC4801390, main149–250/Table3; actualPDFTable3p7 andmethodsp8 rendered. No updatedGitHub or earlierRAD/degraded-library recipe substituted.

🟢 Probe factory: freshgDNA6ul/9ulreaction,MseI1U/SbfI-HF2U37°C3h;2xAMPure10ulTris. AnnealedRAD-P1/P2.5uMeach/T4ligase400U20ul16°C3h;1xAMPure/Pippin270bp tight.10ulQ5HotStartPCR(.2mMdNTPeach/.6uMprimers/.2Upolymerase),30cycles98°C20s60°C30s72°C40s. Actualreplicates empiricallyprovide500–1000ngpercapture;Qubit/FragmentAnalyzerQC aliquots, poolreplicates/1xAMPure. SplitRADlibrarystock: actualaliquotSE300MiSeqreference, remainingstock adapter-stripping MseI5U/SbfI10U50ul37°C3h/65°C20min/1.5xAMPure, BioNicknickbiotinlabeling/1.5xAMPure yieldsDNAbaits. BioNickwholeformulationunassembled; not RNAbaits or syntheticdesignedpanel.

🟢 Separate targetDNA: sourcebutterflyindividualDNA splitintointact andCovaris~300bp aliquots, processedparallel; sourcegrasshopperusesunsonicated branch. Selectedbutterflyintactbranchretainsotheractualsonicatedstockunresolved. DNA8ul95°C10min/ice;T4PNK10U37°C30min65°C20min/2xAMPure. TdT10U/4mM**GTP**/CoCl2.25mM/NEB4 with10ulDNA20ul37°C30min70°C10min afterrepeatdenature. ActualPDFp8confirmsGTP, not silentlydGTP; tailribose/deoxyriboseformunassembled.

🟢 Addsecondstrandmix10ul(1ul10xNEB4/.6ul25mMdNTPeach/1ulP2**15mM**/5.4ulwater/2ulKlenowexo−5Uul) to20ulTdTmix23°C3h75°C20min. ActualPDFp8confirms15mM, not silentlyuM. CarryexistingTdTreaction; do notdiscardandrestartDNA. Bluntmix.5ulNEB4/.35ulBSA10mgml/.2ulT4DNApol3Uul/3.95ulwater12°C15min;2xAMPure10ul. P1adapter1ul25uM/T4ligase400U20ul16°Cnight/1xAMPure20ul; Bst12U/.25mMdNTPeach/ThermoPol40ul37°C20min fill. Selectactual10ulof40ulstockintoPCR, retainrest: add15ulpreparedmixQ5buffer/dNTP/primerpair/Q5HotStart.5U;25cycles98°C10s60°C20s72°C25s;gelQC/.7xAMPure;equimolarsamplelibrarypool.

🟢 Capture40ul:6xSSC/50mMEDTA/1%SDS/2xDenhardt/2uMeachBO1–4/500–1000ngprobes/500–1000ngshotgunlibrary;95°C10min65°C48h. HumanCot1DNA500ng is sourcegrasshopperalternative, not automaticallyaddedtoselectedbutterflycapture.10ulM280beadswashed3xTEN/resuspend200ul;add40ulhybmixRT30min. Retainbeads, removeactualsupernatant;200ul65°C15minwashes2x1xSSC.1%SDS,1x.5xSSC,1x.1xSSC. Elute**targetlibraryoffprobes**30ul80°Cwater10min, retainedprobe/beadstockseparate. Actual15ulof30ul eluate→50ulQ5PCR(1U/.2mMdNTP/.5uMprimers),25cycles/1xAMPure;equimolarlibrarypool/MiSeqPE150. GrasshopperHiSeqPE100 branch not merged.

🟢 Table3 I=inosine (sixI inBO3/4), not randomUMI. X6inlinebarcodes and6ntILLPCR2indexes identify samples, actualXvaluesunassembled. SourceILLPCR1 keptverbatim; wholeprimerlanding/possiblegappedalignmenttoreportedP1adapterbodies notvalidated or silentlychangedtocanonicalprimer. Wholemodifiedorders:

- RAD-P1.1: `ACACTCTTTCCCTACACGACGCTCTTCCGATCTxxxxxxCCTGCA`
- RAD-P1.2: `GGxxxxxxAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
- RAD-P2.1: `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`
- RAD-P2.2: `TAAGATCGGAAGAGCGAGAACAA`
- P1.1: `ACACTCTTTCCCTACACGACGCTCTTCCGATCTxxxxxx`
- P1.2: `xxxxxxAGATCGGAAGAGC`
- P2: `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCCCC`
- ILLPCR1: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTT`
- ILLPCR2_01: `CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_02: `CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_03: `CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_04: `CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_05: `CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_06: `CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_07: `CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_08: `CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_09: `CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_10: `CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_11: `CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCAGACGTGTGC`
- ILLPCR2_12: `CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACGTGTGC`
- BO1.P5.F: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`
- BO2.P5.R: `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT`
- BO3.P7.F: `CAAGCAGAAGACGGCATACGAGATIIIIIIGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`
- BO4.P7.R: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIATCTCGTATGCCGTCTTCTGCTTG`

Source investigation closed before implementation at1published originals; mandatorybundle scraper run. Wholeproduct/linkedsource/absenceaudit unfinished.

Deferred:

- PublishedcitedRAD/degradedDNApredecessors/vendorBioNick instructions ordinaryunread; no inheritedrecipe.
- Actualupdatedprotocol https://github.com/chiasto/hyRAD is explicitlyupdated, currentversionunread/notassumed2016equivalent.
- ExacttailGTPsugar/P2concentrationinterpretation/ILLPCR1occupiedlanding andwholeheterogeneousproductsunassembled, notacceptedabsence.

The private graph keeps the source 24 whole oligo orders and computes the source-local RAD-P1 terminal3′TGCA, RAD-P2 terminal5′TA and shotgun-P1 paired blunt interface. Exact occupied substrate/end modifications and sealed molecules remain unassembled. The printed ILLPCR1 is retained without a canonical replacement or invented exact3′landing. Four narrow source-text exceptions cover Table3's literal `xxxxxx` sample-barcode notation; whole-order fixtures and the originalp7 image verify their flanks and six-position width. BO3/4 retain physical inosine descriptors; the drawing permits literalI with a placeholder flag solely because its DNA alphabet excludesI, not because these are random or unidentified bases. The source guard does not validate those two complete inosine-containing inputs or produced PCR components.

The modeled selected butterfly branch has actual intact target material and a separately retained Covaris300bp branch; the latter is not merged into the selected library. ActualRADreference SE300 is a separate source aliquot alongside finalcapture MiSeqPE150, within this one protocol. PCRreplicate, normalizedshotgun-library and captured-library pools expose contributing physical classes and untransferred stocks. Whole stock counts/amounts remain unassigned where not specified. The10µl second-strand mix,5µl blunting mix and15µl PCR mix expose components from their actual prepared additions; other reaction recipes are inventory groupings with unknown packaging, not invented premixed stocks.

The RAD PCR replicate pool is formed from distinct pre-PCR template allocations: one actual reaction and the other separately run reaction array, with empirical count and doses per source recipe. Qubit and FragmentAnalyzer aliquots leave each product class before physical pooling. A single amplified product is never relabeled or split to manufacture independent PCR replicates.
