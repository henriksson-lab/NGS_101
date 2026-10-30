# ARTR-seq — selected standard PTBP1 preparation

[Profiling of RNA-binding protein binding sites by in situ reverse transcription-based sequencing](https://doi.org/10.1038/s41592-023-02146-w), NatureMethods2024. Inspected published PMC10864177 mainStrategy/Methods186–246/ExtendedData1biotin/RTasecontext; actual mainPDFp2; SI1MOESM1 SupplementaryFigure1/Tables1–3, actualPDFp3. Two selected originals; no whole-source or missing-information audit.

## Distinct physical preparation

🟢 ARTR uses antibodies to localize a proteinA/G–reverse-transcriptase fusion to an RBP, then copies nearby RNA in situ with biotin-bearing dNTPs. It does not substitute UVcrosslink/IP, RT&Tag oligo-dT/tagmentation or a CLIP library constructor. Select standard40kHepG2/PTBP1 cells and the pAG-engineeredMMLV(25–497),30aa-linker fusion chosen in the Results. Its upstream protein-expression/purification recipe remains a prepared-stock boundary. Formaldehyde can retain indirect neighboring RNA/protein interactions; proximity labeling is not proof that every read represents a direct protein–RNA contact.

🟢 Fix1.5%PFA10minRT, quench125mMglycine5min, DPBSwash2×, permeabilize0.5%Tritonice10min, DPBSwash2×. Block1mg/mlBSA/.2U/µlRNaseOUT30min. Primaryantibody inblockingbuffer1h (manufacturerIF dilution or1:200 absentguidance), wash3×DPBS3min; fluorescentsecondary1:50030min/wash3×; pAG-RTase10nMinblockingbuffer30min/wash3×. The fusion protein is already bound when RT starts: no second fresh RT-enzyme addition is invented.

🟢 The source prepares the RT mixture in50µlbuffer (final-volume convention/individual stock volumes unassigned), supplying2µMadapter-N10primer, .05mMbiotin16-dUTP/.05mMbiotin16-dCTP/.05mMdTTP/.05mMdCTP/.1mMdATP/.1mMdGTP/1UµlRNaseOUT inDPBS+3mMMg. React37°C30min, stop20mMEDTA/10mMEGTA3minRT. Modified DNA residues contain biotinylateddU/dC; dU is not an RNA ribonucleotide or a claim about the base reported after PCR. Stainanti-biotinAlexa4881:2001h, Hoechst33342 at1µg/ml15min and image. Physical antibody/dye contributions and carryover remain explicit.

🟢 PK1mg/ml37°C2h precedes phenol/chloroformpH8 extraction/ethanol concentration. DigestRNA withRNaseH and1:20RNaseA/T1 in50µlbuffer37°C1h. Enrich biotin-cDNA on10µlC1streptavidinbeads20minRT; sourcebeads were preblocked30minRT with1µg/µl eachBSA, salmon-spermDNA andyeasttRNA. Those are real molecular blockers, not free fresh capture probes. Complete RNA digestion, capture specificity and yields are unassembled.

🟢 Suspendcapturedbeads inligationmixture50mMTris7.5/10mMMg/25%PEG8000/1mMATP/1UµlT4RNAligase1/5µM5P-N8DNAadapter,25°C16h. Releasebiotin-cDNA in95%formamide/10mMEDTA8 by95°C10min, EtOHprecipitate and dissolve10µlwater. Amplify40µlreaction with1×NEBNextUltraIIQ5MasterMix, all10µlcDNA and.5µMIlluminasequencingprimers; E7335S is the source's example, not proof of the exact chosen whole primer/index sequence. PCR98°C30s,13cycles98°C10s/60°C30s/65°C45s,65°C5min. Select180–400bpusing6%NovexTBEgel. SourceNextSeq550/2000 orNovaSeq6000 platforms; exactstandardrunlength/readprimerlanding remain unassembled.

## Oligo and identifier scopes

🟢 Standard adapter-RTprimer: `AGACGTGTGCTCTTCCGATCTNNNNNNNNNN`. Its N10 is random priming sequence, not assigned a UMI or cell-barcode role. The optional low-input/tissue RTprimer inserts a distinct8ntbarcode beforeN10; that branch is not merged into the standard selected construct or used to invent a baseline cell barcode.

🟢 3primecDNAadapter: `5′Phos-NNNNNNNNAGATCGGAAGAGCGTCGTGT-3′SpC3`. The source8ntUMI is in this adapter, supported by primarydata-processing Methods. SpC3 is a nonnucleotide3endblock, not extra S/C DNA bases. Known5Pdonor and full sourceadapter are retained; completeoccupiedcDNA/ligationproduct geometry remains unassembled.

SI1Table1 consists of analyticalqPCR primers, andTable2 of knockdownsiRNAs. The actualPDFp3 header/rows were inspected; neither is substituted for the source's library amplification primers. SI1Table3 describes samples and comparisonreadcounts. No historical NEB primer/kit recipe, barcode whitelist, complete end geometry or primer landing is supplied by brand-name similarity.

The input control omits the primaryantibody; low-input uses4%PFA and pools two barcoded biologicalreplicates afterPK. Tissue/imaging/stress/strongRNase variants have their own source scope. They are not fake splits of the selected sample or counted as separate protocols. FulltargetRNA, modifiedDNApositions, fusion/antibodyoccupiedcomplexes, all offproducts and exact library/read layouts remain unassembled, not confirmed absent. This is a private source-qualified partial preparation.

The PCR primer supply inventory exposes forward- and reverse-role primer classes with whole bodies unassembled. Their actual molecular species count, stock count and packaging remain unassigned; this does not invent an exact two-tube premix or E7335S sequence construction.
