# sci-L3-WGS

High-Throughput Single-Cell Sequencing with Linear Amplification. Molecular Cell (2019), [doi:10.1016/j.molcel.2019.08.002](https://doi.org/10.1016/j.molcel.2019.08.002).

🟢 Source: associated author manuscript PMC6874760, Methods276–298 and Fig1B/Fig2captions764/770. ONE original counted before opening; final published version, SI and TableS5 remain unreviewed. Mandatory scraper reviewed13hits from main, preserving RNA notation and phosphate ends. The selected cultured-cell DNA-only WGS route is separate from germ-cell-specific preparation, targeted PCR and RNA/DNA coassay. Positive associated-source preparation does not certify final-version equivalence.

🟢 Three combinatorial barcode rounds surround T7 linear amplification. Nucleosome depletion adapts SCI-seq xSDS (Vitak2017,10.1038/nmeth.4154); main278 explicitly modifies the buffer for LIANTI compatibility. Main284/288 describes LIANTI(10.1126/science.aak9787) as the looped-transposon/RT design influence. These support scoped component adaptation, not inherited complete recipes.

🟢 Cultured cells:1.5%formaldehydeRT10min, glycinequenchice5min, lysisbufferwash(60mMTrisAc8.3/2mMEDTA/15mMDTT).0.1%IGEPALice20min gives nuclei; pellet/washNEBuffer2.1. Nucleosome depletion0.3%SDS42°C30min/500rpm; add180µl10%TritonX42°C30min; wash2×lysisbuffer/resuspend20000nuclei/µl. Exact unrecovered fractions/crosslink products are not classified perfectly.

🟢 Both Tn5 oligo5′ ends are phosphorylated:

```
5'/Phos/CTGTCTCTTATACACATCT
5'/Phos/GTCTTG [bc1] AGATGTGTATAAGAGACAG
```

Main's bc1 text is garbled as `xrefXX`; no concrete code or width is fabricated. TableS5 is genuinely gated. Source duplex annealing95°C5min then−0.1°C/cycle/9s/cycle700cycles→25°C; diluted1.5µM. Equimolar assembly with LucigenTn5 atRT30min produces0.2µM transposome stocks; a0.1µM test dilution is a source alternative, not substituted universally.

🟢24first-round wells contain1.5µlnuclei/6.5µlwater/0.7µl50mMMgCl2/1.2µltransposome,55°C20min. Add5µl40mMEDTA/1mMspermidinestop, physicallypool, add1mllysisbuffer/pellet. Resuspend312µl dNTP/MgTrisbuffer/water/lysisbuffer and redistribute4.7µl into64second-round wells. Each receives0.8µl1.5µMhairpin plus1.18µlT4ligasemix/PEG4000,20°C30min. Source hairpin order:

```
CAAGAC [reverse complement bc2,7nt] CAGGAGCGAGCTGCATCCC AATTTAATACGACTCACTATA GGGATGCAGCTCGCTCCTG [bc2,7nt]
```

The T7 promoter is in the loop; both-sided, one-sided and failed ligations remain distinct unassembled classes. Stop4µl, pool~630µl, stainDAPI5µg/ml and sort100–300nuclei per new well into3µllysisbuffer. A chosen well and unresolved other wells/uncollected nuclei are explicit; this is not one nucleus per well.

🟢 Lysis75°C45min, Qiagenprotease2mg/ml55°C8h,75°C30min; optional−80°Cstorage. BstWarmStart2.0 gap extension68°C5min precedes20µlHiScribeT7IVT37°C10–16h. EDTA2.2µl0.5M stops transcription; RCC5/elute18µl0.1×TE. Add0.6µlRNAprimer, dNTP/RNaseinhibitor before70°C1min/90°C20s/ice. LiteralRNAordering is `rArGrArUrGrUrGrUrArUrArArGrArGrArCrArG`, retained asRNA rather than substitutedDNA. SSIV runs55°C15min/60°C10min/65°C12min/70°C8min/75°C5min/80°C10min. Foldback RNA priming is possible with both-ended ligation; added RNA primer supports one-ended templates. The graph tracks the added RNA primer as a co-retained same-tube component after annealing, consumed as a structural primer at RT without a second stock addition. Possible foldback substrate classes remain unassembled; no intramolecular primer geometry is asserted.

🟢 RNaseH/RNaseA37°C30min cleaves transcripts and freeRNAprimers before Q5 second-strand synthesis. SSSprimer `NNNN [UMI] ZZZZZZ [bc3] GGGATGCAGCTCGCTCCTG` adds4ntUMI to individual IVTtranscripts, not original genomic DNA, plus6ntthirdcellbarcode. DCC5 and historicalNEBNextUltraII finishing/minimal3PCRcycles make WGS library. Symmetric finishing permits eitherP5/P7orientation; targeted PCR and RNA/DNA coassay finishing are not inherited.

Two text-matching exceptions cover the numbered hairpin components and adjacent annotated UMI/bc3 windows. Raw source patterns/full fixed-base fixtures remain locally checked; they are not concrete study codes. Needed: associatedSI and actual TableS5 original `https://pmc.ncbi.nlm.nih.gov/articles/instance/6874760/bin/NIHMS1537467-supplement-5.xls`, final published version and historical finishing kit. Exact end/ligation/RT/read geometry remains unassembled. No absent-disclosure or complete validation claim.
