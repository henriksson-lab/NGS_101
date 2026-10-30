# Bind-n-Seq: high-throughput analysis of in vitro protein–DNA interactions using massively parallel sequencing

DefiningDOI: 10.1093/nar/gkp802. Capped source-local pass;no finished public schematic or acceptedabsence.

Sources actually inspected/accesschecked:

- https://pmc-oa-opendata.s3.amazonaws.com/PMC2794170.1/PMC2794170.1.xml — `nar_gkp802_PMC2794170.xml.txt`: Published main Methods25–59,Results63–77 and originalFigure1PDFp3 rendered;proteinpreparation scope,93mer design,resin/gel/QC/GAI-II operations
- https://pmc-oa-opendata.s3.amazonaws.com/PMC2794170.1/gkp802_1.pdf — `nar_gkp802_supp_gkp802_1.pdf.txt`: Published SI TablesS1A/B barcode inventories/condition table scope;other figures not systematically read
- Unavailable: https://pmc-oa-opendata.s3.amazonaws.com/PMC2794170.1/gkp802_index.html — supplementary file (PMC Cloud); no article content (442 characters of text)
- Unavailable: https://pmc.ncbi.nlm.nih.gov/articles/instance/2794170/bin/gkp802_index.html — PMC supplementary file; host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after reCAPTCHA challenge page (5.7 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (5.7 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested)

🟢 Published original: PMC2794170 Methods and Figure1 distinguish DNA target synthesis, protein binding, physical enrichment, normalization and sequencing. The synthetic93mer contains a2ntAA leader,3nt condition barcode,N21 binding region and Illumina primer flanks. Whole primer/flank bodies remain unassembled;no modern canonical library is substituted. N21 is the assayed binding substrate,**not a UMI**. The barcode identifies an experimental condition,not individual molecules. Figure1 originalPDFp3 was rendered and inspected;it is a design diagram,not a base-by-base primer table.

🟢 Source templates become double-stranded by primer extension (25µL:0.88µM template,88µM reverse primer,TaqProComplete;95°C2min/63°C1min/72°C4min). Protein and target then equilibrate2h23°C in the stated herringDNA/Zn/DTT/BSA/KCl environment. New DNA copies and original template/primer/reagent classes are retained;protein binding is noncovalent and no affinity or perfect occupancy is predicted.

🟢 Baseline amylose route:50µL prewashed packed amylose resin captures source MBP-fusion protein/DNA30min23°C;three10min washes →50µL10mM maltose10min →resin pelleting/supernatant transfer →repeat pelleting/transfer. Maltose releases MBP-associated complexes from amylose;the eluate is not silently assigned naked DNA or perfect removal of resin/protein.

🟢 Separate gel route:5%TBEpolyacrylamide120V30min with independentFAMcontrol complex marker →excise protein–DNA band →overnight23°C400µL ammoniumacetate/Mgacetate/0.1%SDS elution →transfer liquid →ethanol precipitation−80°C20min →pellet/wash/dry/resuspend50µL. The FAM marker is not a label on every library molecule. Unbound gel fractions are outside this baseline enrichment;no matched unbound Spec-seq library is invented.

🟢 Main calls its **real-time quantitative DNA PCR** “RT-PCR”;this is not RNA reverse transcription. A1µL diagnostic aliquot receives Pr4-f/Pr3-r/SYBR and the source40-cycle/melt/standard-curve program. The diagnostic tube remains aQC branch. The production remainder feeds a separately explicit normalization/amplification interface,whose exact production cycles/workup remain unassembled. Approximate equal-concentration experiment outputs are pooled for IlluminaGenomeAnalyzerGAI(run1) orGAII(run2),version1sequencing reagents. Exact cycle counts and primer landing remain unreviewed.

🟢 The optional extra-selection run uses a distinct source PCR/enrichment round;it is not made compulsory in baseline graph. Likewise longwash andsalt/protein conditions are source variants,not additional protocol entries. PublishedSI provides3ntbarcode families;local source-exampleAAC is not assigned a new sample or universal condition.

🟡 The local `SourceTarget.encoded_region()` describes only known internal windows,not a full93mer strand. Its barcode role is structured and exact primer chemistry remains unknown/unreviewed in construction. No source absence accepted. Family/name similarity does not establish recipe inheritance from RNA Bind-n-Seq,Spec-seq or modernSELEX.

Deferred actual source targets/questions:

- Full93mer/Pr4-f/Pr3-r primer sequences and exact geometry;Figure1 rendered original inspected as design diagram,not oligo table
- Exact source production amplification/normalization workup and run-specific readcycles/primer landing;analytical real-timePCR not silently substituted
- Source optional longwash,extra selection round/gelconditions,protein manufacturing predecessors and other supplementary figure scopes deferred

Mandatory standard source retrieval/primer scraper ran;hits are pointers,not evidence,andonly stated original scopes support the model. No two-pass absenceaudit performed.
