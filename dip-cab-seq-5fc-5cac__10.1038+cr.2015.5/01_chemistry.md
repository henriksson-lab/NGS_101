# DIP-CAB-seq: 2015 5fC / 5caC channels

Lu et al., *Base-resolution maps of 5-formylcytosine and 5-carboxylcytosine reveal genome-wide DNA demethylation dynamics*, Cell Research (2015), DOI[10.1038/cr.2015.5](https://doi.org/10.1038/cr.2015.5).

Private incomplete operational graph: `tools/chem.py`. The5fC and5caC channels are branches of this protocol, not separate entries. The2021 same-acronym5hmC/deamination assay is distinct and not used as evidence here.

## Bounded source scope

🟢 Exact published supplement, `source_pdf_49222db245fceaa530e3.pdf.txt`, printed title/authors1–17; Methods20–32,62–114,164–181; FigS1A223–227. Mandatory final-bundle scraper returned14candidates in one file; actual source lines were read. Oligos there serve dot-blot, spike-in/qPCR and Sanger validation, not an asserted NGS adapter inventory. No complete irrelevance or missing-sequence audit is claimed.

The official publisher `cr20155.html` is a free-account preview, not a full chemistry article. Its direct published SI link was retrieved by the standard fetcher. The articlePDF returned HTML; EuropePMC XML returned500, PMC page is in challenge cooldown, and official PMCCloud returned no usable file. These are actual route outcomes, not proof the article is universally unavailable. The capped source investigation is closed after the preview and published SI. Full main reading remains blocked; the positive SI chemistry supports the private model.

## Actual source preparation

🟢 Start with100µg sonicated Tdgfl/fl or Tdg−/− mESC genomic DNA, average400bp, in four parallel5fC or5caC IP reactions. The SI explicitly adopts the IP method above: antibody4°C4h in phosphate/NaCl/Triton buffer, prewashed ProteinG Dynabeads4°Covernight, three IP-buffer and three higher-salt washes, PK55°C2h, phenol/chloroform extraction and ethanol precipitation. The four recoveries are combined; FigS1A divides the enriched DNA into two halves for traditionalBS and chemical-assistedBS. Neither antibody enrichment nor wash partition is assumed perfect. PK degradation products remain until actual extraction.

🟢5fC protection:10mM O-ethylhydroxylamine in100mM MESpH5,37°C4h.5caC protection:20mM NHS/2mM EDC in75mM MESpH5,37°C0.5h; buffer exchange to100mM sodium phosphatepH7.5/150mMNaCl;10mM(4-aminomethyl)benzylazide37°C1h. The installed protection classes are kept distinct from native bases. There is no click-biotin capture step in this recipe. Protected DNA receives Zymo DNA Clean&Concentration cleanup followed by MethylCode bisulfite treatment; the source says to run the bisulfite thermal program twice for5fC. TraditionalBS and CAB remain independent physical preparations.

🟡 The local selected-state relation explicitly distinguishes an unprotected target's physical BS-derivedU, its reported referenceT and complementary first-copyA. A selected protected target retains a C-class and reportsC with complementaryG. This is a conditional selected molecular class, not a conversion-yield, chemical-structure or antibody-specificity prediction. Exact adduct bond geometry, heterogeneous conversion/degradation populations and kit recovery remain unassembled.

🟢50ng treatedDNA is allocated to EpiGnome Methyl-Seq adapter tagging/amplification; unused treatedDNA remains an independent unresolved fraction. The historical EpiGnome/MethylCode kit internals, sequences, enzyme choices and full endpoint geometry are not filled from modern kits. BioanalyzerDNA1000 consumes analytical aliquots; the production remainder becomes physical sequencing library. Genotype-matched, **nonenriched** genomic input controls are independent source preparations, not IP flowthrough; their exact chemical/library recipe is left opaque. Physical libraries are pooled with other study libraries and sequenced on HiScan rapid two-lane flowcells, SE100, TruSeqSBSv3-HS.

## Reporting, links and limits

🟢 Source comparison uses converted/unconverted counts, coverage and Fisher testing; those computational criteria do not alter DNA or establish perfect true-target identity. The two biological replicates/genotypes are not instantiated as invented sample-index assignments. Full adapter/index/primer/read geometry remains unassembled.

Directly cited compound synthesis original: Lu et al.2013, *Chemical modification-assisted bisulfite sequencing (CAB-Seq) for5-carboxylcytosine detection in DNA*, JACS135:9315–9317, SIreference9/source354–356. Its full original was not chased after this capped pass. SIreference2, Song et al.Cell2013 (the existingfCAB defining paper), is cited for nucleaseP1/alkaline-phosphatase digestion in analytical nucleoside validation (source126–127), not as proof of inherited NGS library chemistry. The separate genomic-preparation prose cites “Song et al.,2011” (source21–22), whose exact identity/applicability is unreviewed; it is not silently equated to that2013 reference. No module or whole-family recipe is imported. No missing-information claim is accepted.
