# ChIP-exo: source-qualified original SOLiD route

Rhee and Pugh2011, *Comprehensive genome-wide protein-DNA interactions detected at single-nucleotide resolution*, DOI10.1016/j.cell.2011.11.013.

## Originals and scope

🟢 Associated accepted manuscript NIHMS338311/PMC3243364, originalidentity checked; ExperimentalProcedures132, Results design/Fig1 read. Final publishedCell version applicability remains unreviewed. The source explicitly points to full chemistry in SupplementaryInformation DOC; its exact PMCroute was on a challenge cooldown and no chemistryDOC was obtained. The XLS target is marked TDM/gated outside the available cloudobjects. Mandatory scraper1candidate/1mainfile, not evidence of unavailable oligos. No recursively cited sources followed after this cap.

## Source operations

🟢 Original SOLiD branch: start with immunoprecipitated crosslinkedDNA on Sepharose from prior Venters/Pugh2009 ChIP (not reviewed/imported); polish DNA, ligate P2adapter, lambda5′→3′digest onresin, elute DNA, reversecrosslinks/proteinaseK65°C, ethanolprecipitate, denature95°C, extend from P2primer, ligate P1adapter at the lambda-defined boundary, PCR, gelpurify and sequence on SOLiD. Exact adapter/primer bodies, modifications, ligation ends, enzymes/conditions beyond those named, cycles and read layout remain unreviewed. No later Illumina ChIP-exo or Nexus/iCLIP recipe supplies these.

🟡 A caller-selected resection reference retains both removedDNA and remainingDNA and creates a complementary copy without changing the native strand. The selected stop is inferred, not a factor-specific physical barrier/offset/yield prediction. Digestion fragments and protein remnants persist until the actual source elution/precipitation. Denaturation separates strands without destroying DNA. Peak-pair positioning and read-coordinate shifts are computational readouts, not physical editing.

This private incomplete graph includes source-known downstream operations through SOLiD acquisition. Defining published main, chemistryDOC, cited prior ChIP (Venters/Pugh2009 DOI10.1101/gr.084970.108) and historical SOLiD reagent/oligo/run originals remain explicit later targets; no information is claimed absent and no full-source audit completed.
