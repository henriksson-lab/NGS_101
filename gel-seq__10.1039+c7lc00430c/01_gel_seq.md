# Gel-seq: whole-genome and transcriptome sequencing by simultaneous low-input DNA and RNA library preparation using semi-permeable hydrogel barriers

Defining DOI: https://doi.org/10.1039/c7lc00430c

## Published original and scope

🟢 Hoople et al., Lab on a Chip17:2619–2630(2017), DOI10.1039/c7lc00430c. The institutional PDF begins with an eScholarship cover sheet, followed by the publisher-typeset article carrying this DOI, journal pagination and22June2017publication date. This is the published original, not the later2018JoVE protocol. Read Experimental overview and Gel-seq protocol onpp2620–2624; discussion confirms that genomic MDA was tested but not required in the chosen route. The publisher HTML returned403 to the standard fetcher; its abstract was identity-screened via web. Publisher-linked supporting information returned404. The institutionally hosted article remains available despite those failures.

Mandatory oligo scraper was rerun after the institutional PDF became available. The one relevant primer hit was verified in the article. This is a low-input100–1000cell method, not a single-cell barcode protocol or an inferred G&T-seq recipe.

## Actual separation substrate

🟢 Intact cells enter diluted SMART-Seqv4 lysis buffer. Source adds3′SMART-SeqCDSPrimerII plus a random hexamer with SMART handle: `AAGCAGTGGTATCAACGCAGAGTACNNNNNN`. This N6 is a random priming domain, not an asserted UMI. Lyse72°C3min, then SMARTScribe/template-switch reaction42°C90min and70°C10min inactivation. Genomic DNA remains in the same sample. The separation substrates are genomicDNA and RNA/cDNA hybrids: reverse transcription occurs BEFORE separation, not after isolating nativeRNA.

🟢 In the polyacrylamide device, the low-density4%T/3%C layer stops genomicDNA while permitting smaller transcript-derived material to pass; the30%T/5%C layer retains RNA/cDNA hybrids. Apply210V30min, then excise spatially distinct gel regions. The local GelFraction objects encode source-declared populations and positions only; they do not calculate pore size, perfect partition, migration velocity or losses.

## Physically independent branches

🟢 Each gel region enters a separate tube. Crush/soak with40µl water for genomic material or80µl for RNA/cDNA,37°C8–12h. An8µm mesh filter at2600RCF5min removes gel fragments. Only the genomic fraction receives Qiagen protease50°C15min followed by70°C15min. Protein digestion is not a physical purification; nonvolatile reaction products persist during the subsequent vacufuge concentration to5µl. The RNA/cDNA sample is concentrated separately to10µl. Evaporated water is explicit waste; sample concentration does not erase protein/RNA contaminants.

🟢 The transcript fraction undergoes20–30cycle qPCR using one source PCRPrimerIIA species with KAPASYBRFast, choosing cycles from the qPCRcurve, then AMPure cleanup under the historical SMART-Seq manual. Its unknown inward primer sites are explicitly unavailable; symmetry is declared only because the source specifies one species, not because a pool is silently treated as one primer.

🟢 Both branches subsequently use half-volume NexteraXT protocols. These remain unassembled terminal handoffs; no kit adapter topology, selected indices or run schedule are invented. The graph has no finished library or acquired read output.

## Deferrals and lineage

Publisher supporting information, exact historical SMART-Seqv4 oligos/manual and NexteraXT instructions remain follow-up originals. Full step-level kit workup and final strands/read layouts are unassembled, not confirmed absent. Device fabrication materials are source-described but not a computed polymer model. The article explicitly attributes modified Smart-Seq to Ramsköld and cites NexteraXT; these support component-level lineage suggestions only, without copying an unspecified or current kit revision.

## Completion

Bounded original-source pass only. Catalogue remains notes. Private partial graphs retain explicit uncovered source inventory and unresolved terminal handoffs; passing structural/local checks does not complete source review. No two-reader missing-information pass has been performed.

## Bounded publication follow-up (2026-10-30)

The published institutional article's Experimental section explicitly sends both recovered branches to separate half-volume Nextera XT preparations (text492–504). The current model stops at concentrated genomic DNA and amplified/purified cDNA, before those library-producing preparations. The required `unassembled_source` item therefore remains uncovered; a downstream sequencing label alone cannot complete the model. The cited originals are Illumina Nextera XT DNA Library Preparation Guide Part15031942 Rev.E (2015), reference13/text987–988, and Clontech SMART-Seqv4 Ultra Low Input RNA Kit User Manual (2016), reference14/text989–990. Their applicable detailed chemistry remains unreviewed; current guides are not substitutes.

A future bounded implementation should preserve actual SMART and KAPA supply components, independent genomic/cDNA library paths, QC and loading stocks without inventing occupied polymers or kit formulations. This pass reread published text418–438,477–504,987–990 and ran the mandatory existing-bundle scraper. One original; no acquisition, model changes or builds. Notes status remains a model-scope hold with real source dependencies, not an absence finding.
