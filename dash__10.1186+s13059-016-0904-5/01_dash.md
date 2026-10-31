# DASH

Gu et al.(2016), *Depletion of Abundant Sequences by Hybridization (DASH): using Cas9 to remove unwanted high-abundance species in sequencing libraries and molecular counting applications*. DOI https://doi.org/10.1186/s13059-016-0904-5 .

## Originals and selected scope

🟢 Published non-manuscript PMC4778327 (identity/provenance59), Fig1 and Results on mitochondrial-rRNA depletion, Methods191–211. Published Additionalfile1/MOESM1 PDF guide-target table inspected, including the complete20base mt-rRNA-1 target `ATTTTCAGTGTATTGCTTTG`. Mandatory full-bundle scraper run; T7 template oligos checked in Methods197, not taken from its scaffold similarity hits. Other data figures/supplements remain unreviewed. Superseded preprint031708 was automatically fetched but not read or used. One capped pass, not a complete source/absence review.

DASH here is a DNA-library processing module: it cuts adapter-bearing DNA derived from RNA-seq, rather than directly degrading cellular rRNA. Source HeLa/CSF libraries use NuGEN Ovationv2 amplified cDNA and a Nextera partial-adapter frontend. Those exact historical polymers/conditions are unassembled. The KRAS/ddPCR application remains a separate branch in the same defining method and is not substituted for the RNA-seq sequencing route.

## Physical molecular contract

🟢 Cas9/sgRNA RNP creates a selected double-strand break three nucleotides before a compatible NGG PAM (Fig1). Each broken DNA piece retains its own physical sequence and the original outer adapter site on that side; the pieces do not disappear or become RNA. Amplification between the declared adapter pair requires both sites on the same molecule. The local exact-match selected-event helper checks an explicitly supplied20base recognition sequence and NGG context; it does not predict efficiency, mismatch tolerance or off-target activity, and is not a generalized activity simulator.

Source notes expected depletion of a homologous mitochondrial-rRNA pseudogene, so no global no-off-target guarantee is inferred. Digestion mixtures retain cut/uncut DNA, RNP and other material until actual workup. ProteinaseK does not physically remove protein fragments by itself. Phenol/chloroform extractions and ethanol precipitation are separate physical boundaries. PCR leaves original/cut templates and other nonamplified material in the mixture until actual selection.

## Source-owned route through readout

🟢 Prepare RNP at37°C10min in source buffer, expose tagmented DNA to it37°C2h. rRNA-depletion samples then receive ProteinaseK37°C15min, three phenol/chloroform/isoamyl extractions and one chloroform extraction; ethanol precipitation recovers DNA. This differs from the heat95°C15min/Zymo cleanup of the KRAS/ddPCR branch.

RNA-seq DNA then undergoes10–12cycle KAPA amplification with dual-index primers, BluePippin360–540bp selection, Zymo DNA cleanup and another KAPA qPCR-guided amplification. Source pooling/quantification precedes HiSeq4000PE135. The graph carries these known operations through sequencing metadata; exact historical adapter/primer/strand/run geometry remains unassembled and uncovered source items prevent finished publication. No new library-design default is supplied by importing another DASH application.

## Lineage and deferrals

The article explicitly cites an optimized sgRNA scaffold and prior assembly method (refs47/48) and source cDNA/clinical preparation predecessors (refs1/25). Their exact versions are unreviewed; neither matching scaffold fragments nor a Nextera name establishes whole-recipe inheritance. Smart-seq-total explicitly applies DASH later, but its guide sets/reagent schedules are not used here. Historical NuGEN/Nextera/index primers, complete guide-template transcription geometry and acquired primer landings remain open. No absent-information claim.

## Bounded publication review, 30 October 2026

Published original PMC4778327 Methods191–211 was re-read. Its manifest identifies the actual published original, not an author manuscript. Mandatory bundle scraping was repeated, and the exact T7 assembly/target oligos and final pooling paragraph were read directly. The source-reported pooled DNA-library output is now in the literal source inventory. This fixes one coverage omission without claiming that its complete molecules were reconstructed.

The normal publication gate still retains `unreviewed_full_input_geometry`: exact NuGEN/Nextera input polymers, complete source guide-transcription implementation and historical index/primer/read geometry remain unassembled. Some required source facts are present in the available main: Methods197 prints complete guide-template oligos, PCR primers, pooled transcription templates and T7 transcription/cleanup. The remaining work is therefore not a claim that the guide sequence is absent or unavailable. The cited optimized-scaffold and assembly originals, plus input-preparation sources, were not followed in this bounded pass. No superseded preprint or other DASH application was substituted.

Additional model work should preserve actual independently prepared DASH/no-DASH samples and their unknown pooling composition, fresh first/second PCR supplies, analytical ddPCR contributions/remainders and source-positive producer scope. The main's pooling statement does not identify every control's physical pool membership or quantitative ratio. The present generic other-library input does not establish those details. Exact dose/index/whole-strand completion and guide-producer assembly remain deferred; no coverage exclusion or source-absence claim was made.
