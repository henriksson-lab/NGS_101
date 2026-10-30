# dRNA-seq — original TEX/mock primary-transcript libraries

Sharma et al. *The primary transcriptome of the major human pathogen Helicobacter pylori*. Nature464,250–255(2010). [Published article](https://doi.org/10.1038/nature08756).

## Sources and boundary

Read the published main identity/abstract and its subscription-preview notice, actual published Supplementary Methods206 oligo-role statement113–114 and central preparation237–282, companion-library paragraph287–296; actual published Supplementary Tables208 Table2 and rendered originalp6. Three originals, selected before opening; source investigation closed. No specimen-culture instructions, remaining supplement originals, vendor/predecessor recipes or superseded/later protocol versions were adopted. The mandatory oligo scraper was run (600hits in3files) and relevant source pointers reviewed; analytical hybridization/5′RACE primers identified by Methods/Table17 are not assigned as library oligos. The full article remains an actual access barrier, but its published supplement supplies the core preparation recipe. This is not a full source or absence audit.

Differential RNA degradation enriches primary transcripts relative to a matched mock library. The selected model covers the original ML−/ML+ pair as two sequencing outputs, not two protocols. Other sample conditions/replicates are branches. The paper's separate Illumina operon libraries are described below and are not merged with this 454 workflow.

## Published reaction order

🟢 Supplementary Methods237–254: remove residual genomic DNA using DNaseI, then divide equal RNA amounts between Terminator5′-phosphate-dependent exonuclease(TEX, EpicentreTER51020) and buffer-only arms. TEX is1U/1µg RNA,60min30°C; mock receives the same incubation without TEX. Source supports selective enrichment, not complete disappearance of all processed transcripts. Enzyme cleavage products remain in the reaction until actual workup.

🟢 After TEX/mock, extract with phenol/chloroform/isoamylalcohol25:24:1 and precipitate overnight using2.5volumes ethanol/0.1M sodium acetate,pH6.5 mixture. Treat both arms with1U tobacco acid pyrophosphatase(TAP,EpicentreT19100),1h37°C, to generate5′monophosphates for linker ligation; repeat organic extraction/precipitation. This converts eligible pyrophosphate-bearing ends to monophosphate rather than defining5′OH dephosphorylation. Other RNA end classes and incomplete processing remain possible.

🟢 Methods256–282: vertis constructs cDNA libraries using its cited microRNA procedure, explicitly omitting RNA size fractionation before cDNA synthesis. Poly(A) polymerase tailing comes before RNA-adapter ligation to5′P RNA, then an oligo(dT)-adapter primer and M-MLV-RNaseH-minus reverse transcriptase. RT is42°C20min, ramp55°C, then55°C5min. High-fidelity PCR adds a library-specific4nt barcode at the5′cDNA end, targeting20–30ng/µl cDNA. Primer/catalyst identities, cycle count, adapter modifications, exact tail/landing/junction geometry and intermediate recoveries remain unassembled; no current vendor recipe is substituted.

The graph retains the physical RNA stock remainder, both actual arms and organic/precipitation fractions. Post-RT and PCR product populations retain RNA/cDNA/free-primer/reagent companions rather than asserting an unstated purification. A partial5′P RNA class is represented for ligation, while the complete adapter end and sealed product remain explicitly unassembled. It is not a source-absence claim.

## Barcodes and readout

🟢 Published Table2 originalp6: ML− has **ACTA**, ML+ has **AGGT**. The caption locates each4nt barcode in the5′end linker. These fixed sample identifiers are declared as structured features within otherwise unassembled PCR-primer contexts; four known bases are not falsely presented as whole operational primers. The table reuses some codes for other sample conditions, so they are not asserted globally unique across all study libraries.

🟢 The original pair is sequenced using Roche454FLX with variable-length reads, not an Illumina index-read mask. Source maximum mapped-read observations do not define a fixed sequencing cycle count. These sample tags are not i5/i7 sample-sheet indexes. Whole sequencing adapters/primers and read landing geometry remain unassembled.

🟢 The companion operon libraries in Methods287–296 start from total RNA across sample conditions and deliberately fragment RNA95°C5min before polyadenylation/linker/cDNA preparation. Their PCR primers follow Solexa amplicon instructions; sequencing is76single-read cycles on GenomeAnalyserII. That fragmenting Illumina route is separate from the unfragmented TEX/mock454 pair and is not fabricated as another stand-alone method here.

## Deferred work

Obtain the full published article/Methods; geturl.txt records its actual publisher URL. Read the directly cited historical vertis microRNA library original and applicable vendor specifications in a later depth-one pass to establish whole modified adapters, RT/PCR primers, source enzyme identities and physical library products. Relevant source identity/version must be checked before importing any constructor. Other supplement figures/tables and Table17 full analytical scope remain unread. No missing-sequence, reagent or layout claim is accepted. A complete model/public schematic is not claimed.
