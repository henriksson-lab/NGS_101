# RAD-seq — original restriction-site-associated DNA library

Baird et al.(2008), *Rapid SNP Discovery and Genetic Mapping Using Sequenced RAD Markers*, defining DOI https://doi.org/10.1371/journal.pone.0003376 .

## Sources and scope

🟢 Published PMC2557064 JATS identity, Results72/Fig1, Methods141 and computational-processing145 inspected. Actual main includes the adapter/primer sequences and preparative order. Mandatory bundle scraper returned20hits in the main PDF/JATS twins; pointers and phosphate/placeholder notation verified in the JATS original. RFLP validation primers elsewhere in main are not RAD library oligos.

The standard supplementary ZIP request timed out/stalled and produced no inspected supplemental original. The read main does not identify a named chemistry supplement to substitute for it; an incomplete generic archive is not evidence of absent sequences. Historical Solexa/Phusion product protocols and the directly cited Coyne2004 Y-adapter predecessor remain unreviewed. One chemistry original read in the bounded pass; no independent absence audit or public promotion.

## Source components

🟢 P1's fixed58nt body is identical to the canonical `illumina.NEBNEXT_UNIVERSAL_PRIMER`, imported rather than retyped. Its bottom fixed body is the computed reverse complement. SbfI has top `[body]xxxxTGCA` and bottom `Phos-xxxx[reverse complement of body]`; EcoRI has top `[body]xxxxx` and bottom `Phos-AATTxxxxx[reverse complement of body]`. The source's lowercase x placeholders represent designed sample barcodes:4nt for SbfI and5nt for EcoRI, not molecular UMIs. Particular experimental barcode members remain unassigned. Scene computes their antiparallel placement and preserves source bottom5′phosphate and source3′TGCA/5′AATT cohesive ends.

Whole printed P1 strings are checked against original Methods141 locally, including barcode width and phosphate, and skip when the external original is absent. The generic strict source matcher does not join lowercase x windows, so one precise P1-input exemption covers those two strands; no synthetic source transcript or global parser change. An independent reviewer compared all four whole printed P1 strands, both P2 strands and both PCR primers against the model and confirmed transcription and source phosphate/barcode-width distinctions. This is a critical implementation review, not a full original/linked-material absence audit.

🟢 P2 top is `Phos-CTCAGGCATCACTCGATTCCTCCGAGAACAA`; bottom is `CAAGCAGAAGACGGCATACGACGGAGGAATCGAGTGATGCCTGAGT`. The terminal bottom3′T extends beyond their derived common complementary window. Full ordering strands are shown separately; the computed stem window is a component binding reference, not a complete paired Y fold or sealed genomic junction. Forward PCR primer is `AATGATACGGCGACCACCGA`; reverse is `CAAGCAGAAGACGGCATACGA`. These are respectively20nt/21nt canonical prefixes, not full modern P5/P7 primers.

## Preparation and physical fates

🟢 Selected source SbfI branch: digest37°C15min,65°C20min inactivation, P1 ligation at room temperature20min,65°C20min inactivation, then physical pooling of separately tagged samples. Heat does not remove free adapters or enzymes. Random sonication follows tagging/pooling, producing material with and without P1. This differs from ddRAD's two-restriction-end selection; enzyme/species/phenotype applications remain branches of the original RAD method.

Source shearing averages500bp;300–700bp gel selection/MinElute recovery follows. Quick Blunting Kit end polishing, Quick Spin cleanup, Klenow exo-minus/dATP adenylation37°C, another purification, P2 ligation at room temperature and final purification/50µl elution are separate operations. Unknown reaction durations or historical kit internals are not filled from contemporary recipes. Adapters are ligation-ready source components; complete genomic substrate ends and sealed product junctions remain unassembled.

Five microlitres enter18cycle Phusion PCR; the independent remainder has an explicit unresolved fate. Source intended selection occurs because P1-initiated copying generates the complementary P2 reverse-primer site. Direct binding of the reverse primer to the original same-sense divergent arm is not assumed. This is the source design, not a checked whole-library polymerase relation or guarantee of background depletion. New copies coexist with original templates, non-P1 backgrounds and reaction remnants until actual final300–700bp gel recovery and10nM dilution.

## Readout and limits

🟢 Source reports Illumina Genome Analyzer/Solexa sequencing; exact acquired cycles and primer/end geometry remain unassembled. Methods145 trims reads to barcode plus26nt computationally, which is not an assigned physical run length. The discussion's paired-end future application is not an acquired paired-end run for the selected original preparation.

Source Results72/Fig1 cites Coyne2004 (DOI10.1128/AEM.70.9.5298-5304.2004) for the divergent Y design, not a complete inherited library. Earlier RAD microarray methods establish historical development but are not counted as NGS protocols or imported as equivalent recipes. The private graph retains explicit preparative steps and unresolved molecular geometry; no accepted missing-information or complete chemical-validation claim.

## Bounded public-page review, 2026-10-30

Rechecked one published main original, PMC2557064 Methods141 and Results104. The mandatory scraper reported20hits/2files; actual P1/P2/PCR orderings and barcode/phosphate notation were inspected in the main text. No new chemistry supplement or historical vendor recipe was obtained.

P1 ligation now receives actual rATP and buffer/water supplies. Quick Blunting is an opaque source-named molecular-supply interface rather than an entire kit called a catalyst. The A-tail donor remains the source-written adenine supply with chemical form unreviewed. The P2 functional ligase is constructor-inferred and marked in the preamble and step caption; no particular enzyme or dose is asserted from the historical protocol reference.

The actual5µl amplification-primer mix yields its copresent forward/reverse components. The50µl Phusion master mix similarly supplies polymerase, nucleotide and other components without adding fresh independent doses;40µl water is explicit. Known produced primer bodies now have original-source regression checks because the external-only guard excludes these inventory views. Existing whole P1 fixtures and the precise lowercase-barcode matcher exemption remain unchanged.

Literal inspected source coverage does not cover whole Y/genomic junctions, exact historical products or acquired primer/read geometry. Original `historical_readout` remains required and uncovered. No source-absence claim or public completion marker is introduced.
