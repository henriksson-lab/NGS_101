# CRAC

Granneman et al.(2009), *Identification of protein binding sites on U3 snoRNA and pre-rRNA by UV cross-linking and high-throughput analysis of cDNAs*. DOI https://doi.org/10.1073/pnas.0901997106 .

## Actually inspected original

🟢 Published PMC2688437 HTML original: identity/abstract1105, rationale1123, Fig1caption1131, detailed library discussion1143–1149, Solexa branch1171–1179 and Methods1204–1212. The main article explicitly includes both Sanger-clone and direct high-throughput Solexa readouts. The latter establishes an NGS protocol distinct from subsequent CLASH RNA–RNA joining. Mandatory complete-bundle scraper returned zero long-sequence hits in the available main article; this does not establish that TableS5 oligos are absent. The actual linked SI index and PDF returned404 under standard retrieval; neither was inspected.

## Linker versions and noncanonical nucleotide blocks

🟢 Original linker discussion describes a3′DNA linker that is5′adenylated and terminates in blocked ddC; the later5′linker is a DNA–RNA hybrid with an inverted dideoxythymidine at its5′end, designed to suppress concatamers. These are nucleotides with end/backbone orientation properties, not non-nucleotide spacers. The local typed inventory preserves cytidine versus thymidine and normal versus reported-inverted orientation. It refuses to flatten an unassembled inverted/block-containing oligo into ordinary linear DNA.

Those descriptions accompany the initial/Sanger workflow. The article then explicitly states that linkers were modified for compatibility with Illumina Solexa. The exact Solexa sequences/end topology in TableS5 are unreviewed, so the graph does not silently transplant the Sanger linker inventory into the NGS library. It carries the source-known ligation order using distinct unassembled Solexa adapters. Shared unknown-end display factories that assert “not reported” are not used.

## Source-known NGS-route operations

The graph chooses the in-vivo crosslinked HTP-RNP route, beginning at the TEV-eluate boundary. It does not interchange this with the in-vitro branch's later UV exposure. Source partially digests RNA with RNaseA/T1, adds6Mguanidine, immobilizes His-tagged protein/RNA on nickel and washes. RNA2′/3′phosphates are removed, the3′linker is ligated on beads, RNA5′ends receive PNK/[γ-32P]ATP labeling, then the5′linker is ligated on beads. Molecular debris remains until a physical selection/wash; no complete RNA recovery or pure RNP population is predicted.

Source imidazole elution, TCA precipitation, NuPAGE/nitrocellulose selection and proteinaseK release follow. Residual peptide/adduct chemistry and detailed RNA recovery remain unassembled. RNA is then reverse-transcribed/PCR-amplified, gel-purified and subjected to source Solexa sequencing. The known operations are carried through that reported readout; exact RT/PCR reagents, primers, NGS-size selection and acquired read cycles remain uncovered source items. The mapped RNA lengths22–34nt are not assigned as run-cycle lengths. RT deletions/substitutions are inferred copying signatures, not physical deletion of bases from recovered RNA.

## Sources, ancestry and completion

Actual inaccessible SI originals: https://pmc.ncbi.nlm.nih.gov/articles/instance/2688437/bin/supp_106_24_9613__index.html and https://pmc.ncbi.nlm.nih.gov/articles/instance/2688437/bin/0901997106_0901997106SI.pdf . The failed EuropePMC XML route does not make the main article unavailable: its published HTML was read. Modern CRAC protocols, Sanger primers and later CLASH/HTP recipes are not substitutes for the original Solexa component version.

CLASH's defining Methods explicitly cites this CRAC source for purification/library strategy, supporting a scoped CRAC→CLASH edge. Main CRAC discusses earlier CLIP, RL5 and preadenylated-DNA ligation sources; they remain directly linked follow-ups, not imported recipes. The graph has an unassembled read plan and explicit uncovered SI/source inventory, no finished-publication marker and no absent-information claim.

## Bounded physical-reagent corrections with original SI hold retained

The same published main supports molecular guanidine-HCl and imidazole
contributions. TCA is now a molecular precipitation reagent with a physical
cleanup operation, rather than a selection matrix. RT and PCR supplies expose
aggregate copresent nucleotide/buffer functional classes without asserting
stock packaging, exact doses, enzyme identity or any uninspected Solexa recipe.

Nine focused checks pass. The actual `unreviewed_original_si` item and uncovered
NGS RT/PCR/workup/primer scope remain unchanged. Original SI index/PDF HTTP 404
records are prior availability history; no retry or new access claim was made.
The main remains an available published original, so the hold is specific to
the unresolved SI and modified Solexa components, not to all source evidence.

The Sanger DNA/ddC and DNA–RNA/inverted-ddT inventory remains separate from
the explicitly modified Solexa linkers. No Sanger body, block or size window
is assigned to the NGS route. RNA crosslink-associated copying errors are not
physical RNA mutations/deletions. There is no PAGE_FILE builder, complete-page
claim, source item exclusion or accepted missing-information audit.
