# CLASH

Defining original: https://doi.org/10.1073/pnas.1017386108

Kudla et al.(2011), *Cross-linking, ligation, and sequencing of hybrids reveals RNA–RNA interactions in yeast*.

## Originals read

🟢 Published article `pnas.1017386108_PMC3116431.html.txt`: abstract1121, rationale1127–1129, Fig1 caption1137, Results1142–1151 and Methods1218–1220. Article identity/title and exact SI/predecessor links checked. Original SI PDF and datasetXLS remain gated; the failed EuropePMC XML route does not make the article unavailable because its published PMC HTML is readable. Mandatory full-bundle scraper run after resumed retrieval: one non-library Northern-probe hit, checked against its source context and not used as a library primer. Preliminary2026-10-28 reading is history, not an independent second pass.

## Actual strategy and model

🟢 HTP-tagged yeast RNPs were UV-crosslinked, captured with IgG, TEV-cleaved, denatured with6Mguanidine and immobilized on nickel resin. The article explicitly reports that noncovalent protein interactions are disrupted while RNA duplexes appear stable. Linker-ligation reactions16°C or25°C allow hybrid formation. Results describe directly joined RNA arms, not an inserted bridging adapter; the remaining outer ends accept library linkers. Source high-throughput readout is Illumina–Solexa50nt. This pass does not establish exact historical adapters, ligase/end preparation, PCR or primer geometry from that summary.

🟡 The graph begins with the already-purified nickel-bound protein/RNA complex, records the conditional hybrid-forming reaction and stops unresolved. Protein, unjoined RNA and alternative ligations are not silently removed. Ligase/buffer/end details from the inaccessible SI and unreviewed CRAC predecessor remain explicit source items. The local selected-event constructor can join two explicit RNA arms only when the caller supplies compatible5′phosphate/3′OH; it does not assert those states were prepared by an uninspected source step. Joining preserves every RNA nucleotide and consumes the selected end groups. It predicts neither cellular interaction truth nor ligation yield/folding stability.

No base-pair drawing is invented for unknown study arms. No RT template-switching recipe is substituted for RNA ligation. Biological duplex assignment requires the original crosslink/capture evidence and computational/experimental controls; an observed chimera alone is not proof of a unique in-vivo interaction.

## Sources deferred and ancestry

🟢 Main Methods explicitly says CRAC on Nop1/Nop56/Nop58 was described previously, ref2 DOI10.1073/pnas.0901997106, *Identification of protein binding sites on U3 snoRNA and pre-rRNA by UV cross-linking and high-throughput analysis of cDNAs*. This supports an explicit component-lineage candidate for purification/library strategy, limited by the unread predecessor and inaccessible SI. No predecessor constructors are imported. The2013 Argonaute/miRNA CLASH application is not treated as this original yeast library recipe.

Actual blocked SI: https://pmc.ncbi.nlm.nih.gov/articles/instance/3116431/bin/1017386108_pnas.201017386SI.pdf ; dataset https://pmc.ncbi.nlm.nih.gov/articles/instance/3116431/bin/1017386108_sd01.xls . Relevant SI index also gated. Notes status, source-qualified prefix only; no missing-information conclusion or finished output.
