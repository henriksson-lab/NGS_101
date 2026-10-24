# TLA: targeted locus amplification

De Vree et al. Targeted sequencing by proximity ligation for comprehensive variant detection and local haplotyping. Nature Biotechnology 32:1019–1025 (2014). [Published article](https://doi.org/10.1038/nbt.2959).

## Sources and scope

The published main article is available as an [author-institution PDF](https://www.hubrecht.eu/app/uploads/2017/11/DeLaat_Key_2014_DeVree_Targeted-sequencing-by-proximity-ligation-for-comprehensive-variant-deletion-and-local-haplotyping.pdf). The published Supplementary Text and Figures is available from [Springer's supplement](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnbt.2959/MediaObjects/41587_2014_BFnbt2959_MOESM1_ESM.pdf). Both originals and extracted text are external to the checkout; see the reference manifest.

The bounded review covers main Online Methods “TLA sample preparation” and “Bioinformatics: Sequence alignment”, Supplementary Figure 1 and the generic tissue-culture Supplementary Protocol through TLA PCR. This documentation pass directly checked main text lines 1097–1123, the Supplementary Figure 1 caption at lines 87–104, and the Supplementary Protocol heading. Detailed protocol reading was recorded by the initial curator; it is not an independent full-source pass.

The primer scraper was run on both originals and returned 66 candidates. Hits include study-specific primer-table entries; no source pair is selected in this generic locus model. Finding those entries does not assign their sites on an unspecified locus. Full target maps and target-specific applications remain outside the declared scope.

## Molecular definition

🟢 The published generic sequence of molecular classes is crosslinked genomic material, primary restriction fragments, proximity-ligated chimeric DNA, released DNA, secondary restriction fragments, circular chimeric templates and inverse-PCR products. Physical proximity is distinct from a covalent DNA junction before ligation.

🟢 The primary NlaIII recognition motif is CATG; the secondary NspI class RCATGY includes that core. The source explains how the secondary digestion retains more captured-fragment interiors. The module checks the declared motif relation without predicting exact cuts on an unspecified locus.

🟢 Source TLA primers are distinguished from the adapter-tailed primers in the Figure 1 comparison to 4C. TLA PCR products are subsequently randomly fragmented and prepared for sequencing using manufacturer protocols. The module therefore does not inherit a 4C sequencing endpoint or select a Nextera/TruSeq preparation recipe.

🟡 All example fragment order, circle closure positions, anchor primer boundaries and library ends are schematic inferred role segments. A target-specific primer pair and map are needed for base-resolved amplicons. No oligo sequence or modification is transcribed for an unselected target. The combined PCR-product population represents source pooling without inventing independent molecular substrates.

## Sequencing and model boundaries

🟢 Main Methods names Illumina MiSeq, HiSeq 2000 and HiSeq 2500 and paired-end acquisition. In the general analysis the ends were mapped independently because of sequence reshuffling; this is computational handling, not removal of a physical read.

Exact library kit, adapter/index geometry, sequencing primers and read-cycle lengths are unassigned in this generic graph. This is an explicit scope boundary, not a claim of confirmed published absence. The page declares paired Read 1/Read 2 acquisition with a sequencing-plan limitation. No ready-to-copy sample-sheet index is inferred.

Molecular validation checks the graph's declared relations and source scope. It does not establish ligation frequencies, yield, locus coverage, variant sensitivity or source completeness.

## Deferred source work

The main Methods cites Simonis et al., Nature Genetics 38:1348–1354 (2006), reference 10, and Splinter et al., Methods 58:221–230 (2012), reference 37, for predecessor procedures. These have not been followed in this bounded pass and are not judged irrelevant. The manufacturer's library-preparation document cannot be selected until a particular kit is assigned; an unassigned kit is not an identified unavailable document. A complete linked-source inventory and independent second reading remain pending. No missing-information claim is accepted.
