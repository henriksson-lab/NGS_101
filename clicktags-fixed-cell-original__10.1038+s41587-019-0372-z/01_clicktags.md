# Highly multiplexed single-cell RNA-seq by DNA oligonucleotide tagging of cellular proteins

Defining DOI: https://doi.org/10.1038/s41587-019-0372-z

## Sources, versions and actual model boundary

🟢 Gehring etal., Nature Biotechnology38:35–38(2020), published online23December2019, DOI10.1038/s41587-019-0372-z. Published supporting-information identity/title/DOI checked. Read SupplementaryFigure1a–c and caption in MOESM1 (original renderedPDFp2), sample-tag interpretation captionsFig2–3, and MOESM4 modified oligo inventory. Publisher calls MOESM4 “SupplementaryTable4”; the workbook's own header says “SupplementaryTable2”. These are document labels, not a sequence-version substitution. Mandatory scraper was run; BC21/fullmodification verified in the workbook. MOESM3 is marker-gene results, not the primer table.

The published main PDF request returned a subscription-preview HTML page; full main Methods were not read. CaltechAUTHORS recordm07dw-3jq86 labels its315333.full.pdf “Submitted”. That superseded preprint is excluded rather than used as if it were a final original. No later live-cell ClickTags implementation is substituted.

## What the published supplement actually supports

🟢 SupplementaryFig1a/b labels yeast in a one-pot two-step reaction with NHS-TCO and MTZ-Cy5. SupplementaryFig1c is a DIFFERENT activity assay: MTZ-activated DNA ClickTags react with TCO-Cy5 and products are separated by PAGE. Lanes1–12 are3′amine-modified,13–14 are5′amine-modified. The local graph models only the latter selected IEDDA coupling class. It does not merge fluorescent yeast-labeling and DNA-activity assays into an invented full fixed-cell tagging recipe.

The chemical state retains DNA while coupling a Cy5 payload. No base deletion, RNA conversion, specific buffer/time/yield or copper catalyst is supplied. Unreacted dye/DNA remain in reaction; final PAGE separation is source-reported but its physical product partition is not assembled here. The assayed DNA sequence and sequence-to-lane mapping remain unassigned.

## Independent source oligo inventory

🟢 BC21 is printed as `/5AmMC6/TCGTCGGCAGCGTCAGATGTGTAAGCAGTTACAGBAAAAAAAAAAAAAAAAAAAAAAAAA`. The5′amine-C6 modification is preserved. The constrainedB is retained; it is not silently converted to an unrestricted N or a UMI. Other table entries use3′AmMO; these placements are not interchangeable.

🟡 The8base variable region `AGCAGTTA` is interpreted as the sample-tag window by comparison of the published oligo family and sample-specific tag captions. This field interpretation is inferred in the constructor and marked on the private preview. The BC21 inventory is not asserted to be an actual lane in Fig1c. No full capture/read geometry or sequencer index is inferred from a matching canonical motif.

## Deferred originals and lineage

Full published main Methods, exact conjugation/activation and wash conditions, oligo-to-experiment mapping, fixed-cell capture/RT/library preparation, and complete run/primer geometry remain unresolved. Relevant published documents have not all been read; no missing-information claim is made. The primary abstract establishes fixed-cell protein-directed DNA sample labeling, distinct from MULTI-seq's noncovalent membrane lipid anchor. Exact preparation ancestry requires the unavailable main Methods; no CITE/REAP or membrane recipe is imported.
