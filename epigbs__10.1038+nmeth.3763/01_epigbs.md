# epiGBS

van Gurp et al., *epiGBS: reference-free reduced representation bisulfite sequencing*, Nature Methods13,322–324(2016), [DOI10.1038/nmeth.3763](https://doi.org/10.1038/nmeth.3763).

This is a capped source-qualified preparative note/private incomplete graph, not a finished schematic or an accepted missing-information audit. The defining main Methods were not obtained. No earlier GBS constructor is imported.

🟢 Published Supplementary Figures MOESM158, FigS1 originalp1 and captionp1–2, FigS3p4 and FigS9p11 were inspected. Both FigS1 and FigS9 were rendered from the original PDF: image-based oligo rows have no bases in the text extraction. The mandatory scraper initially returned zero hits; that did not establish absent sequences. Supplementary sample-description workbook MOESM159 Species_level header/first sample rows were inspected only for its actual identity/scope. Other result workbooks remain uninspected, not declared chemically irrelevant after a full audit.

🟢 FigS1 attributes adapter design to Elshire2011GBS and explicitly distinguishes child B-adapter barcodes. It shows methylated complementary ordering oligos, unphosphorylated adapters, PstI-compatible TGCA cohesive ends and historical Illumina paired-end PCR primers. The local model selects A1/B1 `AACT` as inventory components only; no study sample, well, pooling ratio or frequency is assigned. A-barcode lengths vary4–6nt, B barcodes shown are4nt. They are sample identifiers, not UMIs. The Csp6I branch uses AT overhangs and its own listed barcodes; it remains within this protocol, outside the selected PstI model.

Every C in the source ordering oligos is 5-methyl-C, including barcode and cohesive-end C. Sequence spelling alone does not express that chemistry; graph inputs retain a separate typed methylated-adapter inventory. FigS1b full strand equations, all5′→3′, are:

- A top: canonical Read1 body + A barcode + TGCA.
- A bottom: reverse-complement A barcode + reverse-complement Read1 body.
- B top: B barcode + reverse-complement historical common body.
- B bottom: historical common body + reverse-complement B barcode + TGCA.

Bodies, complements and historical PE-PCR primers are computed from `illumina`; local `Scene.anneal()` checks both duplexes rather than hand-typing paired columns. Input5′OH and3′TGCA adapter geometry are explicit. The complete genomic precursor and all covalent ligation products remain unassembled. All ordered AA/AB/BA/BB classes persist; no perfect AB selection is asserted.

🟢 FigS9 describes nicks between genomic3′ ends and unphosphorylated adapter5′ ends. Before bisulfite treatment, DNA polymerase I removes adapter nucleotides with5′→3′ exonuclease activity and replaces them by extension from the genomic3′OH, using methylated C replacement nucleotides, proceeding through the adapter end. This is a source-selected repair class, not a prediction of repair efficiency, exact enzyme kinetics or exclusive products. Original unrepaired/partial molecules, excised nucleotides and enzymes remain physical companions until an actual workup. Synthetic adapter methylation does not rewrite native genomic methylation.

🟢 FigS3 establishes actual paired-end fragment acquisition, including overlapping and non-overlapping read classes. Its10-N concatenation is computational; it does not create a physical DNA insert. FigS1 lists historical Illumina PE-PCR primers. The graph keeps bisulfite conversion and PCR descendants separate from original adapter marks; PCR copies are not assumed to retain5mC chemistry. Detailed cleanup, aliquot/QC/pooling, PCR conditions and actual read-cycle count await main Methods, rather than a borrowed GBS recipe.

🟡 The local repair helper preserves source adapter identities/native marks across a conditional selected repair event. It does not simulate complete strand chemistry or conversion probabilities. All new workup boundaries are explicitly unreviewed rather than invented purification procedures.

The standard source tool obtained published SI and result workbooks, but returned no open full defining article. Direct predecessor GBS original, generator documentation, complete main Methods and run/library originals remain ordinary deferred targets. No absent-sequence, absent-condition or publication-error conclusion is accepted. See the source checklist; neither a capped pass nor independent implementation review supplies the two-reader absence gate.

Independent critical review re-rendered FigS1 from the original at2600pixels, checked all six full ordering/PCR rows and FigS9 chemistry, and approved only the two B-adapter text-match exemptions. Four A/PCR strands remain source-guard checked. This selected transcription/implementation review is not an absence audit.
