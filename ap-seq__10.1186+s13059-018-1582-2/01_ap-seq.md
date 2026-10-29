# AP-seq

Defining DOI: https://doi.org/10.1186/s13059-018-1582-2

Poetsch et al. (2018), Genomic landscape of oxidative DNA damage and repair reveals regioselective protection from mutagenesis.

## Sources inspected

🟢 Published open-access PMC6284305 main, AP-site colorimetry/AP-seq and OGG1 enrichment Methods171–179, library preparation181–183, data interpretation thereafter. Publisher Supplemental PDF MOESM1 FigureS1 was rendered and visually inspected on PDFpage2: ARP tags BOTH intact AP aldehydes and AP β-elimination intermediates (bE sites), forming biotin conjugates. Publisher Supplemental DOCX MOESM2 TableS1 identifies synthetic in-vitro pulldown/qPCR oligos; those are not the genomic library adapters. Superseded preprint request failures are excluded from active chemistry evidence and are not source blockers.

## Direct AP and OGG1-AP branches

1. 🟢 Source genomic DNA is kept cold during extraction. For directAP-seq, tag5.7µg with5mMARP inPBS37°C2h, then1.8xAMPure/two70%ethanol washes. ARP is O-(biotinylcarbazoylmethyl)hydroxylamine. The reactive site class includes AP aldehydes and β-elimination products, not a single intact sugar structure only; exact genomic site coordinates/yields are unassigned.
2. 🟢 OGG1-AP is a distinct source pretreatment: recombinantOGG1,0.1µg perµg genomicDNA,NEBbuffer2/BSA37°C1h, AMPure-purify, thenARP. Glycosylase reaction/end-product distribution is not assumed quantitative or identical for all genomic substrates. This branch contains native AP signal in addition to newly generated reactive products; it is not modeled as a perfect subtraction measurement of8-oxoG.
3. 🟢 Only AFTER tagging/cleanup, Covaris-fragment toward300bp. Reserve30µL of130µL as the matched sequencing input; this input has undergone tagging. It is not replaced with untagged genomic DNA from another method.
4. 🟢 Prewash MyOne streptavidin beads3times, resuspend2MNaCl/TE, capture genomicDNA10h at room temperature, wash3times1MNaCl/TE and resuspendTE. No antibody or strand-selective alkali elution is imported.
5. 🟢 Treat BOTH enriched material and input withPreCR. The source then AMPure-purifies input and pulldown supernatant, and recombines recovered pulldownDNA with the beads. Thus both soluble and bead-associated material are retained for library preparation; PreCR is not assumed to leave every adduct or molecule bead-bound. Exact tag/repair/end chemistry and partition remain unassembled, not claimed absent.
6. 🟢 Source KAPA kitKK8504, described as125bp paired-end ChIP-seq preparation, thenHiSeq2000 rapid/high-output runs; combine read data. Adapter/index/PCR polymers, original kit revision and exact primer geometry remain unassembled. No synthetic-controlM13 primer is substituted for genomic reagents.

## Interpretation and boundaries

Enrichment is quantified relative to the matched input over genomic bins; this is a regional damage signal, not a lesion-encoded single-base call. Computational alignment/duplicate/repeat analyses are distinct from molecular purification. The colorimetric assay is source-described as semi-quantitative and was not used to normalize sequencing. DirectAP andOGG1-AP patterns are compared with the source's caution about relative quantification and fragment distributions.

The private graph preserves tagging-before-fragmentation, both source branches, input allocation, PreCR-associated soluble/bead fractions and recombination. It cannot validate full PreCR/ARP chemistry, quantitative enrichment, exact historical library products or run/index geometry. It has an uncovered historical kit/source-assembly item and no finished marker. OG-seq appears as a comparison/reanalysis source, not evidence that its oxidative tagging/complement-release recipe is inherited. Direct reagent originals and independent review remain deferred; no missing-information claim is accepted.
