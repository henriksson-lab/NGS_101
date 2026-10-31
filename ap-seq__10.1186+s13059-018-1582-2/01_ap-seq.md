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


## Bounded publication adjudication — 2026-10-30, batch 11

🟢 DirectAP and independent OGG1-AP genomic preparations are supported. Tagging precedes Covaris;30ul/130ul matched input and remaining100ul pulldown remain paired. PreCR-treated soluble pulldown material is purified and recombined with bead-associated material before library preparation. Lesion enrichment is regional relative to tagged input, not a lesion-specific basecall.

🟡 Original required historical_reagent_and_library_assembly remains uncovered. It explicitly requires historical PreCR/KAPA original reaction details and full library assembly. Existing originals establish actual genomic selection, PreCR handling/recombination and KAPA KK8504/HiSeqPE125 handoff, but not the required historical reaction/polymer/read assembly. A sequencing terminal or narrower label cannot satisfy this original obligation. No source acquisition attempted and no information is claimed absent.

Source scopes inspected: Published PMC6284305 main XML/PDF: Methods163–183: validation controls versus genomic preparation; ARP before shearing, tagged input, OGG1 separate branch, PreCR and pooled bead/supernatant DNA, KK8504/HiSeqPE125; Published MOESM1 PDF: Actual FigureS1 original PDFp2: intact AP and beta-elimination aldehydes with ARP oxime/biotin classes and capture diagram; Published MOESM2 DOCX: Actual document.xml TableS1 six rows; guanine/8oxoG damage controls and M13 qPCR primers, no genomic-library primer substitution.

Bounded remaining physical/model obligations: Document required historical PreCR/KAPA scope hold; retain both actual DNA fractions. Actual ARP is a physical covalent tagging reagent and must be K.MOLECULAR/R.STRUCTURAL, not only condition. Actual PreCR supply is enzyme/reagent inventory, not just K.CONDITION. OGG1 BSA is physical component, not only catalyst label. Prewashed bead stock and PreCR/KAPA supplies are currently shared across independent source branches without actual allocation/provenance. Any future extension should allocate true source branch stocks and QC/loading remainders without inventing study sample counts. Do not substitute synthetic in-vitro/qPCR TableS1 polymers for genomic KAPA library reagents. Whole library/input lesion states may remain unknown, but original required source obligation cannot be excluded.

The existing private preview and normal publication gate remain unchanged. No required source item is excluded or renamed, no unknown occupied molecule is invented, and this is not a confirmed absence claim. No source acquisition or architecture expansion was undertaken.
