# TRIBE

McMahon et al.(2016), *TRIBE: Hijacking an RNA-Editing Enzyme to Identify Cell-Specific Targets of RNA-Binding Proteins*. DOI https://doi.org/10.1016/j.cell.2016.03.007 .

🟢 Associated accepted manuscript PMC5027142 provenance24, abstract/intro36–42, specificity limits126–132 and Methods140–150 inspected. Actual detailed Procedures and TableS1–S5/figures are gated. Final VOR procedure applicability remains unverified. Mandatory bundle scraper run returned no sequences; this does not establish absent primers. One capped pass, no absence audit.

## Molecular signal and source boundaries

TRIBE fuses an RBP with the catalytic domain of Drosophila ADAR, without its two native dsRNA-binding domains. RNA adenosine becomes inosine chemically; sequencing often reports an RNA-oriented G. It is not physical A-to-G substitution of the RNA nucleotide. The local selected-site constructor preserves RNA length and identity through A-to-I marking, and derives a selected complementary cDNA/reference-oriented report separately. It predicts no RBP specificity, exact contact position, edit efficiency or universal response to every adenosine. Source notes ADAR sequence/structure preferences and false-negative risk, so a nearby A alone is not a reaction guarantee.

Source S2 fusion expression is copper-induced24h before RNA/protein harvest. S2 RNAseq libraries use the historical IlluminaTruSeq kit; exact original extraction, kit edition/internal operations/strand/read geometry remain unassembled. The separate nascent-RNA procedure explicitly references earlier extraction and rRNA removal, followed by two rounds of oligo-dT polyA depletion. Those operations are represented only on a source-qualified nascent branch; they are not imposed on all S2 libraries. Detailed applicability and molecular workups remain explicit boundaries.

Fluorescently labeled neurons are dissected/triturated and manually sorted by glass micropipette, not a single-cell barcode/UMI operation. Their low-input RNAseq preparation is source-deferred to Abruzzi2015; no modern scRNA recipe is substituted. Distinct RBP targets, neuronal populations and enzyme variants are branches/inputs, not inflated protocol entries. CLIP comparisons are separate experiments, not extra operations in the TRIBE library.

Source genomic DNA/control RNA sequencing and removal of endogenous edits distinguish SNPs/background from introduced RNA edits computationally; these controls are not physically pooled into RNA libraries. Reads report derived base calls, while original RNA and reaction remnants remain until actual source workup can be assembled. Source kits and predecessor procedures are opaque, including any physical removal fates. No finished schematic.
