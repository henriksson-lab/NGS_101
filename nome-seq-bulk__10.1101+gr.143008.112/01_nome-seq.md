# NOMe-seq — bulk defining implementation

Kelly et al., *Genome-wide mapping of nucleosome positioning and DNA methylation within individual DNA molecules*, Genome Research (2012), [defining article](https://doi.org/10.1101/gr.143008.112).

## Sources inspected

🟢 Published original PMC3514679, DOI and printed Genome Research identity checked. The manifest records an open-access publisher article, not an author manuscript. Read Results on context interpretation (JATS55–75), Methods “Nucleosome footprinting” (119), “Library construction and sequencing” (123), and alignment/context/read filtering (127–135), plus the identities of the directly cited preparative predecessors in its reference list.

The mandatory primer scraper ran on both extracted main-article representations and returned no long oligos. That is not an absence conclusion. The source explicitly states locus-specific PCR primers are available upon request; those locus-specific cloning assays are not this whole-genome library. The exact whole-genome library polymers remain unassembled pending linked-source review.

The fetched Europe PMC supplementary zip contains only main-article JPG/GIF figures, not an inspected chemistry supplement. The publisher's supplemental-material index at `https://genome.cshlp.org/content/suppl/2012/10/04/gr.143008.112.DC1` returned HTTP403 through the standard declared fetcher. The original supplement checklist remains blocked; no specific unknown attachment is invented. See [retrieval manifest](ref/MANIFEST.md).

## What is measured

🟢 M.CviPI exposure marks accessible GpC cytosines, while cytosines at endogenous CpG sites provide a second signal. The native genomic molecule carries both signals. Accessibility is source-described biochemical interpretation, not a computed rate or a binary classifier of all nucleosome/protein occupancy.

🟢 The source separates GCH accessibility, HCG endogenous methylation and ambiguous GCGs; H means A/C/T. GCGs are excluded from its interpretation because either methylation origin can produce the same protected C readout. The main Results also reports limited off-target CCG methylation. Most analyses retain HCG including CCG; the within-read comparison instead uses WCG, where W means A/T. The model does not silently impose one universal exclusion rule on every source analysis.

🟡 The local `MethylatedReference` accepts explicit reference DNA and caller-supplied endogenous/added m5C positions. It validates cytosine coordinates and distinguishes GCH, WCG, ambiguous GCG and source-qualified CCG off-target events. It does not assign methylation to every accessible motif or infer the origin from a sequencing C. Its selected idealized bisulfite reference contains physical U at unmodified C; newly copied DNA has T at that reference position. The helper is neither an efficiency simulator nor a complete modified-cytosine model. Its calculation fixtures are synthetic, not study molecules or source oligos.

## Source IMR90 route

1. 🟢 Begin with isolated IMR90 nuclei. The article refers directly to Miranda et al.2010 for nuclear isolation; that original remains unreviewed.
2. 🟢 Incubate with the first M.CviPI/SAM dose, then add the fresh enzyme/SAM boost. These are separate timed source operations, not a single inferred equivalent incubation. SAM is a molecular methyl donor, distinct from the enzyme catalyst. The graph keeps marked and unmarked populations and unassigned carryover; it predicts no occupancy or labeling fractions.
3. 🟢 Stop the reaction and extract DNA. Exact source stop/extraction chemistry remains unassembled. Physical extraction fractions are separated without assigning quantitative recovery or perfect removal of a particular molecular species.
4. 🟢 Fragment extracted DNA, end-repair with the source Epicentre reagents, and ligate methylated Illumina adapters **before** bisulfite conversion. The fragmenting mechanism, historical adapters and exact ligation-ready ends remain unreviewed; neither a sonication method nor a newer adapter recipe is substituted.
5. 🟢 Apply the named Zymo EZ DNA Methylation kit, then PCR. The historical intervening workup is an explicit unresolved model boundary; no undocumented wash or disappearance of molecular products is asserted. Physical U, protected selected m5C and newly copied DNA C/T are distinguished.

This private graph follows the source operations but retains unresolved molecular relations and an uncovered linked-source inventory item. It has no finished-protocol marker and remains notes-only.

## Readout

🟢 IMR90 used HiSeq2000 single-end {{= nome.READ1_CYCLES }}. The two glioblastoma samples instead had separate paired-end runs of different source lengths. Those are not merged into a paired-end IMR90 plan. Run-primer landing sites, index members/orientation and complete final-library strands remain unassembled.

Context filtering, genotype-strand checks, duplicate collapsing, trimming nonconverted read ends and genomic alignment are computational interpretation. They do not physically remove bases or genomic molecules from the preparation. The source's same-molecule interpretation does not make the bulk preparation a physically single-cell protocol.

## Deferred sources and questions

- Publisher supplemental originals behind the identified index, including the actual condition and off-target evidence panels. No attachment identity is guessed from the failed index.
- Miranda2010 nuclear isolation/single-promoter assay, DOI10.1002/0471142727.mb2117s89.
- Directly cited library predecessors Lister2009, Kelly2010 (DOI10.1016/j.molcel.2010.08.026) and Berman2012 (DOI10.1038/ng.969), plus exact historical reagent manuals. These citations do not authorize constructor reuse or prove identical adapters across protocols.
- Full oligo/end/library/read geometry, residual modified-base classes, source workups and an independent original/linked-material pass. No confirmed missing-information claim has been accepted.
