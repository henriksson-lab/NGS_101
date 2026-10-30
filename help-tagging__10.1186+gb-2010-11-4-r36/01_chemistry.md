# HELP-tagging: published EcoP15I / T7 selection library

Suzuki et al., *Optimized design and data analysis of tag-based cytosine methylation assays*, Genome Biology11:R36 (2010), DOI[10.1186/gb-2010-11-4-r36](https://doi.org/10.1186/gb-2010-11-4-r36).

Private incomplete operational graph in `tools/chem.py`. This bounded pass establishes selected source operations and actual oligo inventory, not full source completeness or accepted missing information.

## Sources read

🟢 Published PMC2884539 article (original title/DOI and publication status verified), Methods103–111 and Results37–51/Fig1; published Additional file1, TableS9 (convertedtext107474–107501), original page218 rendered and inspected. PDF and JATS are the same article, not independent sources. Mandatory scraper returned37candidates in two sequence-positive files. Large SI result tables, other validation primers and every linked historical component original were not exhaustively reviewed.

## Distinct molecular preparation

🟢 Main Results37–51 explicitly distinguishes this assay from MSCC: EcoP15I produces longer flanking tags, and the added AE T7-promoter→IVT→RT route selects against unwanted single-adapter products without the MSCC gel-extraction step. This is an actual library change, not merely substituting an enzyme or normalizing counts. HpaII/MspI are separate source channels of one method.

🟢 Source5µgH1/WA01 genomicDNA is digested in separate200µl HpaII or MspI reactions, then phenol/chloroform extracted and ethanol precipitated.3µg HpaII or1µg MspI digestedDNA feeds AE Quickligation; unused digest is an independent fraction with unresolved fate. AE ligation→AMPure→EcoP15I→end repair→single3′dA→AS Quickligation→MinElute→MEGAshortscript T7IVT→RNeasy→SuperScriptIII first-strandRT→18cyclePCR→QIAquick follows the printed order.

No discarded EcoP15I daughterDNA is inferred immediately on cleavage. IVT creates RNA while physical DNA templates and AS-AS molecules remain in the reaction; transcription preference is not physical waste removal. RNeasy's exact DNA/RNA/remnant partition stays unassembled. RT produces newDNA with originalRNA/unextended primers/reagents retained; no unstated RNA-hydrolysis reaction is inserted. Source end-repair and dA-tailing enzyme identities/conditions are not supplied by another method or a modern kit.

## Exact source oligos

🟢 Original TableS9 uses lowercase to denote T7 sequence, bold for EcoP15I sequence, `p` for phosphorylation and `x` for a phosphorothioate bond:

- AS1 `5′-ACACTCTTTCCCTACACGACGCTCTTCCGATC-x-T`.
- AS2 `5′-p-GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`.
- AE1 `5′-AcagtaatacgactcactatagggagaaggctCAAGCAGAAGACGGCATACGACAGCAG`.
- AE2 `5′-p-CGCTGCTGTCGTATGCCGTCTTCTGCTTGagccttctccctatagtgagtcgtattactg`.
- PCR PS `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
- PCR PE `5′-CAAGCAGAAGACGGCATACGACAGCAG`.
- Sequencing PSeq `5′-ACACTCTTTCCCTACACGACGCTCTTCCGATCT`.

Actual terminalT remains part of AS1; the phosphorothioate is the bond before it, not a replacement base. AS2/AE2 carry printed5′phosphates. No5′phosphate is asserted for the unmarked strands. Canonical fixed bodies are imported where identical; custom AE/PE source sequence is not silently replaced with canonicalP7. No random segment is assigned UMI or sample-index meaning.

🟡 Local caller-explicit IVT/RT reference functions distinguish a DNA template from a newly synthesized RNA complement and the next complementaryDNA copy. Their coordinate/primer registers are explicitly supplied and inferred, not an assigned operational T7 start or RT primer. Full adapter/genomic ligation, IVT start/end, RT and PCR backbone geometry remain unassembled despite knowing the oligo inventory.

## Acquisition, lineage and deferred scope

🟢 Source GenomeAnalyzer36bp single-end acquisition;27nt flanking tags and CG/EcoP15I motif filtering are source readout expectations and downstream computation, not fixed enzyme-coordinate or perfect-yield predictions. HpaII/MspI normalization does not chemically rewrite native genomic cytosines.

Explicit conceptual/component predecessor: existing MSCC (Ball2009), cited in main35–37 and Fig1/Results47–51, motivates restriction-generated ends/adapter-installed distal cutter and tag counting. The HELP library changes the cutter and adds T7IVT/RT; no MSCC constructor or whole recipe is imported, and no adapter or chemistry equality is asserted. Earlier microarrayHELP and the actual historical MEGAshortscript/SuperScriptIII/end-repair component originals remain deferred after this two-original capped pass. This is ordinary unreviewed scope, not proof information is absent.
