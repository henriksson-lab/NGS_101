# DIMPLE: pooled variant-library assembly

Macdonald CB et al. (2023).
[DIMPLE: deep insertion, deletion, and missense mutation libraries for exploring protein variation in evolution, disease, and biology](https://doi.org/10.1186/s13059-023-02880-6).
Genome Biology 24:36.

## Sources inspected

Bounded source reading, 2026-10-21; the implemented scope is the Kir2.1 branch.
Read the defining published article, PMC9951526, own DOI verified: Methods “In
silico library generation”, “Library construction”, “Sequencing library preparation
and genomic DNA extraction and data analysis”, and the landing-pad construction
scope in “Cell line generation and cell culture”; local JATS text lines 137–175.
Read the published supplementary inventory at lines 183–201.

Ran `tools/scrape_primers.py` on the retrieved source bundle. Published Additional file 6, Table S3 (MOESM6 XLSX), assigns nine Kir2.1
sublibrary primer sets; read all 36 oligo-amplification and inverse-PCR rows,
including original letter case. The public page lists these computed from
`tools/dimple.py`, together with the two source genomic-target primers.
Automatically retrieved preprint files are superseded and excluded from this pass.
The linked version-1 protocols.io method is retrieved but not yet reviewed; retrieval
does not count as reading it.

Read published Additional file 3 (MOESM3 DOCX), Figure S2 and caption.
The figure was inspected from a temporary LibreOffice PDF rendering, page 2;
the source text twin was generated with `tools/doctext.py`. The variable red
length buffers lie outside the type IIS sites and the additional four unmutated
bases lie inside the retained variant-bearing region. The schematic alone does
not assign any concrete Kir2.1 variant or complete vector sequence.

## Molecular scope

🟢 DIMPLE separates computational variant design from physical construction.
Target regions are partitioned into sublibraries. Sublibrary-specific primer
handles are outside the transferred variant-bearing sequence, and length buffering
accommodates insertions and deletions without changing the synthesis-pool size.

🟢 A chosen oligo sublibrary and its corresponding backbone PCR product are
assembled by BsaI-mediated Golden Gate cloning. Constructed sublibraries are
propagated, purified and pooled. Subsequent transfer into the experimental
landing-pad backbone is a separate assembly step; it is not the same BsaI reaction.
The article names BsmBI for this transfer and distinguishes target-specific
sequence choices. Those choices must remain local to the represented experiment.

The molecular example uses the paper's Kir2.1 assay branch, with
other target-specific variants represented only after their own source assignments
are verified. No biological effect is inferred from the designed variant class.

## Sequencing path

🟢 Source input differs by experiment: some samples start with genomic DNA after
cell sorting, others with plasmid-library DNA. Target PCR enriches the library
before Nextera XT sequencing-library preparation and sample indexing.
The paper names IDT for Illumina UD indexes and paired-end NovaSeq6000 sequencing.

Keep these concepts separate in the model:

- variant and sublibrary identity;
- synthesis/amplification handles, including buffering sequence;
- retained cloning junctions versus removed type IIS recognition sites;
- sequencing sample indexes;
- physical clone selection versus computational rejection of unexpected reads.

## Implemented graph and validation

🟢 Nine separately amplified oligo/backbone pairs feed nine separate BsaI
assemblies. Source kanamycin selection precedes pooling. Pooling preserves nine
molecular classes; it does not concatenate them. The BsmBI transfer explicitly
produces separate donor and opened-recipient classes before assembly.

🟡 `backbone_fragment()` derives the primer-written BsaI ends. An independent
`plasmid.digest()` check confirms those ends for each primer pair. This checks
known terminal sequence only: the full internal reference interval and cut sites
remain unassigned. Intended variant-side junctions are marked inferred; compatible
insert ends are not independently established by matching a backbone. No complete
cohesive-assembly product claim is made.

The cell assay and genomic recovery are an explicit unmodelled graph boundary.
The readout declares source paired-end acquisition but no programmed cycle count: a
NovaSeq SP300 kit capacity does not establish it. Nextera XT and IDT UD family
names do not establish an exact index product, length or sample value. The two
sample-index features therefore have unknown lengths and values. Sublibrary
amplification handles are structured features that are removed before the
retained construct; they are not sample indexes or retained variant barcodes.

The 46 local checks cover source transcription, source-assigned sublibrary
inventory and independent primer-end cuts. The source guard matches 38 exact
primer inputs and reports four inputs without exact molecular states. Unknown
geometry remains unresolved; passing checks do not prove a fully specified
reaction, library coverage, abundance balance or functional effect.

## Remaining work

- Assign actual Kir2.1 oligo members and the predecessor reference sequence; the
  initial SurePrint pool passage names other genes and is not their identity.
- Check the cited version-1 protocols.io original in a later bounded follow-up;
  preserve the revision identity and do not substitute its latest export.
- Verify the exact source-owned vector and example-target records before deriving
  cassette cuts, PCR boundaries, strand ends or ligation-ready states.
- The source software publication version is identified by Zenodo
  10.5281/zenodo.7574260. The validation model should check declared products,
  not recreate a general variant-design or fitness-optimization pipeline.
- No sequence-compatibility check establishes experimental library coverage,
  abundance balance, cellular phenotype or a clinical variant interpretation.
- A second independent source reader and complete linked-source inventory remain
  pending. No confirmed-missing information is claimed.
