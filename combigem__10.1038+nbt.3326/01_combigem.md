# CombiGEM: barcoded combinations in human cells

Wong ASL et al. (2015).
[Massively parallel high-order combinatorial genetics in human cells](https://doi.org/10.1038/nbt.3326).
Nature Biotechnology 33:952–961.

## Sources inspected

Bounded source readings, 2026-10-21; a scoped assembly/readout model is now implemented.

- Published article author manuscript, PMC4785103, own DOI verified. Read Online
  Methods “Creation of the barcoded single miRNA precursor library”, “Pooled
  combinatorial miRNA library assembly for high-throughput screening” and
  “Barcode sequencing”: local JATS text lines 109–119 and 129–135.
- Published ESM `MOESM33` contains Table 1 on original PDF pages 3–4 and
  the Table 13 construct inventory beginning at page 60. Read Table 1's first
  three primer pairs in the rendered page 3 and follow their text. Read Table
  13's construct headings and storage/destination inventory; this is not a
  reading of the actual vector maps. `MOESM32` Figure 1 concerns expression
  validation and does not substitute for the original main assembly figure.
- Ran `tools/scrape_primers.py` across the retrieved source bundle. The Methods
  passages above establish the assembly order and readout, rather than a scraper hit.

This entry concerns the 2015 human-cell implementation. It does not inherit the
2014 bacterial workflow or the later CombiGEM-CRISPR chemistry by name alone.
This is not a complete source or missing-information audit.

## Molecular distinction

🟢 Each input precursor has a source-assigned eight-base identifying barcode.
The single-insert PCR product is described in this order:

`BglII – precursor – BamHI – EcoRI – barcode – MfeI`.

🟢 Storage-vector inserts are released with BglII/MfeI and ligated to the
BamHI/EcoRI-opened destination. The paper explicitly pairs BamHI with BglII
and EcoRI with MfeI as compatible sticky ends. Repeating the assembly joins
precursors at one region and concatenates their associated barcodes at another.

🟢 Table 1 reverse primers establish a {{= len(combigem.SPACER) }}-base `{{= combigem.SPACER }}` spacer on the construct
strand between BamHI and EcoRI. It is also present in the barcode-PCR forward
primer. The Methods' shorthand order above does not say that these sites touch.

🟡 The declared fragment model derives mixed BamHI/BglII `{{= combigem.PRECURSOR_JUNCTION }}` and MfeI/EcoRI
`{{= combigem.BARCODE_JUNCTION }}` junctions, while retaining the insert's central BamHI–spacer–EcoRI cassette
for the next round. Appending a precursor prepends its barcode: a precursor path
A→B→C has barcode order C→B→A. Circular traversal is essential; the barcode tract
crosses the model's drawing origin. This is a strand-order derivation for the
declared fragments, not proof of a whole-plasmid digest or experimental association.

The public schematic illustrates the first three Table 1 source assignments:

| Member | Construct-strand barcode | Barcode in written reverse primer |
|---|---|---|
| {{= combigem.MEMBERS[0].precursor }} | {{= combigem.MEMBERS[0].barcode }} | {{= combigem.MEMBERS[0].reverse_primer_barcode }} |
| {{= combigem.MEMBERS[1].precursor }} | {{= combigem.MEMBERS[1].barcode }} | {{= combigem.MEMBERS[1].reverse_primer_barcode }} |
| {{= combigem.MEMBERS[2].precursor }} | {{= combigem.MEMBERS[2].barcode }} | {{= combigem.MEMBERS[2].reverse_primer_barcode }} |

The construct-strand values are reverse complements of the written primer windows.
These identify one illustrative trajectory through the pooled rounds. They do not
enumerate all combinations or establish recovered genomic identities.
Barcode-to-precursor identity is distinct from the sample multiplexing index.

## Sequencing path

🟢 The source recovers barcode amplicons from pooled plasmid DNA and from genomic
DNA after the cell-based assay. Barcode PCR adds Illumina anchor sequences and a
sample-specific eight-base index. The barcode amplicons are purified, quantified
and pooled for HiSeq sequencing with custom primers.

🟢 Source-transcribed oligos, 5′→3′, from Methods lines 133–135; printed spacing
is removed here without changing the base order:

| Role | Sequence |
|---|---|
| Barcode PCR forward | {{= combigem.FORWARD }} |
| Barcode PCR reverse | {{= combigem.REVERSE }} |
| CombiGEM barcode sequencing | {{= combigem.BARCODE_PRIMER }} |
| Index sequencing | {{= combigem.INDEX_PRIMER }} |

The `N` window above represents the sample index and is not a known fixed sequence.
Canonical Illumina regions should be imported in the model, with their assignment
to this protocol justified by these source oligos.

## Current model and validation

The protocol owns its storage/destination declarations, selected-site release,
three iterative assemblies, separate plasmid/genomic barcode-PCR inputs, target
band selection, concentration measurement and sample pooling. The generic
cohesive operation supplies no vector, enzyme choice, part identities or barcode
association defaults. It checks only supplied fragment ends and circular product
sequence. Unknown precursor/backbone cores keep this relation conditional.

Source-association bookkeeping checks that every declared member's barcode remains
unchanged and appears in reverse precursor order. Precursor bodies remain labelled
roles: this is not sequence proof of precursor identity, absence of recombination,
or fidelity of barcodes in recovered genomic DNA. The genomic input starts after
the reported assay, with member identities unassigned. No delivery or integration
reaction is manufactured between the plasmid and genomic branches.
Repeated member identities are permitted by the source's Cartesian combinations;
the declaration checks preserve their multiplicity rather than deduplicating them.

The pooled-reaction ligase identity remains unassigned. T4 DNA ligase is named
for a different individual-vector construction passage; this alone does not
establish the reagent used in the pooled assembly.

The library ends derive from the written PCR primers; their full source-template
positions remain unassigned. PCR geometry therefore records that unresolved
scope explicitly. The 35-base custom index primer has a computed 20-base 3′
match supplied by reverse PCR. Its remaining 15-base 5′ region lies in an
unassigned interval, so a full 35-base match is not claimed. Primer sites and
potential read starts are computed separately for the two sample classes; the
final pool retains independent molecular classes rather than concatenating them.

The concentration step is a measurement of aliquots with unchanged retained
library states. Calibration chemistry and qPCR product molecules are not modelled.
Actual sample indexes, pooling proportions and acquisition cycles remain open.
The source's 359-bp band description does not set read cycles or establish that
all assembly orders have the same length. The exact graph-input oligos pass
the source-text transcription guard; the first three Table 1 barcode/spacer
prefixes have an additional source-dependent regression check.

## Remaining work

- Inspect the remaining Table 1 rows and original main graphical assembly diagram;
  obtain the cited storage/destination-vector originals and verify their identities.
- Add the initial genomic-precursor PCR/storage-library creation from the remaining
  source primer pairs, with their own template and vector geometry.
- Assign full source maps before claiming complete digests, restriction-site
  uniqueness in unknown part bodies or exact genomic primer coordinates.
- Verify the pooled-reaction ligase identity within its applicable source scope.
- The current concentration step abstracts measurement; detailed qPCR molecular
  products and calibration remain outside the graph.
  Plasmid-library QC and genomic recovery are separate source-backed paths.
  Keep delivery and cellular phenotype procedures outside the initial cloning and
  sequencing scope; do not turn an unmodelled cellular stage into exact chemistry.
- Verify sequencing acquisition lengths against the published readout sources;
  an amplicon's length does not specify read cycles.
- Experimental representation, packaging, integration or barcode-association fidelity
  do not follow from sequence-compatible assembly declarations.
- A second independent source reader and complete linked-source inventory remain
  pending. No confirmed-missing information is claimed.
