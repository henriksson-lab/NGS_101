# spinDrop

**Defining source:** doi:10.1038/s41467-023-40322-w.

## Molecular boundary

🟢 Cells or nuclei are co-encapsulated with inDrop-v3 microgels, fluorescence-sorted, heat-treated and sent through a second device that picoinjects RT mix.

## Evidence boundary

🟢 Sorting and picoinjection are explicit physical steps; the final identifiers follow the inDrop-v3 molecular scheme.

## Oligo source review

🟢 The published article (PMC10409775 Methods lines281–293) specifies
inDrop-v3 bead preparation, UV primer release, and RT after picoinjection.
It states61/8/8/14 acquisition. The superseded preprint is retained as history
and is not an active source for this migration.

🟡 The spinDrop molecular model is currently role-level. To resolve its
primer-binding geometry, the cited inDrop preparation protocol (Zilionis et
al., doi:10.1038/nprot.2016.154, ref. 28), the cited v3 barcoding source
(Briggs et al., doi:10.1126/science.aar5780, ref. 29), and the exact bead
implementation used here must be compared; an inDrop sequence from another
version cannot establish the bases of these beads. The article's available
workbook supplements are cell counts, device fabrication and cost tables,
not an oligo table.

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC10409775 JATS265–294, microfluidic encapsulation and RT picoinjection; article author heading.

Corrected citation to DeJonghe etal. Added25-minute23°C lysis/7-minuteUV primer release before70°C10-minute treatment; clarified bead primers precede picoinjected RT reagent. Recorded61/8/8/14 acquisition.

Deferred: Detailed inherited inDrop-v3 IVT/library construction and recovery, selected bead oligos and RT-variant conditions; quantitative barcode-like workbook hits not treated as molecular oligos. This capped pass does not establish a complete source review or publication-wide absence.

## Local chemistry ownership migration — 2026-10-18

The local graph now declares the source's physical picoinjection step, followed
by SuperScript III RT, separately from the inherited library route. The bead
primer was already released by UV before picoinjection; it is not newly added
by the RT solution. Published Methods line293 gives the selected 50 °C,2 h
RT followed by70 °C,15 min and the 61/8/8/14 acquisition, now carried by the
graph read plan. The Maxima H− and PEG options remain separate alternatives.

The whole downstream inDrop preparation ends at one explicit source-qualified
boundary: recovery, second strand, IVT, later transformations, adapters and
PCR are not yet mapped for the v3 implementation. A shared default no longer
adds a terminal PCR beyond that boundary. The local cell/UMI/index roles retain
unassigned sequence lengths; eight acquired index cycles do not establish
eight-base index regions. This is an implementation/source-review limit, not
a claim that inherited protocol chemistry was unpublished.

Re-read one original: published PMC10409775 Methods lines265–293. The source
bundle scraper returned no sequence hits; two large extracted tables were
explicitly skipped by its size limit. This result supplies no absence evidence.
The source-defined v3/predecessor methods remain deferred.
