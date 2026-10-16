# Oxford Nanopore Rapid Sequencing V14 (SQK-RAD114)

## 1. What it is

SQK-RAD114 is Oxford Nanopore's PCR-free rapid genomic-DNA workflow: a fragmentation mix prepares DNA for rapid-adapter attachment immediately before loading.

## 2. Sources read

- 🟢 Oxford Nanopore Technologies, *Rapid Sequencing Kit V14 — gDNA (SQK-RAD114)*,
  RSE_9177_v114 Rev Q, 12 December 2025.
- 🟢 Oxford Nanopore Technologies, *Chemistry Technical Document*, rapid-based
  sequencing section and section 15 adapter sequences.

## 3. Molecular steps

1. 🟢 Fragmentation Mix tagmentates 100–150 ng high-molecular-weight genomic DNA for
   2 min at 30 °C, followed by 2 min at 80 °C.
2. 🟢 The transposase simultaneously cuts the DNA and installs transposase adapter
   sequences on the resulting ends.
3. 🟢 Diluted Rapid Adapter attaches to those tags for 5 min at room temperature in an
   enzyme-free reaction.
4. 🟢 A motor protein on the sequencing adapter controls passage of the attached strand
   through the pore. No synthesis sequencing primer is used.

## 4. Evidence boundary

🟢 The Chemistry Technical Document publishes the Rapid Adapter strand as
`CCTGTACTTCGTTCAGTTACGTATTGCT`.

The Fragmentation Mix tag, complete motor-bearing assembly, attachment chemistry and
exact motor-bearing strand are not assembled in the current model. These are open
molecular questions; the old audit row named desired component information rather
than an identified inaccessible document. This is not a confirmed publication-absence
claim.

The Rapid Adapter's source-printed strand is a local graph input. Its bases are not
placed across an unverified attachment junction in the final diagram. Identical
top-strand sequences printed for other ONT kits do not assign their bottom strands,
end geometry or motor assembly to SQK-RAD114.

## 5. Final library and read entry

The sequenceable product is a transposase-tagged genomic fragment bearing a motor-loaded
Rapid Adapter at the entering end. The adapter-bearing end enters the nanopore first; the
current diagram states that relationship at component level. The exact strand layout
and attachment junction remain unassembled; the known RA strand is shown separately.

## 6. Linked QC and readout sources

The linked Input DNA/RNA QC, Flow Cell Check, Flow Cell Wash Kit, MinKNOW and
Data Analysis originals have been read separately. Their roles remain distinct:
input measurement, available-pore checks, removal of a previous library,
acquisition settings and conversion of current signal into basecalls. The
wash know-how document explicitly describes a nuclease digesting the existing
DNA library and its sequencing adapters. That wash is outside the construction
of the incoming SQK-RAD114 library. Its barcoded demonstration does not assign
barcodes to this unbarcoded protocol.

The original [Dorado read-trimming documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_trimming/)
describes kit-aware removal of adapter or primer sequences from basecalled
reads. Its primer-based orientation and UMI features apply to specified primer
contexts; they do not provide a physical sequencing primer or UMI for
SQK-RAD114. The [basecalling overview](https://software-docs.nanoporetech.com/dorado/latest/basecaller/basecall_overview/)
also separates raw-signal processing from sequence processing. These software
snapshots describe processing capabilities; the graph does not choose a Dorado
version, trimming setting or analysis model.

The displayed RA strand therefore documents a reagent. Retained basecall ends
depend on processing. The [BAM](https://software-docs.nanoporetech.com/output-specifications/latest/read_formats/bam/)
and [POD5](https://software-docs.nanoporetech.com/output-specifications/latest/read_formats/pod5/)
specifications describe output records and raw-signal containers, respectively;
they are separate from the unresolved physical adapter assembly.
