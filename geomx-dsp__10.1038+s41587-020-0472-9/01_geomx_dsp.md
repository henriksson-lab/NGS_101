# GeoMx DSP NGS readout — source notes

## Scope

This page covers the named **GeoMx Digital Spatial Profiler (DSP) NGS readout**: a
photocleavable synthetic tag is released from an RNA probe or antibody inside a selected
area of illumination (AOI), collected into one well, PCR-completed with well-specific
GeoMx Seq Code primers and counted by short paired-end Illumina sequencing. It does not
describe the alternative nCounter readout.

## Sources read

- 🟢 Merritt et al., *Nature Biotechnology* (2020),
  [doi:10.1038/s41587-020-0472-9](https://doi.org/10.1038/s41587-020-0472-9), defines DSP
  with antibody or RNA probes attached to oligonucleotide tags through a photocleavable
  linker and spatially patterned UV release.
- 🟢 NanoString, *GeoMx DSP NGS Readout User Manual: Library Preparation & GeoMx NGS
  Pipeline*, MAN-10153-01 (April 2022), is the authoritative commercial library-prep
  source. Figure 2 gives every region and length; Tables 6–9 give PCR, cleanup and run
  settings.
- 🟢 The current Bruker GeoMx knowledge base independently recommends paired-end 2×27
  sequencing and reports an expected 162-bp Seq Code product.

## Published construct facts

🟢 Before PCR, the 5′→3′ photocleaved tag is `SPR1(20 nt) – UMI(14 nt) – RTS ID(12 nt) –
SPR2(20 nt)`. The 12-nt readout-tag sequence identifier names the biological target. The
14-nt random UMI is used to remove PCR duplicates.

🟢 GeoMx Seq Code PCR extends SPR1 from 20 to 33 nt and SPR2 from 20 to 34 nt, and adds
P5(29), i5(8), i7(8) and P7(24). Each well receives a unique dual-index pair, so the index
pair carries the AOI/collection-well identity. Eighteen amplification cycles follow a
30-min 37 °C UDG incubation and 10-min 50 °C deactivation. Products from the plate are
pooled and purified with two 1.2× AMPure XP cleanups.

🟢 The run is paired-end Read 1=27, Index 1=8, Index 2=8, Read 2=27. Read 1 therefore
contains UMI cycles 1–14 and RTS-ID cycles 15–26. Read 2 sees the reverse end: RTS-ID
cycles 1–12 and UMI cycles 13–26. The 27th base on each side enters the opposite PCR arm,
which permits the two reads to be stitched.

## What is not public

🔴 MAN-10153-01 prints region lengths and chemistry, but not the bases of SPR1, SPR2,
UMIs, target RTS-ID whitelist, Seq Code primers or UDI sets. The assay-specific RTS map
is supplied in a Probe Kit Configuration file, while the run-specific index values are
supplied in `SeqCodeIndices.csv`.

🟡 The final arm lengths and standard four-read Illumina geometry match canonical TruSeq
P5, Read 1, Read 2, P7 and index-primer regions. The schematic uses those canonical bases
only to compute and draw sequencing-primer binding, and marks every such segment inferred;
it does not claim that NanoString published the proprietary oligo sequences.
