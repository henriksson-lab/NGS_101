# Ribo-seq — ribosome profiling

## 1. What it is

Ribo-seq converts nuclease-protected ribosome footprints into a sequencing library whose mapped positions report translating ribosomes.

## 2. Sources and evidence

- 🟢 Ingolia NT et al. *Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling.* Science (2009). [doi:10.1126/science.1168978](https://doi.org/10.1126/science.1168978).
- 🟢 Ingolia NT et al. *The ribosome profiling strategy for monitoring translation in vivo by deep sequencing of ribosome-protected mRNA fragments.* Nature Protocols (2012). [doi:10.1038/nprot.2012.086](https://doi.org/10.1038/nprot.2012.086).
- 🟢 The detailed protocol specifies RNase digestion, monosome recovery, footprint size selection, 3′ adapter ligation, reverse transcription, cDNA circularization and indexed PCR.
- 🟡 The public schematic expresses historical experiment-specific index oligos as canonical TruSeq run roles so primer geometry can be computed without claiming an unsupported sample index.

## 3. Molecular path

The defining selection happens before library construction: RNA outside the ribosome is digested. Circularization of first-strand cDNA supplies the second PCR boundary without a second RNA-adapter ligation.
