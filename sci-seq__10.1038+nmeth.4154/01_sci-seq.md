# SCI-seq

## 1. Protocol boundary

SCI-seq is the Vitak et al. single-cell whole-genome method. LAND or crosslink/SDS treatment
depletes nucleosomes while retaining nuclei. Indexed Tn5 complexes give each nucleus a first
bipartite index; nuclei are pooled and resorted; indexed PCR supplies a second bipartite index.
This is not sci-ATAC-seq: nucleosome depletion is intended to make coverage genome-wide.

## 2. Sources read

🟢 Vitak et al., *Sequencing thousands of single-cell genomes with combinatorial indexing*,
*Nature Methods* (2017), DOI `10.1038/nmeth.4154`; PMC author manuscript `PMC5908213`.

🟢 Amini et al., *Haplotype-resolved whole-genome sequencing by contiguity-preserving
transposition and combinatorial indexing*, *Nature Genetics* (2014), DOI `10.1038/ng.3119`,
especially corrected Supplementary Table 4. Vitak explicitly says that its transposome complexes,
associated sequences and custom sequencing primers are those of Amini et al.

🟢 Vitak et al. Supplementary Protocol, Tables 1–3, prints the complete
transposome oligo set, indexed PCR primers and all four custom sequencing primers
for this implementation. The archived files are
`nmeth.4154_41592_2017_BFnmeth4154_MOESM252_ESM.pdf` and its companion
index-layout archive.

## 3. Indexed transposomes

The transferred P5-side oligo has this structure:

`S5 – TCCACGC – [8-base transposase i5] – GCGATCGAGGACGGC – ME`

The transferred P7-side oligo has this structure:

`S7 – CTGTCCCTGTCC – [8-base transposase i7] – CACCGTCTCCGCCTC – ME`

🟢 The constant regions and complete index oligos are printed in Vitak Supplementary
Protocol Table 1. The universal non-transferred strand is
`/5Phos/CTGTCTCTTATACACATCT`. Eight P5-side and twelve P7-side oligos make 96 transposome
combinations.

## 4. Workflow

🟢 Two thousand nuclei are sorted into each of 96 tagmentation wells, receive one indexed
transposome combination and incubate at 55 °C for 15 min.

🟢 Wells are pooled, and 22 nuclei are sorted into each PCR well. Each well already contains two
10-base indexed primers. SDS dissociates Tn5; the xSDS workflow also reverses crosslinks. A 72 °C
extension precedes 15–20 monitored PCR cycles.

🟢 A cell identity is the four-part tuple: 8-base transposase i7, 10-base PCR i7, 8-base
transposase i5 and 10-base PCR i5. One 96-well PCR plate permits 9,216 combinations.

## 5. Sequencing primers and cycle structure

🟢 Read 1: `GCGATCGAGGACGGCAGATGTGTATAAGAGACAG`.

🟢 Read 2: `CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`.

🟢 Index 1: `CTGTCTCTTATACACATCTGAGGCGGAGACGGTG`.

🟢 Index 2: `CTGTCTCTTATACACATCTGCCGTCCTCGATCGC`. Vitak Supplementary Protocol
Table 3 prints this fourth custom primer; it also equals the reverse-complement
C15/ME site derived from the P5-side construct.

🟢 The custom NextSeq recipe is R1: 50 imaged cycles; I1: 8 imaged, 27 dark, 10 imaged; I2:
8 imaged, 21 dark, 10 imaged; R2: 50 imaged. Thus each index read crosses a transposase barcode,
skips the constant adapter/connector, and then reads the PCR-well barcode.

## Bounded curation round 3 (2026-10-16)

Inspected: Published PMC5908213 JATS122–138.

Added55°C/15min tagmentation,22 nuclei per PCR well and15–20-cycle real-time PCR labels. Candidate insertion product recovery is unresolved. Source distinguishes QRP/DOP control preparations from sci-seq.

Deferred: Published indexed transposome predecessor/oligos, LAND versus xSDS handling and reverse-crosslinking graph branches; exact imaged/dark-cycle mask validation. This capped pass does not establish a complete source review or publication-wide absence.
