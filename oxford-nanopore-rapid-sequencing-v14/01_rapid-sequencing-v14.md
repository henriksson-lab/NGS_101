# Oxford Nanopore Rapid Sequencing V14 (SQK-RAD114)

## 1. What it is

SQK-RAD114 is Oxford Nanopore's PCR-free rapid genomic-DNA workflow: a fragmentation mix prepares DNA for rapid-adapter attachment immediately before loading.

## 2. Sources read

- 🟢 Oxford Nanopore Technologies, *Rapid Sequencing Kit V14 — gDNA (SQK-RAD114)*,
  RSE_9177_v114 Rev Q, 12 December 2025.
- 🟢 Oxford Nanopore Technologies, *Chemistry Technical Document*, rapid-based
  sequencing section.

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

🔴 The kit component sequences, attachment chemistry and exact motor-bearing strand are
proprietary. The schematic therefore shows role-labelled structures, not guessed bases.

## 5. Final library and read entry

The sequenceable product is a transposase-tagged genomic fragment bearing a motor-loaded
Rapid Adapter at the entering end. The adapter-bearing end enters the nanopore first; the
vendor sources do not support a more detailed base-level orientation.
