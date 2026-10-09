# PacBio Kinnex full-length RNA

## 1. What it is

Kinnex full-length RNA converts full-length cDNA molecules into ordered multi-insert arrays, then encloses each array in a SMRTbell for PacBio HiFi sequencing and computational segmentation.

## 2. Source read

🟢 Pacific Biosciences, *Preparing Kinnex libraries using the Kinnex full-length RNA
kit*, document 103-238-700 Rev09, August 2026.

## 3. Construct-changing operations

1. 🟢 Poly(A)-primed, template-switch cDNA synthesis is followed by Iso-Seq barcode PCR.
2. 🟢 Eight parallel Kinnex PCRs add orientation-specific segmentation sequences.
3. 🟢 Equal volumes of the eight products are pooled.
4. 🟢 Kinnex enzyme, ligase and a barcoded Kinnex terminal adapter assemble the cDNA
   segments into a linear array at 45 °C for 60 min; DNA repair follows for 30 min.
5. 🟢 Nuclease treatment removes incomplete molecules with exposed ends. The recovered
   product is the SMRTbell library used for sequencing-primer annealing and polymerase
   binding.

## 4. Evidence boundary

🔴 Kinnex oligo sequences are proprietary. Array order, segmentation-adapter roles, nuclease selection and SMRTbell topology can still be represented without guessed bases.

## 5. Read layout

The PacBio sequencing primer binds a proprietary terminal-adapter site. Polymerase circles
the closed array repeatedly to produce HiFi sequence; ordered segmentation junctions split
that consensus molecule back into its eight full-length cDNA segments.
