# scNMT-seq

## 1. What it is

scNMT-seq jointly measures chromatin accessibility, endogenous DNA methylation and RNA expression by GpC marking followed by physical RNA/DNA separation and branch-specific libraries.

## 2. Authoritative source

🟢 Clark et al., *scNMT-seq enables joint profiling of chromatin accessibility, DNA methylation and transcription in single cells*, Nature Communications (2018), DOI `10.1038/s41467-018-03149-4`.

## 3. Protocol boundary

🟢 The RNA branch uses bead-capture Smart-seq2; the DNA branch uses M.CviPI marking followed by scBS-seq. Both branches belong on this protocol page because their pairing is the defining operation.

## 4. Construct path

🟢 Cells are treated with M.CviPI for 15 minutes before lysis. Oligo-dT magnetic beads capture RNA; wash supernatants are transferred to the DNA plate.

🟢 The RNA branch is Smart-seq2 cDNA followed by one-fifth-volume Nextera XT. The DNA branch is bisulfite conversion followed by the published scBS-seq two-round random-priming chemistry and PE1.0/iPCRTag amplification.

🟢 Both branches were sequenced paired-end. The schematic keeps their distinct Nextera and historical scBS-seq primer layouts separate.
