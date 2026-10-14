# Butt-seq

**Defining source:** doi:10.1101/gad.350434.123.

## Molecular boundary

🟢 TGIRT template switching captures the native 3′ end of purified nascent RNA without fragmentation; salt switching, second-strand synthesis and PCR complete the library.

## Evidence boundary

🟢 The supplementary detailed protocol and primer workbook are the oligo authority; the model does not invent a single invariant starter-duplex sequence.

## Primer inventory from the source

🟢 The article's `Supplemental_Primers.xlsx` (PMC10270195) separates
**First Strand Synthesis**, **Second Strand Synthesis** and **PCR Primers**.
The first-strand `SCR2_DNA_A/T/C/G` oligos share
5′-CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT-3′
followed by A, T, C or G respectively. The partner `SCR2R_RNA` is the
RNA oligo
5′-rArArGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrUrCrArC/3SpC3/-3′.
The terminal `/3SpC3/` blocker is explicit in the workbook.

🟢 The second-strand family `R1ME_A/T/C/G` is
5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG-3′ followed by A, T, C or G.
Its RNA partner `R1MER_RNA` is
5′-rCrUrGrUrCrUrCrUrUrArUrArCrArCrArUrCrUrGrArCrGrCrUrGrCrCrGrArCrGrA/3SpC3/-3′.
The detailed protocol (lines 36–69 of the extracted DOCX) directs an
equimolar hand mix of all four terminal variants for **each** DNA family
before annealing with its respective RNA oligo. Consequently there is no
single fixed terminal base for either duplex.

🟢 The PCR sheet lists universal `SCP7`
(5′-CAAGCAGAAGACGGCATACGAG-3′) and the indexed `v2_Ad1` primer
panel. The current model's role-only first strand, second strand and PCR
states omit these source-published sequence families and RNA-end blockers;
its unresolved sequencing-primer placement is a modeling task, not an
unpublished starter-oligo sequence.
