# Butt-seq

**Defining source:** doi:10.1101/gad.350434.123.

## Molecular boundary

🟢 TGIRT template switching captures the native 3′ end of purified nascent RNA without fragmentation; salt switching, second-strand synthesis and PCR complete the library.

## Evidence boundary

🟢 The article Materials and Methods (PMC10270195, extracted XML lines 177–195),
the supplementary detailed protocol, and the primer workbook all print oligos.
The scraper located these strings, but the original source text, not the scraper
output, is the evidence used below. The four terminal DNA variants are a pool,
not one invariant starter sequence.

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
panel. The page now shows these oligos and both source-backed blocked starter
duplexes, with a representative A-ended DNA from each four-member pool.
The post-PCR library remains role-level: its exact sequencing-primer sites
have not yet been established by the model.

## Published second-starter disagreement

🟢 The article Materials and Methods, extracted XML line 181, prints a
27-nt R1ME RNA (`rCrUrGrUrCrUrCrUrUrArUrArCrArCrArUrCrUrGrArCrGrCrUrGrC/3SpC3/`)
and 28-nt DNA (`GCAGCGTCAGATGTGTATAAGAGACAGN`). The primer workbook,
TGIRT Primers lines 11–15, instead prints 33-nt `R1MER_RNA` and 34-nt
`R1ME_A/T/C/G`, each with an additional 5′ `TCGTCG` on the DNA strand.
The detailed supplemental protocol says to use the primer-table sequences
and an equimolar mixture of its four DNA variants; the drawing therefore uses
the workbook version. 🟡 The article does not explain whether the shorter
Methods pair was an earlier variant or a textual truncation. The public page
displays both claims and the model choice.
