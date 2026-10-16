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
🟡 The post-PCR drawing uses `v2_Ad1.60` as one example index. Its P5, i5
and Nextera S5 sequence follow directly from that primer; the second
starter continues S5 into the mosaic end, then its variable one-base
overhang. The other side carries the first starter's variable overhang,
35-base read-2-side donor, eight-base UMI and P7. The actual RNA insert
length and bases vary by molecule. The primer-binding regions and their
orientation are computed from this construct.

🟢 The [author's analysis pipeline](https://github.com/albertdyu/BuTTSeq)
at commit `a2a815f04b2b7a0e0f65c23be5c4b24cffa82fc4` labels its three
expected FASTQs R1 = Read 1, R2 = UMI read, and R3 = Read 2. `Snakefile`
lines 302–311 extract an eight-base UMI from R2 with
`(?P<umi_1>.{8})` and attach it to both biological reads. The detailed
protocol, extracted DOCX lines 730–739, instructs changing the second
`RunInfo.xml` read from an eight-cycle indexed read to a non-indexed read
before demultiplexing. The R1/R3 lengths are not stated. 🟡 The constructed
duplex has sites for standard Nextera Read 1 and Index 2 plus TruSeq Index 1
and Read 2 primers. The index-1 primer reads the UMI; the index-2 primer
reads the example Ad1 i5 index in reverse-complement orientation on NextSeq
500. A canonical TruSeq Read 2 primer would read one fixed T from the 35-base
first-starter donor before the variable starter overhang and RNA-derived insert. These are
computed sequence-compatible primer sites; the paper and author pipeline
do not identify the physical primer oligos loaded on that sequencing run.

🟢 The cited Mohr 2013 TGIRT article (PMC3683930, extracted XML lines
126–134 and 192) describes an equimolar pool of starter duplexes with
different single-base 3′ DNA overhangs. This supports the capture mechanism,
while Butt-seq's own workbook establishes its particular blocked oligos.
Rhee and Burke 2004 describe TCEP as a stable reductant for RNA chemistry;
Butt-seq Methods line 181 independently specifies its annealing-buffer
concentrations. The 25 `v2_Ad1` names and sequences in the Butt-seq workbook
match the original Buenrostro 2015 supplementary primer table
(`nature14590` MOESM36) exactly. None of these linked sources identifies
the physical Butt-seq run-primer inventory. Illumina's
[NextSeq 500/550 custom-primer guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/nextseq/nextseq-500-550-custom-primers-guide-15057456-03.pdf)
(#15057456, v03) identifies separate Index 1 and Index 2 primer positions;
its [indexed-sequencing overview](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-08.pdf)
(#15057455) specifies
the reverse-complement i5 workflow for NextSeq systems. Those guides
support the run orientation, while Butt's oligos establish the library's
primer-site sequences.

## Read 2 first-base interpretation

🟢 The article's Results paragraph beginning “To analyze transcriptional
dynamics” says the first sequenced Read 2 base should reflect the last
nascent-RNA base. The published SCR2 DNA, however, has a 35-base donor
region ending in a fixed T before its variable one-base 3′ overhang; the
canonical [Illumina TruSeq Read 2 primer](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/OligoSeqsMultiplexOligo.htm)
is 34 bases. 🟡 With that standard
primer, the computed Read 2 path starts with fixed T, then the overhang,
then RNA-derived sequence. A 35-base custom primer ending in T could skip
the fixed base, but the article, supplement and pinned author pipeline do
not name a loaded Read 2 primer or provide a run `Primers.xml`. The public
page shows the standard-compatible geometry and this unresolved conflict;
it does not adopt the paper's first-base interpretation as proven.

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

## Published protocol authority

The published article Materials and Methods (PMC10270195, RNA-removal paragraph)
and its detailed supplemental protocol (lines 289–301) specify 3 µL each of
1 M / 1 N NaOH and HCl. The model uses these published values. The superseded
preprint is excluded from current protocol evidence; its revision differences
are not presented as a published-source discrepancy.
