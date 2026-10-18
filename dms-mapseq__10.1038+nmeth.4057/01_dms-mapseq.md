# DMS-MaPseq — research notes

## Sources read

- 🟢 Zubradt et al. 2017, DOI [10.1038/nmeth.4057](https://doi.org/10.1038/nmeth.4057), PMC5508988, Methods and Supplementary Table 1 in the publisher's Supplementary Text and Tables PDF.

## Protocol boundary

The public page follows the genome-wide branch. The paper also describes targeted
gene-specific RT-PCR, a ten-base-UMI targeted branch, and Nextera XT conversion of
targeted amplicons. Those are not silently merged into this library.

## Molecular order

🟢 Cells are treated with 5% DMS, which modifies exposed A/C bases. Total RNA is
Zn2+-fragmented, size selected, dephosphorylated and ligated at its 3-prime end to preadenylated,
3-prime-blocked linker-2. The illustrated route then depletes rRNA before RT. TGIRT-III reverse transcription
at 57 °C reads through DMS adducts and encodes them as substitutions. RNA is hydrolyzed,
the 5-prime-phosphorylated cDNA is circularized with CircLigase, and 9–13 PCR cycles add
TruSeq-indexed Illumina arms. Custom primer oNTI202 produces a 50-nt single-end read.

## Source-transcribed oligos used by the model

- 🟢 linker-2: `5rApp/CACTCGGGCACCAAGGA/3ddC`
- 🟢 oCJ200-link2 RT primer: `5'/5phos/GATCGTCGGACTGTAGAACTCTGAACCTGTCG/iSp18/CAAGCAGAAGACGGCATACGAGATTCCTTGGTGCCCGAGTG`
- 🟢 library PCR primer oNTI231: `CAAGCAGAAGACGGCATACGA`
- 🟢 indexing PCR primer: `AATGATACGGCGACCACCGAGATCTACACGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNCGACAGGTTCAGAGTTC`
- 🟢 custom Read 1 primer oNTI202: `CGACAGGTTCAGAGTTCTACAGTCCGACGATC`

The internal `/iSp18/` is a non-nucleotide chemical spacer, not eighteen bases.
The current model keeps it in the physical RT primer, cDNA and circle as an
explicit role with no nucleotide length. Its PCR boundary uses that physical
single-stranded circle; it does not replace it with a nucleotide-only duplex. The six `N` positions in the indexing primer are the
source-declared sample index, not unknown adapter bases.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5508988 JATS99–105, genome-wide and targeted library preparation.

Added deadenylase/RecJ excess-linker digestion and RiboZero depletion before RT. Replaced invented RNA-removal enzyme with source NaOH chemical hydrolysis; its fragments stay unresolved until unrepresented recovery. Specified TGIRT57°C1.5h and9–13 PhusionHF cycles.

Deferred: Original SupplementaryTableS1 modified oligos and exact circle/PCR geometry, physical RNA/linker recovery, final gel selections; targeted/UMI branches remain excluded. This pass does not establish a complete source review or publication-wide absence.

## Local chemistry ownership migration — 2026-10-18

Removed the unused shared adapter-ligation recipe and its arbitrary 18-base
adapter representation. The source linker state now includes the printed
terminal `/3ddC/` as a C nucleotide beyond the 17-base linker core; that blocked
terminal base is outside the RT primer's binding core and is not automatically
copied into cDNA. The native ligation product uses the supplied linker parts and
an exact named RNA/linker junction, with 5′ preadenylation and blocked3′ end.

The source RT primer and PCR primers are actual molecular inputs. The physical
cDNA circle retains iSp18, while the native PCR geometry explicitly uses that
single-stranded circle. The former duplex view supplied an initial oNTI231 site
on a complementary strand that did not yet exist. Primer sequences and the
primer-derived final-library candidate remain source-positive; modelling the
first indexed-primer copy, oNTI231 binding to the new complement and polymerase
behavior at iSp18 remains an implementation/checker boundary. It is not missing
publication information. Source genomic RNA fragment selection is100–170 nt,
and linker/RNA hydrolysis fragments retain unresolved physical recovery.

Re-read two originals: published PMC5508988 Methods lines95–101 and the
published Supplementary TableS1 modified oligos (text lines140–210). The bundle
scraper was run and the source sequences/modifications inspected directly.
No source-absence audit or targeted-branch inheritance is claimed.
