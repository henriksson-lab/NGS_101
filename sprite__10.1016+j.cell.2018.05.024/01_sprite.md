# SPRITE — research notes

## Sources read

- 🟢 Quinodoz et al. 2018, DOI [10.1016/j.cell.2018.05.024](https://doi.org/10.1016/j.cell.2018.05.024), PMC6548320, main text, Methods and tag-design section.

## Molecular order

🟢 DSG/formaldehyde-crosslinked nuclear complexes are fragmented, sparsely coupled to
NHS beads, washed, blunt-ended, phosphorylated and dA-tailed. Beads undergo five
split–ligate–quench–pool rounds: DPM, Odd, Even, Odd and Terminal. Crosslinks are then
reversed and the DPM/Terminal arms support Illumina library PCR. Reads sharing the entire
tag chain are assigned to one multiway SPRITE cluster.

🟢 Bounded pass 2026-10-17 (round breadth-first-19; PMC6548320 JATS "DNA Repair",
"Split-pool ligation" and "Sequencing library generation", lines 953–1060): end
processing is three bead-bound steps — T4 PNK (NEB) 37 °C 1 h; NEBNext End Repair
cocktail (T4 DNA polymerase + T4 PNK) 20 °C 1 h; Klenow 3′→5′ exo− (NEBNext dA-tailing
Module) 37 °C 1 h — each quenched with Modified RLT and buffer-exchanged. Each tag round
uses 2.4 µL of 45 µM tag with 2X Instant Sticky End Ligation Master Mix (NEB) at 20 °C
for 1 h, Modified RLT inactivation, pooling and three 45 °C SPRITE Wash Buffer washes to
remove free tags. Aliquots are digested with Proteinase K (50 °C 1 h), reversed overnight
at 65 °C, purified on Zymo DNA Clean and Concentrator columns (7x binding buffer),
amplified with Q5 Hot-Start Mastermix and cleaned twice with 0.7X SPRI (AMPure XP). The
graph now declares these catalysts, cleanups and their waste fates; the final PCR cycle
count stays unassigned (Methods gives 9–12 cycles only for the molarity-estimation aliquot).

## Source-owned endpoint and ordering oligos

🟢 Original Table S5 has been recovered from the [author archive](https://authors.library.caltech.edu/records/xrfk6-20h25/files/mmc5.xlsx?download=1).
Its 96-member DPM, Odd, Even and Terminal Top/Bottom families and indexed
PCR primers are the original2018 table, rather than substituted successor
plates. The schematic displays representative first-well oligos with their
phosphates and DPM3′ C3 spacer intact.

🟢 The DPM barcode column contains eight variable bases. The original DPM
schematic (`nihms972930u2.jpg`) highlights nine bases including a fixed A:
first-well `AAACACCCA`. Odd/Even barcode columns have sixteen variable bases;
complete ordering strands also contain invariant ligation sequences.
Terminal Bottom barcode strings are staggered and vary in length. These
column lengths must not be treated as complete physical tag or junction lengths.
The five barcode parts jointly identify a spatial complex.

🟢 The original final-DNA schematic (`nihms972930u4.jpg`) shows DPM at both
ends of genomic DNA and tag extensions outside DPM. The endpoint retains
an opposite-DPM role instead of depicting a generic insert containing only
five contiguous barcode windows. The full base-level junction chain is not
reconstructed in this model; bracketed roles do not define exact lengths.

🟢 “Final Library Amplification Primers” in the article prints the forward
primer modeled as `{{= sprite_protocol.PCR_FORWARD }}` and reverse primer
`{{= sprite_protocol.PCR_REVERSE }}`. They supply standard Read1/Read2
binding handles and one eight-base i7 index. There is no i5 in this source
primer pair. These are library PCR primers, distinct from the shorter loaded
sequencing primers. The page uses this article example without assigning
one index to every experimental library.

## Source discrepancies and read-layout correction

The DPM design paragraph says its barcode is read only in Read 2. The same
article's “Barcode identification” says the first DNA tag is taken from the
beginning of Read 1, and “Alignment to genome” specifies `--trim5 11` for DPM.
The page follows those explicit processing sections and flags the earlier
Read 2 label as an apparent inconsistency. The eleven-base processing trim is
not treated as the variable barcode-column length.

The article reverse PCR primer includes `GCCTAGCC`. Original Table S5 row
`2P_657` has Barcode `GCCTAGCC` but its primer includes `GGCTAGGC`, the reverse
complement. The source declarations use different index orientations; the page
retains the article primer literally and does not equate it to that table
member or silently normalize one sequence to the other.
The article does not identify its example as `2P_657`, so this comparison does
not establish an erroneous primer or the index used in every experiment.

🟢 Original article Methods and the accepted author PDF specify paired reads
of at least 115×100 bases. The former 100/75 model values were repository
assumptions and have been removed. The new model records source minima without
asserting an exact acquired run length or computing cycle spans across the
unassembled junctions.

## Scope and remaining source review

The schematic models the original DNA branch. The article also describes an
RNA/DNA variant using an RPM hybrid RNA/DNA tag, reverse transcription and a
cDNA-specific adapter. These remain a distinct unassembled branch here;
its sequencing-primer geometry is not inferred from the DNA endpoint.
The later detailed protocol DOI10.1038/s41596-021-00633-y is a successor,
not evidence that every later well member was used in2018.

The author-deposited accepted PDF and Tables S1–S6 are now retrievable. Tables
S2–S4 give chromosome/contact/hub analyses; Table S6 gives read-filtering
counts, FISH probe regions and contact/microscopy comparisons. These are not
additional library-adapter tables. Original publisher supplementary figures, movie samples and historical lab
guides now have located independent checks. The expanded linked-source
inventory still needs the pipeline/wiki and procedural-video scope recorded. No publication
absence is accepted.

## Recovered original figures and historical lab protocols

🟢 The [official author-hosted publisher PDF](https://guttmanlab.caltech.edu/wp-content/uploads/2021/08/PIIS0092867418306366.pdf)
contains all seven original Supplementary Figures S1–S7 on PDF pages28–39.
S1 includes tag-completion and ligation-QC results; S4 includes the RNA/DNA
branch cartoon and strand-specificity result. Other panels concern genomic
contact analyses, microscopy and nuclear organization. The publisher's original
Movie S1 is a rotating fluorescence-microscopy reconstruction, separate from
the lab's procedural video. These recovered originals supersede the earlier
blanket statement that the original figures/movie were inaccessible.

🟢 The [current official author protocols page](https://guttmanlab.caltech.edu/lab-resources/protocols/)
links the DNA guide whose title says last edited August9,2018. An ENCODE
deposit supplies a distinct April5,2018 guide. Both print the molecular
duplex examples, explain strand-selective initial PCR, call the article
reverse-primer example `2P_barcoded_85`, and print separate DPM-ligation-QC
primers. This reinforces that the article example should not be identified
as TableS5 member `2P_657`.

Located guide wording needs qualification: §3.1 calls the target overhang
5′A; §3.4 calls the full primers around30nt although the printed sequences
are58/66nt. Its explanatory five further PCR cycles also differ from seven
in the actual §6.2 program. The page retains literal full primer sequences
and does not assign these guide overhang polarities or cycle counts. These observations were independently confirmed against the original guides
and the August PDF pages were visually checked. The affected public page
shows them explicitly; they do not justify inventing end geometry or a PCR cycle count.

The QC primers now enter a separate aliquot branch in the graph and terminate
as QC material, rather than being assigned to final library amplification.
Their precursor duplex remains unassembled; source-positive strands do not by
themselves prove the PCR geometry. Both historical guides are checked by the
local source-transcription suite.

The author procedural video and article-linked pipeline/wiki still require
their recovered source/version scope to be independently recorded.
No absence claim is accepted.

## Located analysis-source scope

The prepublication pipeline snapshot `fcf7662a6b2f3bba32440fcfd2914665857016bc`
uses DPM in Read1 and Y–Odd–Even–Odd in Read2, with spacer matching between
those Read2 windows. All96 members of each family were independently compared
with original TableS5. DPM recognition strings are reverse complements of its
eight-base Barcode cells; Y strings equal the TerminalBottom Barcode cells;
Odd/Even recognition strings omit the first base of their sixteen-base bottom
Barcode cells. These are parser windows, not new ordering strands or a claim
that all physical tag lengths equal their matching windows.

All seven current wiki pages were read independently. They describe barcode
identification, alignment, filtering, clustering and contact analysis. Their later
software configuration is kept separate from the historical experiment.

The successor publisher procedure video was recovered and visually sampled by
two readers. It is distinct from the original microscopy movie. Direct retrieval
of the author-linked YouTube version failed independently; equivalence to the
publisher file and complete narration inspection remain open. Neither sampled
frames nor unavailable video access establish an absence of sequence information.
