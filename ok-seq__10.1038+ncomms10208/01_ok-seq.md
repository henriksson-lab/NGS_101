# OK-seq — research notes

## Sources read

- 🟢 Petryk et al. 2016, DOI `10.1038/ncomms10208`, PMC full text and supplement.
- 🟡 DOI `10.1038/s41596-022-00793-5`, a later Nat. Protocols detailed protocol cited by the earlier reconstruction; not on disk and not inspected. Earlier page statements attributing the adapters and PCR primers to it were corrected in round 14: the defining Petryk Methods print them.

## Construct path

🟢 EdU-labelled genomic DNA is denatured and short single strands are recovered by sucrose-gradient size fractionation. EdU is click-biotinylated, RNA primers are hydrolysed and DNA 5′ ends are phosphorylated. 🟢 Two double-stranded adapters with six-base random overhangs ligate to the single-stranded fragment; biotin capture precedes PCR. The defining human paper explicitly publishes the four adapter oligos shown on the page.

🟢 The defining Petryk paper itself prints the complete forward PCR primer, the reverse primer carrying a six-base i7 index (printed as six asterisks between the P7 and TruSeq index arms), and the A1 Watson oligo used as the custom sequencing primer. Its forward primer begins `ATGATAC…`; the earlier model had silently prepended an extra `A`, which has been corrected. The reverse PCR primer completes the final base of the A2/TruSeq Index-1 binding site; the schematic therefore distinguishes the ligated pre-PCR adapter from the completed library rather than treating them as identical.

## Methods detail (round 2026-10-17-breadth-first-14)

🟢 PMC4729899 JATS Methods (`.xml.txt` lines 122–158; published PDF text lines 1658–1670 agree for the primers): 20 µM EdU for 2 min (EdC gave identical results); proteinase K/phenol extraction and ethanol precipitation; 94 °C 10 min denaturation; 5–30% sucrose gradient (SW28, 17 h, 26,000 g, 20 °C), <200-nt fractions concentrated on Amicon Ultra 10k; click with 2 mM CuSO4, 1 mM biotin-TEG-azide and 10 mM sodium ascorbate (45 min RT, dark), ethanol precipitation; RNA hydrolysis with 10 µL 2.5 M NaOH (30 min, 37 °C), acetic-acid neutralization and Bio-Gel P-30 spin column; 30 U T4 PNK (30 min 37 °C, 75 °C 10 min inactivation). Adapters annealed at 95 °C 5 s then 0.1 °C/s to 4 °C; strands annealed to 80 pM of each adapter (65 °C 10 min, 16 °C 5 min) and ligated at 16 °C overnight with 500 mM betaine, 5% PEG8000 and 4 Weiss U T4 DNA ligase. Capture on 200 µg MyOne Streptavidin T1 with stringent washes removing unligated adapters; **annealing and ligation were repeated on 20 µL bead suspension with 40 pM adapters**, then washed again. PCR: 12 cycles, Phusion, 10 µL bead suspension; MinElute; 150–400 bp excised from a 4% agarose gel. Custom Read 1 primer = A1watson; HiSeq 1000 single read 51.

🟡 The graph now represents each of these steps, the repeated on-bead ligation, MinElute and gel size selection with their discarded fractions, exact-state adapters and PCR primers (the reverse primer's six-base index as a placeholder), and Phusion/12 cycles.

🔴 Published-source note (shown on the page): the forward primer is described as "PE PCR Primer 1.0 truncated by 5 bases", consistent with removal of the 3′ GATCT, but the printed sequence also lacks the leading A of P5 (`ATGATAC…`, in both the JATS and the typeset PDF). The model uses the printed sequence; which reagent was intended remains unresolved, and the Illumina PE PCR Primer 1.0 original was not inspected in this pass. The acquired index-read length is not printed (only the 51-base single read), so the model no longer asserts a six-cycle index read.
