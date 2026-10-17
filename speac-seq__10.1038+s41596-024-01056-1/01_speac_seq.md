# SPEAC-seq

**Defining source:** doi:10.1038/s41596-024-01056-1 — Faust Akl et al., *Nature
Protocols* (2024), read as the PMC author manuscript PMC11805652 (JATS XML; local
`.txt` twin). Sections inspected: Materials → Oligonucleotides (extracted lines
153–175), droplet breakage and gRNA amplification (lines 636–689, including the PCR
mix table at 649–665 and the cycling table at 671–683), library preparation,
sequencing and analysis (lines 691–705), data availability (line 753). Supporting
original: Wheeler et al., *Science* 2023 (doi:10.1126/science.abq4822; PMC10066924).

## Molecular boundary

🟢 One CRISPR-perturbed cell and one fluorescent reporter cell are co-cultured per droplet; responsive droplets are sorted and their integrated guide cassettes PCR-amplified.

🟢 Droplets are frozen at −80 °C, thawed, broken with 1 mL of 20% (vol/vol) PFO in
HFE-7500 and spun at 1,000*g* for 30 s; the aqueous layer goes into the Qiagen Blood
and Tissue kit (two 200 µL elutions), is concentrated on AMPure XP beads at 2×
beads-to-sample into 40 µL, and is used whole as the template of one 100 µL ExTaq PCR
(95 °C 1 min; 28–35 × 95 °C 30 s / 53 °C 30 s / 72 °C 30 s; 72 °C 10 min). Samples
from fewer than 2,000 sorted droplets were amplified at 35 cycles. The product is
purified again on AMPure XP at 2×.

🟢 There is no ligation and no second indexing PCR: that single PCR installs both
Illumina ends onto the integrated lentiCRISPRv2 cassette.

## Oligos (verbatim, Oligonucleotides list, lines 153–175)

🟢 P5 staggered primer mix for Illumina sequencing (IDT), nine species, pooled,
100 µM stock, 0.5 µM final:

- 0 nt stagger — 5′-AT GATACGGCGACCACCGAGAT CTACACT CTTT CCCTACACGACGCT CTT CCGAT CT GAT GT CCACGAGGT CT CT-3′ (spacing as printed)
- 1 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGATGTCCACGAGGTCTCT-3′
- 2 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGCGATGTCCACGAGGTCTCT-3′
- 3 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTAGCGATGTCCACGAGGTCTCT-3′
- 4 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCAACGATGTCCACGAGGTCTCT-3′
- 5 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTGCACCGATGTCCACGAGGTCTCT-3′
- 6 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACGCAACGATGTCCACGAGGTCTCT-3′
- 7 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAGACCCGATGTCCACGAGGTCTCT-3′
- 8 nt stagger — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTGAAGACCCTTGTGGAAAGGACGAAACACCG-3′

🟢 P7 index primer for lentiCRISPRv2 (unique per reaction, 0.5 µM final):
5′-CAAGCAGAAGACGGCATACGAGAT**NNNNNNNN**GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT-3′

🟡 Derived: every species carries P5 + the TruSeq Read 1 sequence
(…CTTCCGATCT) followed by a stagger window and a vector-annealing 3′ end, so the
mix offsets the gRNA position between clusters on this low-diversity amplicon; it is
a phasing device and carries no identifier. The eight N of the P7 primer are the only
sample index in the library, which the model keeps as an unassigned sample-index
feature.

## Published discrepancies (shown on the schematic page)

🟢/🟡 Stagger labels versus inserted bases. Relative to the shared …CTTCCGATCT
end the printed species insert 0, 0, G, AG, CAA, TGCAC, ACGCAA and GAAGACCC bases
for the labels 0–7; the "0 nt" row is printed with internal spaces and lacks the
first A of P5, and the "8 nt" row ends in TTGTGGAAAGGACGAAACACCG rather than the
CGATGTCCACGAGGTCTCT 3′ end shared by the other eight. All nine are transcribed as
printed; only the space-joined "0 nt" reading is marked inferred. Whether these are
author-manuscript typesetting losses is unresolved — the publisher's typeset version
was not inspected and nothing is repaired here.

🟢 PCR table arithmetic. The P7 row prints a "5 mM" stock with 10 µL per 100 µL for
a 0.5 µM final (line 659), which would require a 5 µM stock; the P5 row (line 657)
is internally consistent. The model keeps the stated 0.5 µM working concentration
and leaves the stock unassigned.

## Evidence boundary

🟢 This is a pooled functional screen. Sequencing measures guide enrichment by sorted fraction, not a cell-resolved transcriptome.

🟢 Libraries are sized on a Bioanalyzer high-sensitivity chip (the lentiCRISPRv2
Brie library averages 285 bp), quantified with a KAPA library quantification kit,
pooled at 4 nM each and sequenced to at least 100 reads per gRNA (sorted droplet).
Reads are trimmed with Cutadapt between the anchors ACACCG … GTTTTAG and counted
with MAGeCK.

🔴 The protocol names no sequencer, read length or cycle plan, so the page asserts no
read geometry. 🔴 The lentiCRISPRv2 sequence between the two primer sites is not
printed, so the amplicon interior and the final library stay role-level and the PCR
geometry is declared unavailable rather than computed.

## Open questions

- The Addgene lentiCRISPRv2 (#73632) / Brie library map would locate both primer
  sites and give the exact amplicon; it was not retrieved in this pass.
- The publisher's typeset *Nature Protocols* version (and PMC
  NIHMS2026815-supplement Supp_Data_1/2, both behind a PMC download gate) could
  settle the stagger-label and "5 mM" questions.
- GSE200457 raw reads could support an inference of the invariant amplicon flanks,
  but cannot establish ordered oligo identity.
