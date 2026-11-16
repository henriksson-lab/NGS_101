# Splinkerette PCR

Potter CJ and Luo L. *Splinkerette PCR for Mapping Transposable Elements in
Drosophila*. PLoS ONE 5:e10168 (2010).
[Published article](https://doi.org/10.1371/journal.pone.0010168).

## Sources inspected

The published PDF and JATS original, Introduction (mechanistic comparison with
inverse PCR), Figure 1B–C, Table 1, and Discussion (unmodified adapter design);
published Protocol S1 (`pone.0010168.s002.doc`), Steps 2–7 and Tables 2–14.
The cached originals and their manifest live outside this checkout, under
`$CHEM_DATA/sources/splinkerette-pcr__10.1371+journal.pone.0010168/`.
The primer scraper was run across that bundle; the relevant sequences were
verified directly in Table 1, Protocol S1 Table 13, and rendered Figure 1C.
The insertion-location workbook is downloaded but is not the authority for this
synthetic molecular example. No missing-information/absence claim is made.

## Exact source oligos

🟢 All oligos below are 5′→3′ and unmodified. In particular, the adapter is
**unphosphorylated**, not an ordinary adapter with two ligatable strands.

```text
SPLNK-GATC-TOP GATCCCACTAGTGTCGACACCAGTCTCTAATTTTTTTTTTCAAAAAAA
SPLNK-BOT CGAAGAGTAACCGTTGCTAGGAGAGACCGTGGCTGAATGAGACTGGTGTCGACACTAGTGG
SPLNK#1 CGAAGAGTAACCGTTGCTAGGAGAGACC
SPLNK#2 GTGGCTGAATGAGACTGGTGTCGAC
3′SPLNK#1 CACTCAGACTCAATACGACAC
3′SPLNK#2 GGATGTCTCTTGCCGAC
3′SPLNK-SEQ CGGGACCACCTTATG
```

The implementation exports these sequences; the public page renders that data.
Other P-element ends and piggyBac use the different source primer sets in
Protocol S1 Table 13. Do not transplant this P-element 3′ primer set to an
arbitrary insertion construct.

## Structure and selective copying

🟢 Figure 1C shows the GATC overhang, core duplex, TOP hairpin and BOT arm.
🟡 The exact overlap is computed: TOP[4:27] pairs with BOT[38:61] (zero-based,
half-open coordinates). The hairpin representation pairs TOP[30:37] with
TOP[41:48], leaving TTTC as the loop. This is one explicit complementary fold
consistent with the figure; it is not a thermodynamic folding prediction.

🟢 BstYI is the source's recommended initial digest; BglII, BamHI and BfuCI
also supply compatible 5′ GATC ends. Restriction cleavage gives genomic 5′
phosphate/3′ hydroxyl ends. T4 DNA ligase seals BOT's recessed 3′ OH to genomic
5′ P. The opposite genomic 3′ OH/TOP 5′ OH junction cannot be sealed. Keeping
that nick is essential: TOP's hairpin must not be silently incorporated into
the eventual amplicon.

🟢 After denaturation, the transposon-specific T1 primer copies the ligated
genomic/BOT strand. S1 has the same polarity as BOT's arm; its complement is
created by that first copying event. It binds in the next cycle and copies back
into the genomic region. Two rounds of nested PCR use T1/S1 then T2/S2.

🟡 The simulator walks actual nucleotide predecessors, adds complementary
nucleotides from supplied dNTPs, preserves backbone identities and discovers
eligible 3′ ends. It never receives an expected amplicon sequence. The first
nested extensions still have a physical outer-template boundary on one side;
the subsequent cycle generates the bounded T2/S2 products. This is shown rather
than replacing the substrates with a desired endpoint.

🟢 The source says the TOP hairpin **reduces** end-repair priming/background.
It is not a chemically blocked 3′ end. The simulation does not assume universal
adapter-only immunity, infer kinetic competition, or predict a temperature-
dependent distribution of folds. Thermal separation opens the selected fold
before the demonstrated primer-directed copying.

## Experimental order and endpoint

The page transcribes the source digest, ligation, Phusion HF and nested-PCR
settings. Protocol S1 Table 6 lists two preliminary 98°C/64°C cycles without a
72°C extension entry; the page preserves that printed schedule. Table 7's
water is adjusted to the stated 25 µl total when using less than 1 µl template.

🟢 The amplified DNA is gel-purified or treated with Antarctic phosphatase and
Exonuclease I (37°C 2 h, then 80°C 15 min), and sequenced with T3 by Sanger
sequencing. This page's formal endpoint is the constructed flanking-sequence
amplicon. It does not invent an NGS library, sample index or sequencing run plan.

## Example scope

🟡 The short locus combines the exact T1, T2 and T3 source sequences in the
source-supported nested order, with **invented** intervening DNA, genomic flank
and restriction-site distance. Those synthetic values are not P-element
reference coordinates or a reconstruction of a mapped fly line. All inherited
synthetic template bases remain inferred through extension.

The molecular simulation retains the other restriction fragment, the original
TOP hairpin oligo, parental templates and first-cycle PCR products. Presentation
movement out of the viewport is not cleanup, hydrolysis or selective destruction.
Only representative successful PCR complexes are expanded; reagent concentrations,
whole-vessel population growth, mispriming, yields and end accessibility under
competing conformations are not predicted.

## Method ancestry

🟢 Protocol S1 explicitly says the four adapter/adapter-primer oligos derive
from Horn et al. (2007), *Splinkerette PCR for more efficient characterization
of gene trap events*, Nature Genetics 39:933–934. The Introduction also credits
Devon et al. (1995) for the method. These are method/component ancestry links,
not physical imports of prepared material and not evidence that their complete
protocols are interchangeable. The current chemistry uses the directly inspected
2010 published article and its own supplement.
