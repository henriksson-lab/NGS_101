# scCOOL-seq — research-note provenance

The shared principle, source reading, oligo comparison, and open evidence limits are
retained in the combined
[scNOMe-seq / scCOOL-seq research note](../scnome-seq__10.7554+eLife.23203/01_scnome-seq-sccool-seq.md).
This pointer avoids duplicating that evidence trail; the schematic is tied to Guo
*et al.* 2017.


## Bounded source pass — 2026-10-17

Read published PMC5539349 JATS Methods lines 143–153: single-cell nuclear
marking, the separate bulk variant, and single-cell library preparation and
sequencing. Read published Supplementary Data S1 `cr201782x22.pdf.txt`
identity and day-by-day protocol lines 1–123. The source manifest identifies
both originals through PMC Cloud. Ran the bundle oligo scraper and verified
the modified first nonamer primer and second nonamer primer in JATS153;
workbook motif hits were not treated as synthetic library oligos.

🟢 The single-cell route uses Klenow fragment (3′→5′ exo−), Exonuclease I
and AMPure recovery before M280 streptavidin capture, Klenow second-strand
synthesis, KAPA HiFi HotStart PCR on beads, and two final AMPure purifications
before QC, pooling and paired-end 150-base acquisition. The graph now names
these enzymes and represents primer digestion and the two recovery boundaries.
The exact PCR primer/index family remains unresolved, and the outer-arm
illustration retains inferred styling.

🟢 Published settings differ: the main article's library Methods (JATS153)
state 13 PCR cycles with Forward PE1.0 and reverse indexed primer; Data S1
p3 steps 1–2 (text114–120) specify NEB Universal/Index primers and 12 cycles.
The graph follows the main article's named-primer 13-cycle route; it does
not claim a corrected universal count or complete sequence equivalence of
the differently named primers. The schematic now displays both sources,
the choice and its uncertainty.

🟢 The main article's single-cell marking Methods (JATS145) specify heat
inactivation at 65 C for 20 min; Data S1 p1 step3 (text35) specifies 25 min.
The graph leaves that duration unassigned and displays the discrepancy.
Data S1 also describes six first-priming rounds and two alkaline bead washes;
these details do not establish that a condensed main-article route supplies
all identical conditions. This pass does not silently repair the supplement's
printed pM/nM wording in repeat-priming mixtures.

Deferred: individual rendered supplement figures and remaining tables;
scBS/PBAT predecessor originals and historical vendor primers/kit manuals;
source-specific lysis, lambda spike-in, SAM, protease and bisulfite-purification
states; repeated first-priming and full alkaline separation geometry; exact
PCR primer identities/index widths and modified intermediates. The graph's
shared marking front end still summarizes handling and does not establish
that scCOOL inherits the complete scNOMe chemistry. This is a capped pass,
not complete experimental validation or accepted absence of information.
