# oxBS-seq — oxidative bisulfite sequencing

## Sources read

- 🟢 Booth *et al.*, “Quantitative sequencing of 5-methylcytosine and 5-hydroxymethylcytosine at single-base resolution,” *Science* (2012), DOI `10.1126/science.1220671` (defining paper).
- 🟢 Booth *et al.*, “Oxidative bisulfite sequencing of 5-methylcytosine and 5-hydroxymethylcytosine,” *Nature Protocols* (2013), DOI `10.1038/nprot.2013.115` (open author manuscript, PMCID: PMC3919000).

## Molecular interpretation

- 🟢 Genomic DNA is sonicated, supplied with conversion controls and prepared with the Illumina TruSeq DNA kit before conversion.
- 🟢 Potassium perruthenate selectively oxidizes 5hmC to 5fC. Bisulfite then converts unmodified C and 5fC-derived material to U, while 5mC remains protected.
- 🟢 Uracil-tolerant PCR copies original C and 5hmC as T; original 5mC reads as C.
- 🟢 A matched ordinary BS-seq aliquot reads 5mC + 5hmC as C. 5hmC abundance is calculated from the BS-minus-oxBS difference; it is not directly distinguished within one oxBS molecule.

## Final readout

The separately labelled 2013 reference endpoint has native TruSeq end roles whose
selected historical strands and index architecture are unassigned. Step 35 names MiSeq
or HiSeq, without assigning an acquired read layout in this bounded reconstruction.
No index length or acquisition follows from compatibility alone. The final construct
illustrates genomic PCR signal; control-derived products remain separately bounded.

## Implementation and recovered control oligos

The schematic's sonication/TruSeq preparation follows the authors' later 2013
detailed protocol. The defining 2012 Science article demonstrates oxRRBS. The
original Science supplementary library details are still inaccessible; the
later implementation must not establish the original study's exact kit/index
selection.

🟢 The recovered 2013 Supplementary Table 2 (original Springer MOESM16, p. 1)
prints seven sequencing-control DNAs plus a digestion control, and the seven
primers Q3, Q6, Q7, Q8, DFP1, DRP1 and post-bisulfite DFP2. These synthesize
or assay control molecules; they are not the selected genomic-library adapter
primers. The table defines `5` as 5-methylcytosine. Its caption places 5mC
through modified primers and 5hmC through 5hmCTP incorporation during PCR.
The caption's literal statement that other “Cs and Gs” are 5hmC is ambiguous
for a displayed single strand; no guanine-conversion rule is inferred from it.
The control modifications need strand-aware reconstruction before a molecular
control model is added.

🟢 **Published caption discrepancy:** 2013 Supplementary Fig. 1b (original
Springer MOESM11, p. 1) describes DNA made by PCR with d5fCTP, then names
5hmC as the substrate in one conversion clause, while the opening and final
clauses describe 5fC. The model follows the main protocol and surrounding
caption, representing 5fC → U; the repeated 5hmC in that clause is an apparent
caption error; its intended wording remains unresolved.

The remaining recovered originals show Sanger conversion controls (Fig. 1),
HPLC confirmation of oxidation (Fig. 2), oxidant solution colors (Fig. 3),
TaqI digestion-control traces (Fig. 4), and a statistical coverage example
(Table 1). They do not replace Supplementary Table 2's control oligo identities
or establish the selected genomic-library PCR pair.

## Bounded curation round 2 (2026-10-16)

Inspected: Original detailed2013 protocol PMC3919000 JATS43,55,265–319, Steps8–19.

Added alkaline denaturation before ice-cooled perruthenate oxidation and corrected post-bisulfite state from T to U before PCR. Public notice records contamination-sensitive preparation and need for parallel nonoxidized BS control. Preserved2012oxRRBS versus2013TruSeq distinction.

Deferred: Defining2012original/SOM; full2013supplements and library primers; parallel BS branch, buffer-exchange/recovery graph stages. This capped pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch03)

Inspected three originals: defining Science 2012 `22539555.pdf.txt` pp.934–937,
oxidation discussion and oxRRBS experiment (text 110–145, 470–488); published detailed
2013 PMC3919000 JATS Development 55, Steps 6–22 (255–378) and 32–36 (440–481);
2013 Supplementary Fig. 1b original MOESM11 p.1 caption. Mandatory bundle scraper
returned zero hits; this does not close image-based tables or establish absence.
The original Science SOM `booth-som.revision.1.pdf` remains a named deferred target.

The native graph explicitly selects the later sonicated 2013 reference implementation.
It does not stand in for the original 2012 oxRRBS construction. Step 8 delegates native
TruSeq preparation to a historical kit: this is an explicit unresolved boundary rather
than invented end-repair/dA/Y geometry. Sequencing controls (Step 7) remain separate DNA
molecules, bead cleanup reserves an ordinary BS aliquot (Step 9), digestion control is
added to the oxidation aliquot (Step 10), and water-equilibrated P-6 eluate recovery is
shown before and after alkaline denaturation/oxidation. Both denatured strand populations
feed oxidation. The digestion control persists through the oxidized/converted states;
control-derived amplification/QC fates remain separately unresolved, not discarded.

🟢 Both aliquots receive modified EpiTect FFPE bisulfite conversion (Steps 20–22).
For oxBS, original C and 5fC become U before PCR while 5mC is retained. For ordinary BS,
5hmC forms the stable methyl-sulfonate adduct described in Development 55; it does not
remain physically unchanged or become U. PCR then reports that adduct as C. The model
uses actual converted-strand scenes for native TruSeq primer-cocktail/PfuTurbo Cx PCR
(2 µl input in 50 µl, 12–15 cycles), followed by AMPure recovery (Steps 32–34).
No canonical genomic adapter/PCR primer or index width is inferred from kit compatibility.
The published Supplementary Fig.1b caption discrepancy remains visible publicly.

Deferred: original Science SOM and RRBS construction; exact historical TruSeq strands,
indices, native primer sites and acquired read layout; strand-specific control reconstruction
and Step 23–31 digestion QC. The prior recovered control-table history is preserved,
but that table was not re-read in this capped pass. No absence claim is accepted.
