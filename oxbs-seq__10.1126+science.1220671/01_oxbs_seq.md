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

The final construct is a single-index TruSeq bisulfite library. Paired or single-end genomic reads carry the converted-base signal; the index read identifies the library.

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
