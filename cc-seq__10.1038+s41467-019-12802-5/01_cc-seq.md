# CC-seq: covalent cleavage-complex sequencing

Gittens et al. (2019), *A nucleotide resolution map of Top2-linked DNA breaks in the yeast and human genome*. [Defining article](https://doi.org/10.1038/s41467-019-12802-5).

## Selected source scope

🟢 Published PMC6813358 main Methods205–225 and calibrated-library paragraph261; Results53–63; published Supplementary Information Table5 actualPDFp17 and SupplementaryFig1 caption321–327. Two originals inspected in this bounded pass; mandatory oligo scraper run. Source identities, access and deferred work are in the source review. The source-data workbook exceeded the fetch limit; main identifies it as experimental source data, not the adapter table. Neither peer-review nor reporting-summary PDFs were read as chemistry authority.

## Preparation

🟢 Nonproteolysed DNA is isolated by ethanol fixation, SDS/phenol-chloroform extraction and ethanol precipitation, with RNaseA treatment. For the selected human component, source cells are distributed over three400µl aliquots; subsequently aliquots are combined to1ml. Physical remainder/allocation beyond the described amounts is unassigned. Yeast additionally uses Zymolyase spheroplasting; this is a separate source branch, not inferred nuclease digestion.

🟢 Sonicate300–400bp, dilute withTE, add high-salt Triton/sarcosyl buffer and distribute approximately20µg per silica column. Reapply flowthrough; sixTEN washes then100µlTES elution. This antigen-independent protein-affinity boundary differs from antibody capture. Eluates are pooled to500µl, proteolysed with1mg/mlProteinaseK30min60°C, precipitated with glycogen/NaOAc/ethanol, and resuspended52µl. Actual2µlQubit aliquot is separate from50µl production.

🟢 Historical NEBNextUltraII end repair/dA-tailing and customP7 ligation precede AMPure cleanup. TDP2 is then added to the50µl eluate plus50µl reaction buffer:3µl10µM recombinanthumanTDP2,1h37°C. After cleanup/elution52µl, customP5 is ligated WITHOUT another end-repair step. The source explains this avoids dephosphorylating3′phosphate termini from3′CCs. Whole historical kit composition/aliquot bridge and sealed products remain unassembled.

🟢 Final17µl eluate splits2µlQubit/15µlPCR. Source universalP5/indexedP7 primers and historical kit PCR have unassembled whole-orderings/cycles in the selected scope. AMPure cleanup,200–600bpBluePippin and Bioanalyzer precede multiplex pooling. SourceMiSeqKitv3–150cycles PE75 or NextSeq500Kitv2–75cycles PE42 are separate alternatives, not silently rounded. Computational aligner/terminalMapper and filtering are not physical library selection.

## Adapter structures

🟢 Table5 gives annealed P5/P7 and defines Ti as inverted dT joined through3′–3′phosphodiester. The opposite5′terminalA(P5) orG(P7) is a nucleoside with3′OH, as source states. The source structure, including all fixed bases and terminalmodifications, must remain distinct from an ordinary phosphorylated duplex. Model-owned constructors compute complementary placements; whole reaction products are unassembled.

## Scope and questions

🟡 This is a private partial molecular model, not full validation or a complete source audit. Main identifies the source mechanism, but exact historical kit internals and whole PCR/index/read landing remain unassembled. Calibrated yeast experiments add separately prepared human DNA after sonication at1:100 molar ratio; not imported into the selected human preparation. Spo11/Top2 are source branches of this single defining method, not additional collection entries. No confirmed publication-absence claim.

Computed modified-ordering inventory (actual5′→3′; nonstandardterminalnucleosidesdescribedabove):

- P5top: {{= chem.SOURCE_ORDERINGS["p5top"] }}
- P5bottom: {{= chem.SOURCE_ORDERINGS["p5bottom"] }}
- P7top: {{= chem.SOURCE_ORDERINGS["p7top"] }}
- P7bottom: {{= chem.SOURCE_ORDERINGS["p7bottom"] }}

Bothsourceadapters have a3′Tcohesiveligationend andoppositeblockedend. The invertedterminalTisretained asanucleotidewithnonstandard3′–3′linkage, not discardedtofit apair. SourceTable5pairs32bpP5/33bpP7; wholeproductsunassembled.
