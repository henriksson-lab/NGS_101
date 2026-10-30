# RGEN-seq

Kuzin et al., [RGEN-seq for highly sensitive amplification-free screen of off-target sites of gene editors](https://doi.org/10.1038/s41598-021-03160-8), Scientific Reports (2021).

Three published originals inspected: definingPMC8654851 mainMethods104–110/Fig1actualPDFp2; SupplementaryProtocolMOESM2 steps0–H, original adapter tablep1; SupplementaryTableS1MOESM3 original workbook. Mandatory bundle scraper434hits/4files followed to these originals; no preprint recipes adopted. SI1figures and other workbooks remain ordinary unread targets, not confirmed irrelevant or missing.

🟢 Selected enzymatic-fragmentation route: optionalCIP/TdT/ddATP blocks pre-existing ends; then Cas9/sgRNA cuts create new sites. Actual full guide RNA is unassembled: the table's20-base targeting sequence is not a whole100-base sgRNA. Single-guide CD151 is a source-labeled example, not a new chemical variant. After RNase/proteinase treatment andSPRI purification, Klenow(exo−)/dNTP processing and p7L ligation precede a second purification. Select150–500ng in26µl for UltraII FS fragmentation/end-repair/A-tailing, then p5L ligation andSPRI cleanup. Actual elution transfers/stockremainder/QC/loadingallocations stay separate.

🟢 p7L comprises `rc_p7` and `p7_index`; p5L comprises `p5_index` and `rc_p5`. The original whole orders are retained in the module, with p7 donor5′phosphate/3′amino blocker, terminal phosphorothioate-protected3′T on rc_p7 and p5_index, and dephosphorylated p5L indicated in mainFig1. Source-known local3′T cohesive ends are typed; whole genomic occupancy/products remain unassembled. N6/N8 adapter windows are sample-index windows, not UMIs; actual study index values/read orientations remain unassigned.

🟢 Source treatment/mass/volume boundaries include37µl retained from39µl optional-block eluate;31µl retained from33µl after cleavage cleanup, then30.25µl for Klenow;30µl retained from32µl after first ligation, then26µl/150–500ng for fragmentation;23µl retained from25µl final eluate. The source selection recipe uses actual70µlSPRI aliquots/85%EtOH rather than an invented perfect-specificity cleanup. SI0 retains37µlTE eluate, while SIA6 separately specifies1–4µg DNA in40µl1xCutSmart. Recovered DNA mass/aliquot selection, the buffer-preparation interface and final-buffer convention remain unassembled. This is not a calculated3µl buffer addition or a new genomic DNA input.

🟢 ProductiveP5/P7 versus nonproductiveP5/P5 populations coexist after preparation. Nonproductive molecules are not declared removed by library cleanup. Flowcell clustering/washing is the later reporting/physical fate boundary described by the main paper. AligningR2/end-counting/off-target calling is computational, not a DNA purification reaction.

The main reportsR1=26/R2=151, whereas SupplementaryProtocolH saysR1=25/R2=100–300. Both locations/values are retained; the selected declaration uses the main reported26/151 run, without claiming that it is the exact SIHrun or resolving the difference. PCR-free applies to library construction; source qPCR remains an analytical QC operation.

Whole guide scaffolds, genomic end distributions, occupied ligation junctions, historical kit compositions and complete sequencing products remain unassembled. Adapter-only geometry is not proof of full ligation. This is a private positive partial graph, not a full reaction or absence audit.
