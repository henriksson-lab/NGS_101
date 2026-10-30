# Cappable-seq

🟢 Ettwiller et al., 2016, [defining publication](https://doi.org/10.1186/s12864-016-2539-z). Published final PMC4782308, Methods150–190 and selected Results44–70; published Supplemental text Note A opening read. Two originals, capped positive review rather than complete source audit.

The front end distinguishes RNA 5′ end classes by adding a desthiobiotin-tagged cap to reactive di/triphosphorylated RNA with Vaccinia capping enzyme. This is introduced GTP-derived labeling, distinct from capturing an endogenous NAD cap. The source uses 3′-desthiobiotin-TEG-GTP; its 2′ positional alternative is a distinct tested reagent. SAM is present; the exact cap structure is not assembled here.

🟢 E.coli total RNA is isolated with FastRNA Blue, DNase I and Megaclear. A separate Bioanalyzer measurement reports RIN9. The selected preparation caps3µg RNA in50µl with50U VCE,0.5mM DTB-GTP and0.1mM SAM,37°C30min. Zymo Clean and Concentrator-5 retains≥200nt RNA, four washes,100µl eluate. Fragmentation adds2.5µl10×T4PNK buffer,94°C5min; AMPure/ethanol recovery returns100µl.

🟢 Seventy-five microlitres enters a15min 3′-phosphate removal reaction with8.2µl buffer and4µl **ATP-free T4PNK**, prepared by dialysis. The paragraph does not assign temperature; the local model does not default it. The unused25µl and subsequent unused PNK mixture are explicit unresolved physical remainders.

🟢 Fifty microlitres enters the first streptavidin capture. Beads are prepared with two400µl low-saltB washes followed by two400µl high-saltA washes. Capture lasts20min at room temperature; post-binding washes are twice200µlA then twice200µlB. Competitive biotin elution uses30µlB/1mM biotin for20min. AMPure/ethanol recovery gives60µl; **30µl enters a second capture**, with the other half retained outside the selected path. The second eluate is recovered into30µl. Physical washes/elution are not ideal biological class filters.

🟢 RppH15U/3.3µl Thermopol buffer,37°C60min removes the introduced cap to leave a5′ monophosphate;0.5µl0.5M EDTA and94°C2min terminate the reaction. AMPure recovery gives20µl. Historical NEBNext Small RNA library preparation and15PCR cycles precede E.coli MiSeqV3 single-end100nt acquisition. No current kit adapter is substituted. The source's non-enriched control omits streptavidin capture; its local route leaves intermediate allocation/workup unassembled rather than borrowing the enriched capture fractions.

🟡 The private graph records this order, physical side fractions and an independently scoped control. End structures, historical kit internals and primer/read geometry remain unassembled. The analysis trimming string is not treated as an ordering oligo. No missing-information conclusion or public schematic is claimed.
