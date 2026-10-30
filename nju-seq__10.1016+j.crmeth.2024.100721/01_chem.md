# NJU-seq: selected HeLa MgR/no-MgR preparation

Tang et al., **An integrative platform for detection of RNA 2′-O-methylation reveals its broad distribution on mRNA**, DOI10.1016/j.crmeth.2024.100721.

## Sources and selected boundary

🟢 Published PMC10985248, Fig1 and Methods504–508 were inspected; one original, investigation closed. The source-bundle scraper was run (351candidates/10files; a16MB result table was skipped by its normal size limit). Scraper hits from uninspected data tables are not source evidence. Source analytical substrates, qPCR/RTL-P and Nm-VAQ chimera-probe assays are distinct from the production sequencing route. Supplementary FigS1 and TableS5 were identified but not read as production-oligo authorities. This is not a full linked-source or missing-information audit.

🟢 NJU-seq uses Mycoplasma genitalium RNaseR (MgR) to generate3′ endpoints enriched at Nm+1. Its important preparation boundary is ligation of an RNA tail before MgR digestion, followed by post-digestion5′ phosphorylation and NGS library preparation. The no-MgR NC is a matched preparation, not another protocol; a source-described pause does not establish perfect specificity, quantitative recovery or an exact molecular product for every RNA.

## Literal reagent and end limits

🟢 Methods508 prints the poly(A)-ssRNA adaptor as **`pAGCUAAAAAAAAAAAAp`**, added at3µl25µM per1µgRNA and ligated16°Covernight. The inventory retains the RNA letters, twelve-A tract and both terminal p marks as printed; it does not silently replace the reagent with an unblocked oligo(A) or an Illumina adapter.

The complete fragmentation-end state, tail-ligation catalyst/conditions, terminal-phosphate handling and MgR entry/product topology remain unassembled. Both printed tail phosphate marks are displayed in the partial end inventory, while the ligation geometry is explicitly unavailable. The source T4PNK is **3′-phosphatase minus** and is used after MgR/control workup for5′ phosphorylation; a preceding or following3′-end-healing reaction is not invented. No claim that an end-processing detail was unpublished is accepted.

## Order and physical fates

🟢 Selected HeLa route uses15µgRNA; the source's other samples use10µg. Methods504 describes H/M/R rRNA-probe hybridization95°C2min→22°C0.1°C/s/5min, RNaseH37°C30min and repetition “5/3 times.” The assignment of those two repeat counts and total fresh-stock additions is left unassembled, rather than guessed. DNaseI2U37°C30min follows, then VAHTS RNA Clean Beads. Cleavage products remain until physical workup; successful depletion is not asserted complete.

🟢 Methods508: fragment95°C10min →literal NaOAc/glycogen/ethanol workup →select1µg for each TR/NC preparation →ligate poly(A)-ssRNA tail →MgR30µg/60µl37°C30min or the same procedure without MgR →RNA Clean&Concentrator5 →5′PNK →historical NEBNext SmallRNA LibraryPrep Set for Illumina →6%PAGE →X-Ten sequencing. Untransferred fragment stock and actual workup/gel fractions remain explicit; no downstream count filter is presented as physical removal of RNA.

🟢 PNK mixture is39µlRNA,5µl10×buffer,5µl10mMATP and1µl3′-phosphatase-minusT4PNK. Source kit adapter, RT/PCR orderings and historical version remain unassembled; no current NEB constructor is imported to fill them. The kit is represented as contributing molecular stock, not merely a condition.

## Endpoint and open scope

🟢 The literal source phrase is “Fragments smaller than50nt” extracted from6%PAGE after the kit boundary. The model preserves that wording and leaves insert-versus-whole-library sizing unassigned; it does not silently correct the threshold or claim a sequencing-ready complete molecule below50nt. Source names Illumina X-Ten but not an assigned run/index mask here. Complete end geometry, kit internals and linked precursor/vendor versions remain unreviewed; no absence claim or public schematic is made.
