# DSB-Seq (2014)

[Published defining article](https://doi.org/10.3390/ijms150713111), *DNA Break Mapping Reveals Topoisomerase II Activity Genome-Wide* (2014).

🟢 One original: published PMC4139894 main3.4, common3.5, source3.1/3.2, Figure1b/5 and3.7. XML/PDF are representations of the same article. OriginalPDFp9 was visually checked for the literal200mL unit and wholePCR primer ordering. Mandatory bundle scraper run; mainExperimentalSection pointers verified. This is the article's distinct TdT/biotin/streptavidin/S1 DSB route, separately inspected rather than imported from the SSB nick-translation/DIG route or the later2018chapter.

## DSB-specific chemistry

🟢 Input: source500µg purified high molecular weight HCT116gDNA. Upstream common3.1/3.2 extraction remains outside the graph. TdT tailing uses3mlRochebuffer/24,000UTdT/0.5mMdCTP/0.02mMbiotin-16-dUTP,37°C30min. EDTA20mM stops reaction. A separate500µgDNA negative-label control substitutes TTP for biotin-16-dUTP; it is an actual independent preparation, not sequenced production DNA or a zero-background assumption.

🟢 Phenol/chloroform extraction, two2Mammoniumacetate/ethanol precipitations remove unincorporated label. The source then literally says dissolve in **200mL TE**. ActualoriginalPDFp9 confirms mL, not µL. We retain that printed value and leave intended volume unresolved; we do not silently substitute200µl or use an invented effective reaction concentration. This is a source-unit concern, not an OCR error or an established correction.

🟢 Sonicate with sourceBioruptor30×30smedium/30sice intervals to200–400bp. DynabeadskilobaseBINDER streptavidin captures biotinylated fragments. Wash4timesat50° and4timesRT in10mMTrispH7.5/1mMEDTA/2MNaCl. Disrupt complexes with10mMTris/1mMEDTA/1MNaCl/**2Mβ-mercaptoethanol75°C4h**; purify free DNA byQIAquick. S1nuclease30U in110µlrecommendedbuffer,37°C30min, removes biotinylated tails; QIAquick purification follows. Whole tail lengths, end distributions and incomplete/off-target products remain populations, not guaranteed exact break coordinates or quantitative tail removal.

## Locally declared common library preparation

🟢 Section3.5: EpicentreEND-Repair45minRT with33mMTrisacetate/66mMKacetate/10mMMgacetate/0.5mMDTT,0.25mMeachdNTP,1mMATP,1µlT4pol+PNKmix. Purify, then15UKlenowexo−/0.2mMdATP37°C30min creates a source3′A substrate end. Ligate0.1µlIlluminaadaptermix/1000UT4ligaseRT30min. Whole historical adapter orderings/ends remain unassembled; known substrate3′A is preserved. QIAquick;2%E-Gel200–400bpband/MinElute elute12µl.18PCRcycles98°C30s/65°C30s/72°C30s. Whole primers: Fw `ACACTCTTTCCCTACACGACGC`; Rv `CAAGCAGAAGACGGCATACGAGC` (source23nt ordering retained, not replaced by full canonicalP7).2.5%agarose220–500bpband/Qiagengel cleanup, GenomeAnalyzer cluster generation/sequencing. Source3.7 aligned36bp reads does not establish the complete physical run mask.

## Source limits and discrepancies

🟡 Methods3.1 states20µMetoposide37°C20min; Figure5caption states20µM5min for displayedSSB/DSBprofiles. Both source locations/values are retained; intended per-figure exposure remains unresolved. A future public DSB schematic must show this comparison and the literal200mL concern. The graph starts at preparedHMWgDNA rather than inventing a specific precursor preparation resolving the exposure. Negative-label control recovery/quantification fromFigure2 is distinct from productionNGS; exact downstream control processing beyond shared selection remains opaque. Source predecessor24 and historicalvendor originals are ordinary unread. No complete source/reaction or accepted absence audit.
