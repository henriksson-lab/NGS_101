# ssGap-seq

[Published defining article](https://doi.org/10.1093/nar/gkad263), *RecA and SSB genome-wide distribution in ssDNA gaps and ends in Escherichia coli* (2023).

🟢 Two originals: publishedPMC10287960 Methods59–85 and publishedSI FigS1workflow. Mandatory scraper returns genomicmotifs, not adopted as library primers. Mutant/protein/dose branches remain one defining method; no automatic extra protocols.

## Preparation and chemical meaning

🟢 Formaldehyde1%10minRT then125mMglycine10min, twiceTBSwash/store−80°C. Selected non-UV branch uses a crosslinked pellet from40mlOD~.5,2mlChIPlysisbuffer andCovarisS2 8×60s. Main says200–300bp whereas SI FigS1 says300–500bp; both are retained, with main selected for this graph. Clarify16000g5min4°C; select250µllysate plus250µlbuffer and5µlantibody overnight4°C, retaining remaining clarifiedstock.

🟢 Chicken antiRecA orantiSSB uses100µlPrecipHen2hRT; rabbitantiRecA has a distinctProteinGDynabead branch. Selected graph uses chicken antiRecA, without treating these binding proteins as DNA covalent tags or inheriting another ChIP recipe. Wash2×eachlysis/highsalt/LiCl/TE buffers. Elute twice250µl0.1MNaHCO3/1%SDS5min, pool actual500µl eluates. Amicon10kTE×3 thenBioSpinP6 bufferexchange gives~40–50µl eluted chromatin.

🟢 Add460µlfresh5Mbisulfite/20mMhydroquinone and incubate37°C18h under non-denaturing conditions BEFORE reversecrosslinking. Fourwaterwashes,0.3MNaOH20minRTdesulphonation, fourTEwashes andBioSpinP6. Source claims substantiallygreater conversion ofssDNA than dsDNA; do not assert perfect binary conversion. ChemicalC→U remains distinct from Treportedbysequencing and opposite-strandG→A computationalcalls.

🟢 Add1µl10µg/µlRNaseA37°C10min; finalreversecrosslinkreaction200µl withProteinaseK.25mg/ml/1%SDS65°C4h. SI FigS1 instead labelsRNaseA.2mg/ml; maindose per200µl implies.05mg/ml ifallvolumeisassignedthere, so preserve separate sourceamounts. Do not silently average or treat the supplement as automatically overriding main. MonarchPCR/DNAcleanup elutes20µlEB; select10µl forxGenMethyl-SeqAdaptase library with16PCRcycles, leaving remainingIPDNAstock. No secondbisulfite step is inferred from the kit name. TwoSPRI rounds, actualQubit/qPCRQC aliquots, MiniSeqPE150MidOutput300cycles. Computationaldedup/BSmapping/windowfilters are not physical molecular discard.

## Limits

🟡 Whole adapter/orderings, historical kit formulation, full molecular products/read landing and cited predecessor protocols remain unassembled. This capped pass is not publicationabsence or full molecular validation.
