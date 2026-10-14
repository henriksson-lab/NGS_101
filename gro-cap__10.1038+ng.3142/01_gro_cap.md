# GRO-cap

**Defining paper:** Core et al., *Nature Genetics* (2014), doi:10.1038/ng.3142. The article, methods and supplements were read from PMCID PMC4254663.

## Molecular boundary

🟢 GRO-cap begins with BrUTP-labelled global nuclear run-on RNA. It then enriches capped nascent 5′ ends: uncapped 5′ phosphates are removed, tobacco acid pyrophosphatase exposes a ligatable monophosphate only at capped ends, and a no-TAP library controls background.

🟢 The paper states that GRO-cap libraries use TruSeq Small RNA adapters. The informative coordinate is the first transcribed nucleotide, not the run-on-labelled 3′ end.

## Evidence boundary

🟡 The source names the adapter kit but does not publish one immutable index set. The diagram retains the cap-selection chemistry and read direction while leaving kit-variable bases as roles.

## Oligo source audit

🟢 Checked the open author manuscript (`ng.3142_PMC4254663.pdf.txt`, lines 618–630) and publisher supplementary PDF (`ng.3142_41588_2014_BFng3142_MOESM77_ESM.pdf.txt`). The methods explicitly specify TruSeq Small RNA adapters and a cap-selected 5′ RNA end. The supplement supplies procedure detail but no printed adapter/index or sequencing-primer oligo table. The relevant 5′/3′ adapters, RT primer, PCR/index primers and read primer therefore remain kit-version dependent; no bases were promoted from a different TruSeq release.
