# CUT&RUN.ChIP — bounded source notes

Defining DOI: https://doi.org/10.1016/j.molcel.2018.10.046

Associated author manuscript PMC6475595 main107–129; one original. Published-final equivalence not established.

🟢 Native yeast nuclei factory:500mlYPD30°C/OD~.7, ~5e9cells; Zymolyase/Ficoll, ten frozen aliquots~5e8cells each. The selected graph begins with one thawed aliquot. Formaldehyde belongs to the first referenceDNA preparation, not the experimental yeast nuclei.

🟢 ~200ulConA captureRT≥5min, DigWash+2mMEDTA wash, mouseantiFLAG2h–overnight4°C. AssociatedAM extraction literally says1:2 00 final dilution; ratio is left unassigned, not silently normalized or called a publication error. RemoveEDTA thereafter; rabbitantiMouseIgG1:1001h4°C; pA-MN source360ugml,1:400 descriptor1h4°C with intervening washes.1.2mlDigWash/0°C then2mMfinalCaCl2 onice5min.

🟢 Add1.2ml2xSTOP(150mMNaCl/20mMEDTA/4mMEGTA/50ugmlRNaseA/40ugmlglycogen) including5–50pgml mostlymononucleosome DNA from formaldehyde-crosslinked MNase Drosophila chromatin. Descriptor concentrations are retained without silently imposing stock/final conventions.37°C20min;FLAGpeptide13ugml;5min13000rpm4°C; collect released supernatant on magnet. Nuclear/bead residue is excluded physically.

🟢 Split supernatant into7–8ChIP aliquots and save input4°C. Source lists IgG/H2B/H2A.Z/H4ac/H3K4me3/H3K36me3 without a complete mapping of7–8allocations; graph models selectedH3K4me3 and input, preserves other real aliquots unresolved. HistoneAb1:100overnight4°C;20ulProteinA for rabbitpolyclonals orProteinG for mouseH3K36me3,4°C15min; wash once with Wash(noDigitonin)+.05%Tween20. This is a second **physical immunoprecipitation**, not ordinary CUT&RUN or computational histone assignment.

🟢 ChIP extractionbuffer150mMNaCl/10mMEDTA/2mMEGTA/.1%SDS/.2mgmlPK/40ugmlglycogen plus second5–50pgml mostlymononucleosomeS.pombeDNA. Input receives matched SDS/PK/secondreference concentrations and bypasses serialIP.50°C20min;bufferedPCI thenphase-lockCHCl3; fresh2ul2mgmlglycogen(4ug) incollectiontube;2–2.5volEtOHovernight−20°C/13000rpm10min4°C/100%EtOHwash/airdry/resuspend25ul.1xTE. Drosophila and S.pombe can exchange first/second roles; graph selects source-described order without treating spike genomes as identifier oligos.

🟢 Historical KAPAHyperprep/TruSeq adapters **without size-selection**, AMPureXP removal of unligated adapters;25cyclepaired-endHiSeq2500. Cited Skene2017/2018 and actual linked KAPAguide are ordinary unread originals; no generic kit recipe has been imported.25cycles are acquisition cycles, not PCR cycles.

🟡 Private graph preserves the targetMNase-released supernatant→actual histoneChIP→DNArecovery boundary, savedinput, two referenceDNA additions and carrier glycogen stages. No perfect capture or referenceDNA chromatin occupancy is implied. Complete native ends, ligation geometry, whole adapters/primers/NT inventory and run indexmask remain unassembled.

Source investigation closed within the selected three-original ceiling before model construction. This capped positive read is not a complete source/reaction/absence audit. Missing-information claims are empty.

Deferred originals and applicability questions:

- Cited Skene2017/2018 detailed CUTRUN original methods and source-linked KAPA HyperPrep guide unread within bounded scope; no inherited whole kit recipe.
