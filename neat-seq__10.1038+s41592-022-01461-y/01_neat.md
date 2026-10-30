# NEAT-seq

Defining original: **NEAT-seq: simultaneous profiling of intra-nuclear proteins, chromatin accessibility and gene expression in single cells**, DOI[10.1038/s41592-022-01461-y](https://doi.org/10.1038/s41592-022-01461-y).

Capped source/model pass,2026-10-29. Private and incomplete;structural graph validity is separate from source completeness. No missing-information claim has passed the two-reader gate.

## Originals inspected

- `s41592-022-01461-y_41592_2022_1461_MOESM1_ESM.pdf`: Published SI Protocolv1.0 pp4–6, AppendixIVp10: fixation/permeabilization/block/wash, fiveNPC hashing aliquots/pool/TF stain, Multiome transposition/continuation, preamp primer spike,22uL reporterPCR,0.6x/75uL transfer/1.2xSPRI, terminalddC block, reporter/index oligos and phosphorothioates. Rendered originalp10 verifies all displayed oligo notation.
- `s41592-022-01461-y_41592_2022_1461_MOESM3_ESM.xlsx`: Published TableS4: antibody barcode36nt members and separate HTO/antibody-concentration contexts. Table values not assigned to generic protocolv1.0 unknown hash members.
- `s41592-022-01461-y_PMC11192021.xml`: Explicit author manuscript (not finalVOR); Methods94–144: source antibody–streptavidin/biotin versus BDdirect conjugation, EcoSSB/no intermediatepurification, blocked DNA ddC, MultiomeRevA/35uL continuation, barnyardversusCD4 handles/run contexts. Used only separately qualified,not to override publishedSI22uL orHTONextera design.

## Chemistry and scope

🟢 The default local graph follows the published supplementary protocolv1.0. It fixes cells1.6%PFA2min,quenches withglycine5min,washes,permeabilizes with NP40/Tween/digitonin and performs actual pellet/wash separations. DNA blocking uses source random30mer material with a terminal dideoxycytidine. The local strand retains an actual terminalC nucleotide with blocked sugar/no3′OH,not a C3 spacer or backbone deletion. Whether the source30mer designation includes the finalddC is left unassembled; no31nt orN29C whole string is silently invented.

🟢 Antibody DNA is pre-associated with EcoSSB. The local state carries that qualitative protein association and physical source oligo;it does not calculate saturating occupancy,nuclease resistance or perfect retention. Source incubation/staining buffer is carried directly into staining without an invented intermediate purification. Initial free DNA block is washed before subsequent staining;it is not replenished where the source warns it may compete forSSB.

🟢 Five disjoint nuclei aliquots receive five distinct NPC hashing antibody preparations,then are washed and pooled. They are not five copies of every nucleus,nor five aliquots receiving one identical hash reagent. TF staining follows pooling,then washes,nuclei buffer and Flowmi filtration. The generic sourcev1.0 protocol is not assigned a specific accepted-manuscript cell sample or concentration experiment.

🟢 Source Multiome transposition precedes its historical partition/RT/recovery continuation. ADT/HTO preamp primers are added during step4.1. Separate manufacturer GEX/ATAC continuations remain coarse and source-qualified. ReporterPCR uses22µl preamp+25µlNEBNextHiFi+3µlprimers in the publishedSI,not35µl from the accepted manuscript. Source first0.6xSPRI keeps supernatant;75µl is transferred,untransferred material remains explicit. Added28µlSPRI yields source1.2xfinalselection;ethanol wash and20µlEB elution recover reporter libraries. No ideal class purity or universal molecule recovery is assumed.

🟢 Published ADT oligo `CAGACGTGTGCTCTTCCGATCT[12bp UMI][36bp Barcode]AAAAAAAAAAAAAAAAAAAAAAA*A*A`;HTO `GGCTCGGAGATGTGTATAAGAGACAG[12bp UMI][36bp Barcode]AAAAAAAAAAAAAAAAAAAAAAA*A*A`. The25A tail has the two printed terminal phosphorothioate bonds. Reporter UMI12,antibody identity36,HTO sample identity,library8bp indexes and kitcellbarcode are distinct meanings. N36 represents an unassigned source member,not random de novo antibody identity.

🟢 Published preamp primers `CAGACGTGTGCTCTTCCGATC` and `GGCTCGGAGATGTGTATAAGAGACAG`. Index primers: P5 `AATGATACGGCGACCACCGAGATCTACAC[8bp index]ACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T`;ADT P7 `CAAGCAGAAGACGGCATACGAGAT[8bp index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGAT*C*T`;HTO P7 `CAAGCAGAAGACGGCATACGAGAT[8bp index]GTCTCGTGGGCTCGGAGATGTGTATAAGAGAC*A*G`. These modifications belong to supplied oligos,not all amplified copies. Full PCR junctions,kit strand geometries and read-primer landings remain unassembled.

Version scope: the inspected main is an accepted author manuscript. It describes35µl preamp for reporterPCR and Truseq handles for barnyardHTOs,whereas publishedSIv1.0 uses22µl and generic NexteraHTOs. Those are separately qualified versions/contexts,not a confirmed publication error and not a merged recipe. Its barnyardADT MiSeq28/48/8/8 and NextSeq550 ATAC/RNA contexts are recorded as conditional source information;no run is assigned to the generic publishedSI graph without verifying applicability. Published main and historicalMultiomeRevA remain pending.

Component ancestry: publishedSI explicitly directs continuation of10xMultiome with local reporter modifications. Source EcoSSB–oligo association and intracellular blocking/staining are independently source-owned; no DOGMA/TEA/ASAP recipe is imported from endpoint similarity.

## Deferred originals and questions

- Full published main finalMethods and accepted-manuscript applicability; no whole-version equivalence assumed
- Historical 10xMultiome RevA original and source-specific whole capture/product/end/read geometry
- Actual experimental index members and exact randomblock/terminalddC count convention
- Remaining directly linked sources and independent source review

## Validation

Local source transcription and physical branch/state checks are recorded in the source-review JSON. Full chemical products,library completion and complete source coverage are not established.
