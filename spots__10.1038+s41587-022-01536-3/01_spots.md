# SPOTS

Defining published original: **Integration of whole transcriptome spatial profiling with protein markers**, DOI[10.1038/s41587-022-01536-3](https://doi.org/10.1038/s41587-022-01536-3).

Private,incomplete capped source pass (2026-10-29);no missing-information claim accepted. Structural graph validity is separate from source completeness.

## Sources inspected

- `s41587-022-01536-3_41587_2022_1536_MOESM1_ESM.pdf`: Published SI Library Preparation Workflow pp14–21; antibody pool buffer exchange, staining,RT handoff,second strand/amplification,0.6xbranchselection,two1.9xADTcleanups,ADTindexPCR/final1xcleanup. Rendered originalp21 verifies protectedoligos andRPI_B1.
- `s41587-022-01536-3_41587_2022_1536_MOESM3_ESM.xlsx`: PublishedSupplementaryTable1:CD4RM4-5 identity AACAAGACCCTTGAG;panelclone/barcode context.
- `s41587-022-01536-3_PMC10272089.xml`: Accepted manuscript Methods86–120;source-version identity,SPOTS workflow and NextSeq500/NovaSeq6000 runs. Different manuscript title has sameDOI;final-version equivalence not established.

## Chemistry and scope

🟢 The published SI is the operational authority here. This route represents the optional buffer exchange recommended for panels over20 antibodies. After fixation,blocking,saponin permeabilization and90min antibody staining,the slide is washed and imaged; tissue-removal enzyme/SDS permeabilization and SSC washes precede the historical Visium reverse-transcription handoff. Whole capture probes,antibody-conjugated polymers and historical guide recipe are unassembled.

🟢 Published SI adds2µl100µM FB additive primer to second-strand synthesis (77µl total),then1µl0.2µM to cDNA amplification (102µl total). The accepted manuscript instead gives8.8µl and4µl. These are different source versions,not an established publication error;the model uses the published SI without merging values. qPCR cycle determination is an analytical allocation,not all cDNA consumed before amplification.

🟢 A0.6xSPRI selection retains both fractions: RNA-derived cDNA enriched on beads;ADT-enriched supernatant retained. The latter undergoes two1.9xSPRI cleanups (130µl then95µl reagent),ethanol washes and45µl elution;SI-PCR/RPI PCR6–10cycles;final1xSPRI cleanup/30µl elution. GEX continues through the historical Visium library procedure,whose detailed original remains unread.

🟢 Published oligos verbatim: dT25block `TTTTTTTTTTTTTTTTTTTTTTT*T*T*/3InvdT/`;FB additive `CCTTGGCACCCGAGAATT*C*C*A`;SI-PCR `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`;RPI_B1 `CAAGCAGAAGACGGCATACGAGAT AAGATTAC GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. Its table Index field is `GTAATCTT`,distinct from the printed primer orientation. No sequencer orientation inference is made. The actual blocker strand carries an explicit3′inverted-dT block;both oligos display their source modification notation. The inverted terminaldT and phosphorothioates belong to the physical oligos;PCR copies do not inherit their backbone modifications.

🟢 CD4 antibody barcode example is `AACAAGACCCTTGAG`;this inventory is not the full conjugated oligo. Spatial barcode,UMI,antibody identity and sample index have distinct meanings. Source reports GEX/ADT NextSeq500 (spleen) orNovaSeq6000 (breast cancer);exact run cycles/read-primer binding unassembled.

Component ancestry: the published SI explicitly continues Visium Spatial Gene Expression UserGuide with local modifications;matching endpoint does not establish all source recipes. Staining/blocking adapted from source reference17 in accepted Methods;predecessor recipe not read. No automatic CITE-seq/DBiT recipe inheritance.

## Deferred

- Historical Visium Spatial Gene Expression User Guide exact cited revision and complete downstream GEX recipe
- Final published main Methods and applicability of accepted-manuscript-only values
- Whole capture/ADT polymers,all identifier widths/placements and read-primer geometry
- Remaining SI/linked source review and independent absence-audit passes

## Validation

Protocol-local source transcription and branch/order checks;sourceguard and scoped molecular audit recorded in source review. No full base-resolved products or chemistry completion claimed.
