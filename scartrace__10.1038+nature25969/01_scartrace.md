# ScarTrace

Defining published original: **Whole-organism clone tracing using single-cell sequencing**, DOI[10.1038/nature25969](https://doi.org/10.1038/nature25969).

Private,incomplete capped source pass (2026-10-29);no missing-information claim accepted. Structural graph validity is separate from source completeness.

## Sources inspected

- `vanOudenaarden_Key_2018_Alemany_Whole-organism-clone-tracing-using-single-cell-sequencing.pdf`: Published institutional mirror identity/title/DOI;Methods1385–1402 andFig1d:RT/secondstrand beforeproteinaseK/nestedGFP PCR,wellpool/oilseparation/split,SORT-seq andadapterPCR,PE75NextSeq. Fig1d labelsIVT transcriptomebranch.
- `nature25969_41586_2018_BFnature25969_MOESM6_ESM.xlsx`: PublishedSupplementaryTable1 firstCELSeq2-BC-001 row andtablecolumnidentity;fullRToligo,8Nsegment,8ntcellbarcode,25T+VNanchor.
- `nature25969_41586_2018_BFnature25969_MOESM7_ESM.xlsx`: PublishedSupplementaryTable2 first384sheetrow;sourceT7-eGFPprimercomponent,GTCACCGGwellbarcode andAGGACGACGGCAACTACAAGGFParm. Wholeprimerjunctionunassembled,notinventedfromadjacentcolumns.

## Chemistry and scope

🟢 Published Methods: DAPI-negative live single cells are sorted into384wellplates containing mineral oil,uniquely barcoded RT primers,dNTPs,spike-ins and RNase inhibitor. Plates are spun/storedminus80°C;cells lyse65°C5min. RT42°C and second-strand synthesis16°C precede proteinaseK access to genomicDNA and nestedGFP PCR. The secondPCR introduces separate scarwellbarcodes. Original gDNA-only branch omits step1 (RT/secondstrand).

🟢 All wells are physically pooled,aqueous material is separated from oil,and the recovered mixture is allocated to scar and transcriptome libraries. This is an aliquot split,not DNA/RNA purification. Transcriptomes follow SORT-seq;Fig1d explicitly labelsIVT on pooledcDNA. Exact linked SORT-seq operational continuation remains unreviewed. Scar adapterPCR introduces IlluminaTruSeqadapters;source reports paired-end75bp NextSeq. No new cleanups invented between these stated operations.

🟢 PublishedTable1 first RTprimer is `GCCGGTAATACGACTCACTATAGGCCTTGGCACCCGAGAATTCCANNNNNNNNCGTCTAATTTTTTTTTTTTTTTTTTTTTTTTTTVN`;the explicit Barcode column is `CGTCTAAT`. Main transcriptomeanalysis identifiesRead1cellbarcode/UMI andRead2biological information. The8Nsegment is assigned UMI semantics with that source context;the25T+VNanchor is not another identifier. Table2 first scarbarcode `GTCACCGG` is separate from RNAbarcode and inherited GFPscar;its fixed GFP arm is `AGGACGACGGCAACTACAAG`. The adjacent Table2 pieces are one source primer reagent with an explicitly incomplete component inventory. They are not separately added oligos,and are not silently joined into a validated wholeprimer.

🟢 Histone-GFP recorder scars can occur in multiple transgene copies per cell. The local lineage feature describes inherited allele information,not a unique currentcell identifier,randommolecule UMI or a predicted indel string. Recorder creation remains upstream of actual NGS recovery. Computational clone filtering does not physically remove molecules in the graph.

Ancestry: published Methods expressly directs the transcriptome branch to SORT-seq (reference18) and gives a detailed Single-cellScarTrace ProtocolExchange original (reference33,10.1038/protex.2018.017). This is branch-specific operational inheritance,not identity with LINNAEUS/scGESTALT. Both directly linked recipes remain deferred under the three-original cap.

## Deferred

- Direct detailed Single-cell ScarTrace original https://doi.org/10.1038/protex.2018.017
- Cited SORT-seq source reference18 and exact transcriptome continuation
- Full nestedPCR/adapterPCR oligos,Table2 component junctions and product/read geometry
- Remaining supplements and independent source review

## Validation

Protocol-local source transcription and branch/order checks;sourceguard and scoped molecular audit recorded in source review. No full base-resolved products or chemistry completion claimed.
