# PLATO: ribosome display and recovered-RNA sequencing

Zhu et al. (2013), **Protein interaction discovery using parallel analysis of translated ORFs (PLATO)**. DOI: [10.1038/nbt.2539](https://doi.org/10.1038/nbt.2539).

This is a capped, incomplete source pass and private operational model. The defining published article and published supplement control the recipe. Associated PMC4110636 is an author manuscript: selected Methods were encountered during retrieval, but are not used to fill published-source gaps. Full equivalence was not established.

## Sources inspected

🟢 Published first-party PDF `source_pdf_a5d474dd04edc863d82f.pdf.txt`: printed title/DOI and main approach, especially lines43–65. 🟢 Published `nbt.2539_41587_2013_BFnbt2539_MOESM14_ESM.pdf.txt`: Online Methods, ribosome display/precipitation and mRNA preparation/Illumina sequencing/analysis (100–251). The mandatory scraper found69 candidates in the published SI; pointers were checked in those Methods. No primer body is adopted from a collapsed hit.

The published MOESM15 workbook containing Table4 is available but uninspected after the three-original cap. Its sequences are pending reading, not absent or retrieval-blocked. Source URLs and cache identity remain in the external manifest; see [retrieval](ref/MANIFEST.md).

## Molecular and operational scope

🟢 Stop-codon-free full-length ORFs are displayed through noncovalent mRNA–ribosome–protein association. This differs from covalent puromycin display. A source display population retains RNA as its encoding polymer when EDTA disrupts the ribosome; newly synthesized cDNA is a separate molecular class.

1. Source Gateway LR constructs and T7B/TolAk PCR produce linear ribosome-display templates. Purified18super-pools are combined equally by DNA weight, transcribed with T7Ribomax, and RNA purified with RNeasy. Detailed recovered bacterial/backbone geometry is unassembled.
2. RTS100E.coliHY translation is30°C15min. An actual12.5µL aliquot is diluted in cold RD selection buffer; the remaining translation preparation has an explicit unresolved allocation. After14000g5min4°C, the supernatant is transferred and sediment is removed.
3. Selected prepared/blocked bait beads bind the displayed population at4°C6h. Six500µL washes physically remove wash fractions. No perfect binding or retained-purity prediction is made. The model keeps prepared bait manufacture coarse instead of assigning one source application's chemistry to all baits.
4. EB20 containing20mMEDTA disrupts complexes at37°C10min. Release itself does not discard beads, proteins or RNA. Source RNeasy purification follows.
5. Recovered RNA is fragmented with NEBNext magnesium at94°C3min, Spin-50 cleaned and copied with SuperScriptIII/TolA100RT. ExonucleaseI37°C30min acts before RNaseH37°C30min; digestion remnants stay until QIAquick purification. These are source operations, not a computation of nuclease specificity or yield.
6. TdT polyadenylates **cDNA3′ ends**, followed by QIAquick purification. Adaptor-(dT)24 receives one extension cycle before P5-PRDREV and the additional PCR mixture are added for30cycles. This order is distinct from adding both primers initially.
7. Purified firstPCR material supplies100ng to a second10cycle P5-PRDREV/INDEX PCR. The secondPCR is purified, a separate analytical qPCR aliquot quantified with P5/P7_2, and production libraries combined equally. Analytical material does not become the sequencing substrate.

🟢 Published acquisition is HiSeq2000 multiplex single-end50 cycles with custom PRDREV-attB2-SP. Seven-nucleotide barcodes distinguish source samples; they are not UMIs. ORF assignment and normalization are computational analyses, not physical DNA changes. Full primer landing/library strand geometry is unassembled.

## Open reading and ancestry

Published Table4, directly cited RD backbone methods and historical kit details are pending. The publication explicitly builds a Gateway-compatible vector from a ribosome-display backbone (SI OnlineMethods); this is a component-level inheritance claim, not evidence that an existing PhIP-Seq or modern ribosome-display recipe applies. No unverified catalogue ancestry edge is proposed. No missing-information or whole-protocol completion claim is made.
