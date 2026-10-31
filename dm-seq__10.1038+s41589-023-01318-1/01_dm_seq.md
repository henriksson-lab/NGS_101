# DM-Seq

## Sources read and boundaries

🟢 `s41589-023-01318-1_PMC10763687.xml.txt`: Associated accepted manuscriptNIHMS1948388: Results91–109, optimizedMethods199 versus unoptimized185, and GBMacquisition227.

🟢 `s41589-023-01318-1_41589_2023_1318_MOESM1_ESM.pdf.txt`: PublishedSI identity/DOI, Table1 originalp2 rendered (/tmp/dm16-table.png), detailedfinalprotocol572–800 and copy-strand discussion630–635.

Cappedpass: at most10minutes/3originals; definingmain/chemistrysupplements before linkedoriginals. Ordinary deferredhistoricalgeometry is not an accessbarrier.

## Selected chemistry and strand identities

🟢 DM-Seq requires distinguishing original unmodifiedC, original5mC, original5hmC and a newly synthesized methylated copy. OriginalDNA is ligated to5-propynylcytosine-containing (5pyC) adapters. Fully methylated copy primer(v2) and Bstlargefragment create a copy using dmCTP. This is synthetic methylation, not a second native methylome. Original unmodified CpGs are protected by MBP-M.MpeI-N374K/CxSAM as5-carboxymethylcytosine(5cxmC), and original5hmC is protected as glucosyl5hmC. A3A deaminates original5mC into physicalT and unprotectedC into physicalU; protected5cxmC/5gmC/5pyC remain substitutedcytosines. The local class preserves these adduct identities and distinguishes the actual chemical product from copied/reference bases.

The source reports the copy strand does not readily sequence; the model does not physically discard it at deamination or assume perfect suppression of all copy-derived reads. Possible original complement/displacedDNA, copied primers, enzymes and other products remain until a source physical separation. The insert-only complementary relation in local tests is not a complete primer-derived Bst product or a guaranteed reaction register.

## Source operation order and branches

🟢 ExtractedDNA receives the source1:100control mixture (unmethylatedpUC19/CpG-methylatedlambda/T4hmCDNA), then Covaris/SPRI and separate Qubit aliquots. UltraIIendprep20°C30min/65°C30min precedes annealed5pyCYadapter ligation20°C15min and1.2xSPRI. The source mass allocated to the copy reaction is10ng in main or~1–10ng in SI; an unallocated remainder is not silently consumed. Bst65°C30min follows copy-primer annealing. βGT/UDPglucose37°C1h precedes MungBean30°C30min and1.2xSPRI. Digestion products remain until cleanup. MBP-M.MpeI-N374K/CxSAM37°Covernight/95°C5min is followed by ProteinaseK37°C15min and1.2xSPRI.

The **optimizedGBM branch** uses source main50µlA3A after snap-cooling, with detailed A3A conditions directly citing an ACEprotocol not independently followed here. Its actual GBM libraries were NextSeqHighOutputv2.5single-end75bp. The **publishedSI commercial alternative** uses20%formamide85°C10min, cooling, commercialAPOBEC100µl total and a slow4°C→50°Cramp. The main associates this alternative with ExtendedData7/8; no GBMSE75 or MiSeq-material-list item is silently assigned to its actual acquisition. The unoptimized v1-primer/Klenow workflow is not merged into either selected branch.

After deamination1.2xSPRI, qPCR consumes a separate analytical aliquot to choose indexing cycles; that material is not the preparative PCR template. KAPAUracil+ indexingPCR,0.8xSPRI and Qubit/Bioanalyzer aliquots precede known acquisition. Table2's10-cycle example and detailed adaptive cycles are retained as source scopes. SI elutes10µl before qPCR but later indexing table asks19µlDNA; exact scaling/dilution is unassembled, not silently created or declared a publication error.

## Ordering oligo and geometry limits

🟢 PublishedTable1 originalp2 prints v2 copy primer `GA/5mC/TGGAGTT/5mC/AGA/5mC/GTGTG/5mC/T/5mC/TT/5mC/GAT/5mC/`. Every displayedC is explicitly5mC. The exact body and typed methyl marks are modeled; no extra C/T is inserted to match a familiar TruSeq primer. Full5pyCY-adapter bodies, index members, complete Bst/product junctions and historical kit geometry remain unreviewed. The source guard originally exposed a literal `/5mC/` parsing failure. A narrowly scoped exact-token parser alias now preserves all seven nucleotide positions and methyl-C marks, and the actual source strand passes with no exemption. The chemistry, exact state and source text were preserved.

## Scope and ancestry

This two-original cappedpass used the published SI and associated main with finalVOR applicability explicit. Other publisherreporting/data files were fetched but not chemistry-read. The associated main explicitly invokes ACE-derived A3A conditions; the cited original is Wang et al.2021, DOI10.1007/978-1-0716-0876-0_27 (mainreference25), whose detailed chemistry remains unreviewed here. This supports a component-level ancestry candidate to that detailed protocol, not full recipe equality with the existing2018 ACE-seq defining entry. Similarity to TAPS/EM-seq does not establish recipe inheritance. No full review or accepted missing-information claim.

## Actual retrieval history

- https://pmc.ncbi.nlm.nih.gov/articles/instance/10763687/bin/NIHMS1948388-supplement-Supplementary_Tables_and_Note.pdf — PMC supplementary file; host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after reCAPTCHA challenge page (5.7 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (5.7 h left; not requested) (6.0 h left; not requested); the publisher's ESM files were fetched, so probably redundant

Availablepublishedsupplements were read where stated; failedPMC aliases and supersededpreprints are retrievalhistory, not requests to reacquire available/superseded originals. No unsupported source-absence conclusion.

## Local validation

Privategraph-onlyimplementation; no page builder or completionmarker. Focusedtests cover selected chemical identity/copying/fates. Fullhistoricalmolecular/end/readgeometry remains unassembled. Exact sequenceguard results and independent critical review will be recorded in the source-review JSON.


## Bounded publication review — 2026-10-30

Spike addition now explicitly changes no covalent molecule, preserving the three distinct control species. The optimized Bst copy reaction receives actual dATP, dGTP, dTTP and dmCTP contributions rather than hiding all nucleotides inside an enzyme label; the published v2 fully methylated primer remains exact and physically distinct from native genomic methylation. Actual β-GT/UDP-glucose, MBP-M.MpeI-N374K/CxSAM and branch-specific A3A or commercial APOBEC catalyst contributions are explicit. The optimized main A3A dose/conditions remain delegated to the actual ACE reference; no commercial APOBEC/formamide substitution is made. Index PCR consumes enzyme, nucleotide and buffer components from one actual named KAPA Uracil+ ReadyMix allocation, not fresh duplicate component additions. The SI states a 25µL master allocation; the optimized-main physical allocation volume is left unassigned rather than declaring full branch equivalence. Selected loading material and an untransferred production remainder follow actual analytical QC. The SI 10µL deamination elution versus 19µL indexing-DNA allocation remains unassembled and is not silently normalized. Original `pending` remains exact and uncovered: full 5pyC adapter/indexing-primer originals, ACE applicability and full copy/end/read geometry still require review. No copy-strand disappearance is inferred from poor amplification or computational filtering.

No additional originals were read or acquired during implementation. This capped pass accepts no missing-information claim. Existing full source obligations prevent publication.
