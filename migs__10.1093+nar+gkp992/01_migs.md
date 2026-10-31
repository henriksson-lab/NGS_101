# MiGS

Serre et al. (online 2009; issue 2010), MBD-isolated Genome Sequencing provides a high-throughput and comprehensive survey of DNA methylation in the human genome. DOI: https://doi.org/10.1093/nar/gkp992

MiGS enriches methylated DNA using recombinant MBD2's methyl-CpG-binding domain, followed by ordinary genomic sequencing. The sequencing letters do not undergo bisulfite conversion in the represented MiGS preparation; enrichment is not a base-resolution methylation call.

## Sources read

🟢 Published open-access article PMC2811030 (manifest: author manuscript False), Materials and Methods, DNA preparation / Isolation of methylated DNA by recombinant MBD / Library preparation and sequencing / Mapping of sequence reads, text lines 38–52. Online publication is 11 November 2009; the journal issue is January 2010. The selected defining DOI is unchanged.

🟢 Published supplement `gkp992_1.pdf`, Supplemental Figure 1–5 captions, Supplemental Table 1 header and primer rows, Table 2 context. Table 1 explicitly contains PCR primers for validation bisulfite sequencing. Those oligos are not library adapters or MiGS amplification primers.

## Preparation

1. 🟢 Purify genomic DNA with Proteinase K/phenol:chloroform/ethanol. Sonicate to approximately 150–600 bp in BW buffer. QIAquick purification is explicitly used to exclude fragments smaller than 100 bp. Use 5 µg cleaned fragmented DNA per precipitation.
2. 🟢 Prepare affinity reagent with 41.7 µl/ml protein G magnetic bead slurry, 10 µg/ml anti-His RGS antibody, 1 µg/ml self-ligated pCR4-TOPO blocking DNA and 160 nM His-tagged recombinant MBD2 domain in BW buffer, overnight at 4°C. Physically wash the bead conjugates with one volume of cold 1× BW buffer, remove that wash fraction, then resuspend them in one volume of fresh 1× BW buffer before adding genomic DNA. This physical conjugate wash is distinct from the five post-capture washes. The blocker is a separate molecular species, never inserted into genomic strands or represented as an index.
3. 🟢 Incubate sheared DNA with conjugates overnight at 4°C, and wash five times in BW buffer. Captured DNA is enriched by the domain's methyl-CpG affinity; the model does not assume every captured fragment is fully methylated or every unbound fragment unmethylated.
4. 🟢 SDS/Proteinase K elution at 60°C for 2 h, phenol:chloroform extraction, isopropanol precipitation; resuspend in 45 µl water. 🟡 Possible residual blocking-vector DNA is retained separately through Proteinase K release and the extraction/precipitation workup: protease does not digest DNA, and no complete blocker separation is established here. Exact blocker partition, recovery and eventual library/read contribution remain unresolved; this is a conservative possible-carryover state, not a source-measured contamination claim. The article prints “20 ml Tris (pH 8)” in its elution buffer composition. That wording is retained as a source anomaly/open formulation question; the model does not silently change it to a concentration.
5. 🟢 Process 10 ng MBD-isolated DNA with the original Illumina ChIP-Seq Sample Prep Kit following manufacturer instructions. Libraries have a tight insert distribution near 120 bp. GAII generates 36 bp reads.
6. 🟢 Unique alignment and collapsing identical strand/orientation start positions occur computationally. These are not pre-sequencing molecular depletion or UMI processing.

## Model and source boundaries

The partial graph owns the domain/antibody/protein-G capture sequence and separates physical conjugate/capture wash waste from computational duplicate handling. Possible blocking DNA remains a distinct molecular species through release and recovered DNA, and remains an unresolved possible component at the kit boundary. Exact historical kit components, adapter polymers, amplification details and sequencing-primer landing sites remain unassigned. No later Illumina or antibody-MeDIP recipe is substituted. Fragment methylation state is an enrichment interpretation, not a C/T chemical readout.

The cited expression/purification sources (ref22, COMPARE-MS, DOI 10.1093/nar/gnj022, and ref26, DOI 10.1074/jbc.M505593200) are direct linked originals, deferred after the defining article and chemistry-relevant supplement. Supplementary index export failed validation, while the identified gkp992_1.pdf was obtained and read. Neither failure nor scraper outcome establishes missing information. Full historical vendor library instructions and an independent reader remain open.

## Lineage

Explicit component ancestry: the article says recombinant MBD2 domain expression/purification follows refs22/26. That does not establish inheritance of an NGS library recipe, nor identity with MIRA's MBD2b/MBD3L1 complex. The pCR4-TOPO blocker is a reagent; it is not evidence that a cloning protocol is the parent NGS chemistry.

## Bounded publication assessment, 30 October 2026 — source hold

Private MiGS model retains original washed conjugate provenance and possible independent blockingDNA species. Actual gelQC, selected5µggenomicDNA/5mlwashedconjugate/10ngkit contributions, residualstocks, BioAnalyzerQC and selectedloadingstock are explicit. The mixed washed-conjugate STOCK is molecular material before an actual5ml aliquot is assigned selection-reagent function. Aggregate source kit oligo/nucleotide/catalyst inventory classes remain unassembled in identity/dose/packaging, not imported Illumina polymers or an inferred specific PCR recipe. Literal20mlTris is preserved after verifying actual publishedPDFp2, not repaired as an extraction artifact.

The required SourceItem unreviewed_original remains unchanged: the manufacturer protocol actually cited for the Illumina ChIP-Seq Sample Prep Kit remains relevant/unreviewed. Its study-used revision and direct documentURL remain unassigned, not a fabricated unavailable disclosure document or an observed download failure. No excluded/not-applicable coverage was added.

The normal publication gate still rejects the source obligation. The generated page is a private source-scope preview, with no protocol-complete marker; this pass does not recommend catalogue promotion. Source acquisition was not expanded.
