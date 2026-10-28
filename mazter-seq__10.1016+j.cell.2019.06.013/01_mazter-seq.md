# MAZTER-seq

Garcia-Campos et al., *Deciphering the “m6A Code” via Antibody-Independent Quantitative Profiling*, Cell (2019), [defining article](https://doi.org/10.1016/j.cell.2019.06.013).

## Sources inspected

🟢 The [institutional copy of the publisher PDF](https://weizmann.elsevierpure.com/ws/portalfiles/portal/111306132/ss_Cell_DecipheringThe-m6ACode_PV2019.pdf) identifies itself as the version of record and prints the defining DOI, article, STAR Methods and supplemental figure captions. Read Results p732 (text lines188–219), discussion limitations (2755–2784), the resource-table oligo referrals (3148–3183), Methods “mRNA preparation”, “MazF digestion”, “RNA barcoding”, “mRNA m6A immunoprecipitation”, and “Library preparation” (3330–3376), sequence-bias controls (3377–3386), and cleavage-efficiency analysis (3464–3486). These are specific inspected sections, not a complete source review.

The primer scraper was run on the entire cached text. Its validation-assay hits were inspected in their Methods context; they are not the genome-wide library adapters. The library polymers are explicitly referred to Table S7. The [published article/supplement index](https://www.sciencedirect.com/science/article/pii/S0092867419306762) returned HTTP403 through the standard fetcher. The exact Table S7 attachment URL remains unresolved. A blocked index does not establish that the table is unavailable everywhere or that its information was omitted.

The linked preprint is superseded and excluded. No preprint oligos are substituted. External originals and the retrieval manifest remain outside this checkout; see [retrieval manifest](ref/MANIFEST.md).

## What is new

🟢 MazF preferentially cuts immediately before an unmodified ACA motif; first-position m6A protects that site. The readout is fragment-start/end and read-through coverage, rather than a chemically converted sequencing base. The source reports incomplete, sequence/structure-dependent ACA cleavage and low non-ACA cleavage. Thus an uncut ACA is not by itself a methylation call, and the model does not predict all cuts or molecular fractions.

🟡 A local selected-event constructor validates an internal ACA, returns the two RNA substrings without loss or reordering, and rejects calling a protected first-A m6A site that susceptible event. It deliberately assigns no newly formed phosphate/hydroxyl chemistry. The drawn fragments are symbolic; they do not represent a study transcript. Uncut and protected RNA remains in the physical mixture.

## No-IP library route

1. 🟢 Enrich poly(A) RNA from extracted total RNA using the source input-dependent mRNA enrichment kit. Heat-denature, then digest with source MazF at the printed incubation conditions. Recover RNA with RLT/ethanol/MyOne SILANE cleanup. Purification does not prove that all uncut RNA is removed.
2. 🟢 Treat with FastAP, T4 PNK and Turbo DNase, then clean up. The source calls the product dephosphorylated RNA. Exact MazF terminal chemistry and the complete end-repair mechanism remain unassembled here.
3. 🟢 Combine repaired RNA with the RNA adapter named 3iLL, heat and snap-cool together, then add T4 RNA ligase1. This thermal step does not establish an annealed duplex geometry. Purify and pool source barcoded samples. Exact barcode identities, modifications and adapter polymers await Table S7. No UMI is inferred from the word barcoding.
4. 🟢 For this no-IP route, proceed directly to library preparation. An optional branch instead allocates input/IP material and performs two antibody-IP rounds; that branch is not represented in this graph.
5. 🟢 Prime cDNA synthesis from the RNA adapter with the source rTd DNA primer and SuperScript III. Omit heat inactivation to preserve the RNA/cDNA hybrids. Treat residual primers with ExoSAP-IT, stop with EDTA, then hydrolyse RNA with NaOH and neutralize. RNA/primer-derived products remain molecular material until the physical cDNA cleanup.
6. 🟢 Ligate the source 5iLL-22 DNA adapter to the cDNA 3′ end with T4 RNA ligase1; clean up. Amplify with source primers F/R and KAPA HiFi HotStart ReadyMix for {{= mazter.PCR_CYCLES }} cycles, then AMPure-clean and assess the library.

The graph follows the source operations through the library but retains unresolved ligation, primer and final-product geometry. No historical adapter recipe is imported by similarity. It has no finished-protocol marker and remains notes-only.

## Sequencing and interpretation

🟢 Published Results specify paired-end sequencing; Methods name NextSeq500. Exact run cycles, physical sequencing primers, index polymers and their acquisition orientation remain unassigned. ACA-site start/end versus read-through counting is computational; it does not physically remove uncut molecules. The quantifiable site-distance threshold is library-dependent in the published analysis, not a universal fixed biological cutoff.

## Deferred questions

- Obtain and inspect Table S7 library/IVT oligos verbatim with modifications, distinguish assay-specific primers, and assemble exact ligation-ready ends, RT hybrid and PCR product.
- Inspect the directly cited Imanishi et al.2017 MazF methyl-sensitivity original and historical reagent documentation for the terminal chemistry; no mechanism is asserted from the reagent name alone.
- Inspect remaining original figures/tables and the optional IP/calibration routes and their source-linked predecessors. Results describe the synthetic calibration RNA as88nt, whereas the IVT Methods describe103nt; the model does not select a value or reconstruct those controls before the Table S7 comparison.
- Reconstruct run primers, actual index membership and cycle layout. A second independent source/linked-material review has not occurred; no missing-information claim is accepted.
