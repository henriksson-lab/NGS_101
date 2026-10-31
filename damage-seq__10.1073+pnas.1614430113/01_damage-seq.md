# Damage-seq: sequence a polymerase-stall boundary

Hu et al. (2016), [Cisplatin DNA damage and repair maps of the human genome at single-nucleotide resolution](https://doi.org/10.1073/pnas.1614430113).

## Capped evidence and current boundary

One reader,28October2026; published main article PMC5068337, Results developing Damage-seq, Figure1 caption, lesion-coordinate interpretation and Materials/Methods. The full main HTML/text is available; generic EuropePMC XML returnedHTTP500. **SIAppendix pnas.1614430113.sapp.pdf remains unavailable through the fetched PMC download gate**, with exact locator in externalMANIFEST.tsv. No challenge bypass attempted. Mandatory full-bundle scraper returned no exact oligos; this is not evidence that the missing SI lacks sequences. No abstract-only claim is substituted for a detailed chemistry review: this pass reads the main full text but remains limited by the chemistry-relevant SI.

## Positive main-text workflow

Main Methods reports GM12878 cells,200µMcisplatin/oxaliplatin and1.5/3h treatments; exact drug/assay/run matching remains to be checked in SI rather than inferred from the compressed sentence. The separate XR-seq assay has additional chemistry and is not the same library. Main Results/figure establish the following order:

1. Isolate genomic DNA, sonicate, end-repair and ligate the first adapter.
2. Denature and immunoprecipitate Pt-adduct-bearing single-stranded DNA using anti-Pt antibodies. No perfect selection of damaged strands is implied.
3. Anneal a biotinylated primer to adapter1 and extend with NEBNextQ5. Bulky lesions block extension. The **new primer-derived strand stops adjacent to the lesion; the Pt adduct remains on the original template**, rather than being copied into the new DNA or chemically removed by polymerase. Undamaged/full-length extension products can also occur.
4. Streptavidin purification enriches biotinylated extension products, then adapter2 is ligated and indexedPCR forms the library. Exact oligo sequences, biotin position, ligation ends, strand-separation conditions and PCR geometry remain unassembled pendingSI; these are not claimed unpublished.
5. Source reserves1%ligation products for twoPCR quality tests. Different primer sets distinguish products associated with full-length/unblocked extension versus all products. Exact assay primers and library run mask remain unresolved in this pass.

## Molecular fate versus computational filtering

The model retains both stalled and unblocked extension classes through the source library endpoint. Full-length/undamaged-derived reads are discarded **after sequencing** when the read begins with adapter1's5′sequence; computational filtering is not a physical molecular discard during ligation. The source places the cisplatin lesion immediately5′upstream of the aligned read start, enrichingGG(and lesserAG) contexts. The reported boundary is an inferred lesion location, not the chemical identity of a nucleotide in the new strand. No fixed genomic position or universal single-nucleotide blocking efficiency is assumed.

The source notes this implementation does not detect rare Pt interstrand crosslinks. Its prediction of applicability to other bulky lesions is not equivalent to an experimentally modeled conversion for every adduct class.

## Related branch and deferreds

The same paper modifies UVXR-seq for platinum adducts by reversingPt withNaCN before amplification. This is an explicit **Pt-XR branch** inheriting components from UVXR-seq, not evidence that Damage-seq inherits the TFIIH capture/adapter/photolyase recipe. No constructor is imported from XR-seq.

Deferred: SIAppendix original and its linked chemistry sources, exact adapter/primer polymers and reagent/end details, molecular binding/stall/strand-retention validation, run configuration, and independent second reader. The partial graph is source-positive preparation order with unresolved chemistry, not a completed protocol review or missing-information acceptance.

## 2026-10-30 public preparation review

Reinspected the same published main original: Results developing Damage-seq, Figure 1 caption, lesion-coordinate interpretation and Methods. The bundle scraper was rerun; zero hits do not establish missing oligos. The chemistry-relevant SI access hold remains unchanged. The public page presents the positive main-text preparation with unresolved end chemistry, strand separation, exact oligos and run settings.

The extension now has an explicit unassigned molecular reagent supply. The PCR supply is conceptually partitioned into catalytic, nucleotide-bearing and nonincorporated portions of one unresolved supply; this is not a source-disclosed formulation, fresh extra additions or a verified kit recipe. Both stalled and unblocked molecular classes remain through library formation. The reserved 1% ligation-product QC aliquot stays analytical, and adapter-1-starting background read removal remains post-sequencing. No new absence claim or complete review is asserted.
