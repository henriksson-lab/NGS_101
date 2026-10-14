# cfMeDIP-seq — research notes

Authoritative protocol: Shen *et al.*, [10.1038/s41596-019-0202-2](https://doi.org/10.1038/s41596-019-0202-2).

## Evidence used

- 🟢 Naturally fragmented cfDNA is end-repaired, dA-tailed and adapter-ligated before immunoprecipitation.
- 🟢 Methylated and unmethylated spike-ins measure recovery; λ filler DNA permits 1–10-ng input immunoprecipitation.
- 🟢 Heat-denatured library molecules are enriched with a 5mC antibody, then PCR-amplified for paired-end sequencing.
- 🟡 The defining [Nature Protocols article](https://doi.org/10.1038/s41596-019-0202-2) names an adapter-ligation step, but the accessible article page does not establish a sequence-exact adapter and index-primer set. The protocol's genomic library shell is therefore a platform model pending review of the complete protocol and any vendor oligo manual for the specified kit version.

## Protocol boundary

cfMeDIP-seq is an enrichment assay and does not infer single-base methylation as bisulfite or enzymatic-conversion libraries do.
