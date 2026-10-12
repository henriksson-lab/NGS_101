# DR-seq — research notes

## Sources read

- 🟢 Dey et al. 2015, DOI `10.1038/nbt.3129`, PMC full text, Online Methods and supplement.
- 🟢 Teichlab DR-seq reconstruction, used as a secondary check for the complete Ad-1x sequence and construct order.
- 🟢 Zong et al. MALBAC paper, DOI `10.1126/science.1229164`, for the looping/quasilinear mechanism.

## Shared amplification before the split

🟢 Ad-1x is `CGATTGAGGCCGG` + T7 promoter `TAATACGACTCACTATAGGG` + RA5 `GTTCAGAGTTCTACAGTCCGACGATC` + an 8-nt cell barcode + dT24V. It reverse-transcribes poly(A) RNA while genomic DNA remains in the same lysate.

🟢 Ad-2 is `GTGAGTGATGGTTGAGGTAGTGTGGAG` + N8. Seven quasilinear cycles randomly prime both genomic DNA and cDNA. Most short genome-derived products have Ad-2 at both ends; productive transcript-derived products have Ad-2 at one end and retain Ad-1x at the other. Full Ad-2/Ad-2 products can self-anneal through complementary ends, suppressing their continued amplification.

🟢 The random priming coordinate—not a conventional UMI—is used as a length-based molecular identifier for RNA deduplication. The old page incorrectly presented the 8-nt Ad-1x cell barcode as an 8-nt UMI.

## Branch-specific finishing

🟢 The pool is split only after the shared amplification. The RNA aliquot receives P1 (`CGATTGAGGCCGGTAATAC`) for second-strand synthesis across the Ad-1x/T7 end, then 13 hours of T7 IVT. Genomic Ad-2/Ad-2 products lack that transcript-selective T7 architecture.

🟢 The DNA aliquot is exponentially amplified with P2 (`GTGAGTGATGGTTGAGGTAGTGTGGAG`), re-amplified with 5′-biotin P3 (`GTGAGCTGGAGTTGAGGTAGTGTGGAG`), sonicated, and depleted of the biotin-bearing Ad-2 ends on streptavidin beads before NEBNext Ultra library preparation.

🟡 The paper refers the RNA branch to the previously published CEL-Seq small-RNA finishing workflow, and does not print those outer library oligos in its accessible primary-source package. Likewise, the particular NEBNext library oligos are named by kit rather than printed. The public page therefore marks the final Illumina shells inferred while keeping all source-printed inner intermediates exact.
