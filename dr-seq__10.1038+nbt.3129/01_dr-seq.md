# DR-seq — research notes

## Sources read

- 🟢 Dey et al. 2015, DOI `10.1038/nbt.3129`, PMC full text, Online Methods and supplement.
- 🔴 The original CEL-Seq article and its supplement are not in the checked DR-seq source bundle. Teichlab's DR-seq reconstruction is useful for orientation, but it is not source evidence for Ad-1x bases or construct order.
- 🟢 Zong et al. MALBAC paper, DOI `10.1126/science.1229164`, for the looping/quasilinear mechanism.

## Shared amplification before the split

🟡 The working Ad-1x model is `CGATTGAGGCCGG` + T7 promoter `TAATACGACTCACTATAGGG` + RA5 `GTTCAGAGTTCTACAGTCCGACGATC` + an 8-nt cell barcode + dT24V. The DR-seq paper names Ad-1x and describes its cell-specific barcode and reverse transcription role, but it does not print this complete oligo. Its bases and length require verification in the cited original CEL-Seq material; the upstream reconstruction does not establish them.

🟡 The DR-seq paper prints `GTGAGTGATGGTTGAGGTAGTGTGGAG` as P2, which primes Ad-2-derived products. The working Ad-2 model adds N8, inferred from the cited MALBAC-style random priming rather than a complete Ad-2 oligo printed in DR-seq. Seven quasilinear cycles randomly prime both genomic DNA and cDNA. Most short genome-derived products have Ad-2 at both ends; productive transcript-derived products have Ad-2 at one end and retain Ad-1x at the other.

🟢 The random priming coordinate—not a conventional UMI—is used as a length-based molecular identifier for RNA deduplication. The old page incorrectly presented the 8-nt Ad-1x cell barcode as an 8-nt UMI.

## Branch-specific finishing

🟢 The pool is split only after the shared amplification. The RNA aliquot receives P1 (`CGATTGAGGCCGGTAATAC`) for second-strand synthesis across the Ad-1x/T7 end, then 13 hours of T7 IVT. Genomic Ad-2/Ad-2 products lack that transcript-selective T7 architecture.

🟢 The DNA aliquot is exponentially amplified with P2 (`GTGAGTGATGGTTGAGGTAGTGTGGAG`), re-amplified with 5′-biotin P3 (`GTGAGCTGGAGTTGAGGTAGTGTGGAG`), sonicated, and depleted of the biotin-bearing Ad-2 ends on streptavidin beads before NEBNext Ultra library preparation.

🟡 The paper refers RNA library finishing to CEL-Seq and names NEBNext Ultra for DNA preparation. The cited CEL-Seq original, its supplement, and a version-matched NEBNext guide remain unreviewed. The full Ad-1x and Ad-2 oligos and final Illumina shells therefore remain provisional; only P1, P2 and P3 are exact source-printed primers in the checked DR-seq article.
