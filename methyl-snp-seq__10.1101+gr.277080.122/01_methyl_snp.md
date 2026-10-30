# Methyl-SNP-seq

Yan et al.2022, DOI https://doi.org/10.1101/gr.277080.122 . Private source-qualified capped review; no finished schematic or accepted absence claim.

🟢 Published PMC9808626 main Methods148–164, published Supplemental Materials protocol steps1–7 and Supplemental Tables workbook selected inventory read. The original main PDF p2088 and SI p8 were rendered to verify names: main prints EndoVI (M0304); SI prints Endo IV (M0304). The selected explicit SI procedure controls the local model, with vendor identity unresolved. Source text also names methylated adapter E7536 in main versus E7165/E7140 in SI; no historical edition equivalence inferred.

🟢 Hairpins close both fragment ends; uracil excision/incision generates nicks. PolI translation with 5-methyl-dCTP replaces one original strand with a protected sequence copy, while displaced DNA products remain in the reaction until cleanup. Original-strand methylation and synthetic copy methylation have distinct meanings; the opposite original strand's methylation is not retained as an independently measured original methylome. Physical bisulfite uracil is distinct from amplified/read T. Source-known operations continue through PCR, optional target enrichment and sequencing; full hairpin/Y adapter/cycle geometry remains unassembled.

## Oligo and workflow

🟢 Main p2088 custom hairpin as written:

```
5′-(p)CCACGACGACGACGACGAGCGTTAGGCTCGTCGTCGTCGTCGUGGT-3′
```

The sequence is retained with5′phosphate and internal uracil in a DNA-oligo context. Bare U is not used to assert RNA sugar. Exact intramolecular stem/loop/nick and Y junction geometry remains unassembled; no source-unverified hairpin fold is drawn.

🟢 Human DNA is Covaris fragmented (WGS250bp; exome400–500bp), end-repaired/dA-tailed with historical NEB E7546, column recovered, hairpin ligated and0.9x SPRI recovered. ExoIII/ExoVII digest incompletely protected classes before another column recovery. UDG/printed SI EndoIV incision, PolI translation and Taq/dATP tailing occur by direct additions without intervening purification. Digestion/translation fragments remain until the actual column cleanup. Methylated adapter ligation and0.8x SPRI precede Fast Bisulfite conversion (SI95°C25min) and Q5U amplification (SI8cycles). The graph carries unknown kit/primer geometry explicitly through these operations.

🟢 Optional target enrichment starts with the indexed library pool, uses Twist Human Core Exome102025 and then NEBNext Q5 M0544/E7603 PCR. Q5 enrichment is distinct from the preceding Q5U conversion-library PCR. Human WGS/exome acquisition is NovaSeq6000 PE100; source bacterial comparison uses NextSeq550 PE150. The local default is the human WGS branch; target capture is a branch, not another method entry. The independent EM-seq validation recipe is not imported.

🟢 Additional printed catalogue boundary: main names Abcam ab117127, SI p8 prints ab11727; both are preserved pending historical vendor identity review. Main methylated adapterE7536 versus SI E7165/E7140 similarly remains unresolved as a supplier/revision issue. These differences do not justify inventing a resolved reagent identity. No publication-wide chemical error is inferred from catalogue-number similarity.

## Physical and computational limits

The selected `OriginalCopyPair` relation accepts an explicit original sequence and selected original5mC positions. Its newly copied opposite strand has synthetic5mC at every copiedC; its original-strand bisulfite product contains physicalU at unprotectedC. Its complementary DNA descendant is a distinct object. Synthetic test strings are expectations for this relation, not study DNA sequences or assembled sequencing libraries. No exact per-cycle pairing or hairpin orientation is inferred from matching primer motifs.

Source processing discards poor-quality cycle100 computationally; joint original/protected-copy observations support methylation and variant calls. Those readout operations do not erase DNA and do not restore the replaced original opposite strand's methylation. Historical NEB/Abcam/Twist guides and linked Read_Processing/README.pdf remain deferred within the cap; available Supplemental Tables were selected-inventory read, not fully numerical reviewed. Full source checklist/ref manifest records all scope limits. No source information is declared confirmed absent.
