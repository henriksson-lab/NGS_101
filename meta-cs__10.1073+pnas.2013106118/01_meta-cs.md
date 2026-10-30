# META-CS

Accurate SNV detection in single cells by transposon-based whole-genome amplification of complementary strands. PNAS (2021), [doi:10.1073/pnas.2013106118](https://doi.org/10.1073/pnas.2013106118).

🟢 Sources read: published PMC7923680 main71/173–195; whole DatasetS1 oligo workbook; SI Appendix FigS1–S4 captions. THREE originals, counted before opening; source investigation closed at cap. DatasetS2 and cited predecessor originals deferred. Mandatory scraper reviewed8candidatefamilies/2files; the underlying48oligos were checked in DatasetS1, not accepted from scraper output. Original SI strand graphics were not used to supply unknown end geometry.

🟢 META-CS derives from multiplexed end-tagging amplification: main71 calls this its previous WGA method; main177 cites ref8, Dip-C (10.1126/science.aat5641), specifically for META transposon DNA preparation. This supports a scoped component-adaptation relation, not inheritance of the Dip-C chromosome-contact assay. Carrier ssDNA is explicitly the LIANTI carrier (10.1126/science.aak9787); this establishes only that oligo, not LIANTI T7 chemistry. The graph assembles its own source-specific operations.

🟢 Sixteen tags are annealed individually5µM in10mMTrispH7.5/50mMNaCl/1mMEDTA, pooled equally, mixed equimolarly withTn5 atRT30min and stored125nM/−80°C. Tags supply finite fragment-specific tag combinations, not synthesized N-base UMIs. DatasetS1 contains16tag,16Adp1 and16Adp2 sequences. Representative literal oligos (all5′→3′; no modification specified in these rows):

```
META-CS-1 GGCACCGAAAAAGATGTGTATAAGAGACAG
META-CS-1-adp1 ACACTCTTTCCCTACACGACGCTCTTCCGATCTGGCACCGAAAAAGATGTGTATAAG
META-CS-1-adp2 GACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGCACCGAAAAAGATGTGTATAAG
carrier TCAGGTTTTCCTGAA
```

🟡 The local constructor derives48sequence inventories from the16published prefixes, canonical mosaic end and explicit source-compatible Illumina arm intervals. Every inventory matches local source text. Finite tag identity is structured as molecular-complex barcode; it does not establish a random-base UMI or guaranteed unique molecule identity. Tn5 loaded duplex partners, genomic ends and final complete constructs remain unassembled.

🟢 One mouth-pipetted/FACS cell is lysed in2µlbuffer at50°C1h/70°C15min; carrier500nM and Qiagenprotease1.5mg/ml are present. Source optional storage is−80°C. Add8µltranspositionmix for0.38nMTn5,55°C10min. Add2µlstopbuffer/Qiagenprotease,50°C30min/70°C15min. The source literally prints300nMNaCl in this stopbuffer; it is retained without normalizing to300mM.

🟢 First strand-tagging adds13µlQ5mix including0.8µl50µMtotalAdp1;72°C3min/98°C30s/62°C5min/72°C1min. Add1µlExoI37°C30min/80°C20min. Second-tagging adds4µlQ5mix including1µl50µMAdp2;98°C30s/62°C5min/72°C1min. Another1µlExoI37°C30min/80°C20min follows. ExoI cleavage products remain in the same tube; it is not a physical cleanup. PCR adds19µl including5µlNEBMultiplexUniversalPrimer and5µlNEBIndexPrimer;98°C20s then10cycles (11haploid)98°C10s/72°C2min,72°C2minfinal. Exact historicalNEBprimer/index identities remain unassembled.

🟢 DCC5 cleanup precedes separateBioanalyzerQC and pooling of single-cell libraries. Sequential AMPures retain long fractions at0.6×, middle fractions by adding0.15× to the firstsupernatant(total0.75×), then short fractions by another0.15×. These are different size fractions, not three equal aliquots. The last supernatant's fate is unassigned. Long fractions useHiSeq2500PE250; middle/short useHiSeq4000PE150. Short fractions were sequenced only for340-2/3/4/5/13/14/15, so their selected-cell fate remains unresolved here. Source20%PhiX is a physical loading addition. Strand-agreement SNV filtering is computational, not physical removal or proof of perfect strand specificity.

Open scope: source-specific Tn5/end/primer/index geometry, DatasetS2 run assignments, precursor originals and complete intermediate molecules. This is a private incomplete graph, not complete validation or accepted missing-information review.
