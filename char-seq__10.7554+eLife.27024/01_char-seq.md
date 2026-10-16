# ChAR-seq — research notes

## Sources read

- 🟢 Bell et al. 2018, DOI `10.7554/eLife.27024`, original eLife v2 Figure 1—figure supplement 1 and caption (source image; JATS lines 88–90). The same bridge strings and modifications are printed in the [original bioRxiv supplement](https://doi.org/10.1101/118786), PDF lines 36–40 and 101–109.
- 🟢 Limouse et al. 2020, PMCID `PMC7607860`, *Methods in Molecular Biology* §2.3 (JATS lines 132–134), prints a **later bridge variant** with a three-base UMI and different geometry.

## Construct path

🟢 Crosslinked chromatin-associated RNA is fragmented and ligated in situ to the RNA-compatible end of a biotinylated bridge. The bridge’s opposite end is joined to nearby genomic DNA. RNA is converted to cDNA, material is fragmented and bridge-bearing junctions are affinity enriched before single-end library construction. The internal bridge provides the orientation used to distinguish RNA-derived from DNA-derived sequence.

## Bridge geometry

🟢 The defining 2018 figure prints the bridge top as `5′-App-AAAAAAACCGGCGTCCAAGGATCTAG-3′-Sp3`. The first six A bases form the single-stranded RNA-ligation tail; the seventh A enters the duplex. The activated 5′ end ligates specifically to fragmented RNA 3′-OH ends. This version has **no NNN UMI**. The 3′ Sp3 blocks circularization and is removed when DpnII opens the DNA-facing end.

🟢 The 2018 figure draws the bottom strand left to right as `3′-TGGCCGCAGGTTCCTAGATCCTAG-5′`, which must be reversed to order it 5′→3′ as `GATCCTAGATCCTTGGACGCCGGT`. Twenty bases pair with the top; the bottom has a four-base 5′ `GATC` extension. The figure labels an **internal desthiobiotin-TEG** attachment, but does not establish its exact nucleotide or internucleotide bond. The model marks the paired bottom segment as modified without assigning the attachment to a particular base.

🟢 **Version difference visible on the public page:** the later 2020 *Methods in Molecular Biology* §2.3 instead prints top `/5rApp/AANNNAAACCGGCGTCCAAGGATCTTTAATTAAGTCGCAG/3SpC3/` and bottom `/5Phos/GATCTGCGACTTAATTAAAGATCCTTGGACGCCGG/iBiodT/T`. Its `NNN` UMI, PacI segment and `/iBiodT/` belong to that later implementation. This page and model use the defining 2018 figure; neither source shows that the later oligos were used in the 2018 experiment.

🟢 DpnII digestion removes the blocked distal bridge piece and cuts nearby genomic DNA. T4 DNA ligase then joins compatible DpnII cohesive ends. RNase H and DNA polymerase I complete the cDNA-containing duplex, and sonication followed by streptavidin capture enriches only molecules retaining the bridge.

🟢 The defining paper used NEBNext end repair/dA tailing, NEBNext hairpin adapters, USER opening and TruSeq-compatible indexing, then single-end 152-base sequencing (eLife v2 JATS Methods lines 159–161 and extended protocol lines 774–864). The source names NEBNext Multiplex Oligos for Illumina but gives no catalogue number, kit revision or selected index. 🟡 E7335/E7600 vendor manuals document candidate family strands; none establishes which member was used in 2018. The page therefore marks the standard TruSeq-compatible shell inferred and leaves primer landing and exact index assignment unresolved.
