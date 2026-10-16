# ORGANIC profiling — research notes

Kasinathan et al., *Nature Methods* (2014),
[doi:10.1038/nmeth.2766](https://doi.org/10.1038/nmeth.2766).

## Sources read

- 🟢 The defining article and Supplementary Figure 1 workflow (PMCID: PMC3929178).
- 🟢 Orsi et al., *Mapping regulatory factors by immunoprecipitation from native
  chromatin*, Current Protocols (2015),
  [doi:10.1002/0471142727.mb2131s110](https://doi.org/10.1002/0471142727.mb2131s110),
  PMCID: PMC4410783. This later same-author stepwise ORGANIC procedure publishes a
  concrete short-insert library oligo implementation. It was **not** cited by the 2014 article.
- 🔴 The defining article cites Henikoff et al. (2011), [doi:10.1073/pnas.1110731108](https://doi.org/10.1073/pnas.1110731108), reference 19, for its modified library preparation. The predecessor original (PMCID: PMC3215028) and linked 4.3-MB Supporting Information remain unavailable for full independent checking in the current source bundle. Its search-indexed abstract and Results describe short-fragment paired-end library preparation but are not a substitute for reading the Methods and supplement.

## Protocol boundary

ORGANIC means *Occupied Regions of Genomes from Affinity-purified Naturally Isolated
Chromatin*. It is not ChEC-seq: MNase is free in the nuclear digest rather than fused to
the factor, and the factor–DNA complex is selected afterward by immunoprecipitation.

1. 🟢 Nuclei are isolated without formaldehyde. Low-salt, cold conditions preserve native
   noncovalent protein–DNA complexes, and no sonication boundary is introduced.
2. 🟢 MNase concentration and duration are titrated. The recommended QC ladder is at least
   roughly 75% mononucleosome-sized DNA with a clear smear of shorter fragments while
   retaining small amounts of di- and trinucleosomes. Light digestion retains local
   accessibility and sequence bias; excessive digestion can damage the complex.
3. 🟢 MNase first nicks/cuts exposed DNA and then trims exonucleolytically until a bound
   protein blocks it. A single transcription factor commonly protects 20–40 bp, whereas
   jointly occupied regions can protect 60–80 bp and nucleosomes protect longer fragments.
4. 🟢 Digested chromatin is solubilized at low ionic strength. The epitope-tagged factor
   or endogenous target is affinity-purified with its still-associated DNA. Input DNA from
   the total digest measures MNase/accessibility background.
5. 🟢 Proteinase treatment and DNA purification recover 3–100 ng of footprint DNA. There
   is no crosslink-reversal step because the input chromatin was never fixed.
6. 🟢 End repair, 3′ dA addition, indexed-adapter ligation and KAPA PCR use the authors'
   low-input paired-end workflow, which deliberately retains very short inserts.

## Exact library oligos

The later 2015 same-author ORGANIC procedure publishes this implementation:

- Universal Adapter:
  `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`
- Indexed Adapter:
  `5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′`
- P5 enrichment primer: `5′-AATGATACGGCGACCACCGAG-3′`
- P7 enrichment primer: `5′-CAAGCAGAAGACGGCATACGAG-3′`

🟢 `*` is a terminal phosphorothioate linkage, `/5Phos/` is the indexed adapter's 5′
phosphate, and `NNNNNN` is the six-base i7 sample index. The construct is single-indexed.
🟡 The 2014 article's Supplementary Fig. 1c names TruSeq adapter ligation, but does not
print the oligo bases or identify which six-base index was used for each sample. The
2011 library predecessor cited by that article still needs its own source check before
these later exact arms can be attributed to the 2014 experiments.

## Why paired ends matter

🟢 The defining work used 25-cycle paired-end HiSeq 2000 sequencing. Each pair reconstructs
the complete MNase-bounded fragment. Analyses that collapse mates into independent read
starts lose the distinction between compact transcription-factor footprints, co-occupied
sites and nucleosome-sized fragments. The method paper therefore treats fragment length
and both endpoints as molecular observations, not merely library QC.

## Interpretation limits

🟢 ORGANIC avoids formaldehyde and sonication artifacts and can give base-pair-scale direct
binding maps. Its tradeoff is dependence on a complex surviving native preparation and
affinity purification. Weak/fast-exchanging binding can be lost; MNase base preference,
digestion depth and incomplete solubilization remain possible biases. Replicates, input
digests and digestion-series agreement are the appropriate controls.
