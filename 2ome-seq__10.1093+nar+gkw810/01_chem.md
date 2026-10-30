# 2OMe-seq: Method1 without fragmentation

Incarnato et al., **High-throughput single-base resolution mapping of RNA 2΄-O-methylated residues**, DOI10.1093/nar/gkw810.

## Source scope and selected branch

🟢 Published PMC5388417 mainMethods70/Fig1 and the published supplement File006, Method2 lines215–277, were inspected. Two originals, investigation closed. Mandatory source-bundle scraper ran (54candidates/3files), with only selected original pointers adopted. The supplemental comparison is a different preparation route, not another protocol entry or a source of convenient steps to merge into Method1. Other supplementary tables/analytical methods remain uninspected. No full linked-source or absence audit is claimed.

🟢 Method1 compares AMV reverse transcription at high1mM versus low0.004mM dNTP, using random hexamers bearing a common sequencing handle. RNA degradation and cDNA-end library preparation preserve RT-stop positions. These effects are source-described enrichment/readout behavior, not a complete or perfectly quantitative modification-conversion claim. The two dNTP populations and their molecular contribution are explicitly separate from the AMV enzyme. The source's per-nucleotide-versus-total dNTP convention is not supplied by assumption.

## Source oligos and known local ends

🟢 RT primer: **`CCTTGGCACCCGAGAATTCCANNNNNN`**. The fixed body is the verified reverse complement of canonical TruSeq SmallRNA RA3. The six Ns are random priming degeneracy, not a UMI.

🟢 DNA adapter: **`5′-P-GATCGTCGGACTGTAGAACTCTGAAC-3′-C3 spacer`**. Its fixed body is the verified reverse complement of canonical SmallRNA RA5. CircLigaseII ligates its5′ phosphate to cDNA3′OH; the adapter's C3 terminus is a non-base blocker, not an extra C nucleotide or an unmodified3′OH. A partial Scene declares the known single-strand adapter5′P/cDNA3′OH interface. Complete cDNA sequences, sealed products and primer landing remain unassembled.

## Source sequence and physical fates

🟢 Selected spike-containing Method1 case:2µgRNA plus100ng in-vitro RNA spike mixture →20µlAMVRT0.5U/µl at42°C30min under either final dNTP condition →95°C10min →addRNaseH1µl and RNaseCocktail1µl/37°C30min →AMPure recovery →10%TBE-urea PAGE70–300nt →passive diffusionTE37°C16h →ethanol precipitation →CircLigaseII200U65°C6h with the modified DNA adapter →15cycle TruSeq PCR →AMPure →normalized aliquots for an equimolar pool →NextSeq500. Spike mixture sequences/proportions and actual sample indexes remain unassigned locally.

Actual unretained workup and gel material, residual gel recovery, sourceRNA stock and untransferred normalized-PCR stocks are recorded. Equimolar pooling does not implicitly consume the entire high/low PCR stocks; selected aliquots have unspecified normalization volumes, and remainders retain unresolved fates. Primer/enzyme/buffer components come from one source PCR bulk addition, avoiding double addition of both stock and free primers.

## Separate Method2 and remaining limits

🟢 Supplement215–277 describes fragmentation95°C5min in FirstStrandBuffer, T4PNK end repair, preadenylated3′RNA-adapter ligation,200–250nt RNA selection, adapter-specific RT and50–150nt cDNA selection, followed by the same class of DNA-end ligation and **18PCRcycles**. Its RNA adapter, gel thresholds and cycle count are not substituted into Method1. The source's RNaseH and alkaline-hydrolysis comparison experiments also remain separate.

Source names NextSeq500; precise run/read/index masks and full historical TruSeq PCR primer orderings are not assigned. No modern standard primer or whole-library endpoint is invented from a kit name. Complete target/cDNA/adapter/PCR structures, spike composition and linked CIRS-seq/vendor versions remain unassembled or unreviewed, not confirmed absent. No public schematic or finished status is claimed.
