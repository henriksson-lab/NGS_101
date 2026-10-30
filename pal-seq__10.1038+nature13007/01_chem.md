# PAL-seq — retained-tail sequencing-template preparation

Subtelny and colleagues, *Poly(A)-tail profiling reveals an embryonic switch in translational control*, Nature (2014), DOI [10.1038/nature13007](https://doi.org/10.1038/nature13007).

This capped pass uses the associated accepted manuscript PMC4086860, Methods lines425–482 and rendered original p8, as explicitly qualified preparation evidence. The final first-party Bartel-lab PDF was requested through the declared standard route but failed certificate verification (`curl error 60`); no insecure transport or challenge workaround was used. Final main-version applicability remains unreviewed. This is a private incomplete graph, not a public schematic, complete source audit, or accepted absence review.

🟢 The method has a positive preparation boundary separate from post-sequencing fluorescence. An all-DNA 5′-phosphorylated, 3′-biotin capture adapter is splint-ligated to RNA poly(A) ends. Partial RNase T1 digestion, extraction/precipitation and gel recovery precede streptavidin capture. On beads, a 3′-phosphatase-deficient T4 PNK phosphorylates RNA 5′ ends; T4 Rnl1 adds a distinct mixed DNA/RNA 5′ adapter. SuperScript II then copies the retained native poly(A) tract and upstream RNA into cDNA. Base hydrolysis liberates cDNA, which is gel-purified, denatured, neutralized and applied directly to a flow cell. A preparative PCR step is not invented.

This retained-tail copying and direct flow-cell preparation differ from the 3P-seq tail-only RT/RNase H tail-shortening/release path. Shared analytical goals or an RNase T1 reference do not authorize importing the full 3P-seq recipe. Organisms, developmental stages and fluorescence settings are branches, not extra methods.

🟢 Selected associated-main oligos with physical modifications:

| Component | As written |
|---|---|
| Capture DNA | `p-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGACACATAC-biotin` |
| DNA splint | `TTCCGATCTTTTTTTTT` |
| 5′ adapter, uppercase DNA/lowercase RNA | `C3.spacer-CAAGCAGAAGACGGCATACGAGTTCAGAGTTCTAcaguccgacgauc` |
| RT primer | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACG` |

Mandatory scraper output was verified in the original. In the trace-control RNA, an asterisk denotes an internal 32P label, as the text explicitly says; the scraper's phosphorothioate interpretation is not adopted. The model preserves the 5′ adapter's C3 block and physical RNA/DNA arms rather than flattening it into an ordinary DNA oligo.

🟡 A selected-copy helper computes the full dT tract and upstream reverse complement from an explicit native RNA body and tail length. This conditional sequence relation does not establish exact whole RT-primer landing or a quantitative tail-recovery claim. Canonical Illumina sequence pieces are imported locally only where the selected original actually supplies the matching sequence.

🟢 The source states cBot cluster generation, linearization, 3′ end blocking and primer hybridization before on-flow-cell Klenow extension with dTTP/biotin-16-dUTP. Genome Analyzer II then performs 36 cycles of ordinary sequencing by synthesis. Subsequent fluorescent streptavidin imaging provides tail-length quantitation calibrated against standards; this fluorescence is not a base-called tail sequence. The source describes a 35-base Arabidopsis analysis branch after excluding the first low-quality cycle. The incomplete graph retains this readout limit and does not fabricate full cluster geometry or fluorescence-to-length calibration.

Deferred work includes the failed final PDF, the directly cited 3P-seq partial-RNase-T1 conditions, historical SSII/flow-cell/primer originals, exact substrate/adapter geometry and post-acquisition imaging material fates. Physical fragments and free oligos stay in the reaction until an actual separation; no perfect specificity, complete recovery or accepted missing-information conclusion is asserted.
