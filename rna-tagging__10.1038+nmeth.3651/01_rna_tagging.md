# RNA Tagging

Lapointe et al. (2015), *Protein-RNA networks revealed through covalent RNA marks*. DOI https://doi.org/10.1038/nmeth.3651 .

## Sources and scope

🟢 Associated accepted-manuscript PMC4707952 identity and Methods119–181 were read, plus manuscript supplement Figure1–2 captions and published Figure1E calibration workbook. The main is an author manuscript, not the final published original. Its operational values remain explicitly version-qualified. Published-main full procedure and complete published supplementary applicability remain pending. Standard retrieval was stopped while publisher supplements were throttled; partial external manifest/files are preserved. That interruption does not establish an access barrier or absent information. Mandatory bundle scraper inspected actual oligo contexts; no absence conclusion.

## Physical chemistry

The RBP–PUP2 fusion adds 3′ RNA uridines in vivo. This terminal addition is not an internal UMI, RNA base substitution, universal binding predictor or measured U-tag length for every molecule. The local constructor receives an explicit original RNA, U count and G/I tail, retains riboinosine as I and rejects DNA T as a G/I tail residue. The selected U-select register requires an explicitly supplied nine-residue G/I tail and at least three terminal U residues. Primer C9 and AAA cover those G/I9 and U3 positions; they are primer-derived, not reverse-transcribed tail residues. The firstDNA contains the actual source U-select primer plus a complementary copy only of upstream RNA. No priming register or efficiency is predicted.

The active transcriptome-wide route is distinct from targeted Terminator/RT-PCR/Sanger assays and separate comparison transcriptome profiling. Source yeast harvest/SDS/PCA/glass-bead lysis precedes actual bead/aqueous workup and ethanol recovery. TURBO DNase1h37°C creates DNA digestion products; GeneJetRNA cleanup follows. PolyA selection (Dynabeads mRNA) precedes yeast RiboZeroGold depletion. Historical kit internals and companion fates remain opaque, with no guessed bead/capture mechanism or immediate molecular discard.

Source G/I tailing uses yeastPAP1200U, GTP0.5mM, ITP0.15mM90min37°C, then another1200U30min37°C. PCA/chloroform/ethanol recovery precedes U-select SSIII50°C60min and85°C5min. RNaseH20min37°C retains digestion fragments until GeneJetPCR purification. Exonuclease-minus Klenow37°C30min uses an adapter-tailed random hexamer for second-strand synthesis; random priming is not a UMI. Two1:1RNA CleanXP cleanups precede GoTaqPCR, approximately nine pooled reactions, then twice0.8:1bead size selection.

## Original oligos and readout

🟢 Associated main137: U-select `GCCTTGGCACCCGAGAATTCCACCCCCCCCCAAA`; main163: second-strand `GTTCAGAGTTCTACAGTCCGACGATCNNNNNN`; main167: RP1 `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`; indexed PCR primer `CAAGCAGAAGACGGCATACGAGATXXXXXXGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. XXXXXX is a supplied six-base sample index, not a made-up study value. Structured sample-index identity is attached in the parameterized primer constructor.

Main187 prints the trimming string `TTTGGGGGGGGGTGGAATTCTCGGGTGCCAAGG`. Preserve this written analysis string separately from main137's operational primer; it is not used to silently repair the primer. Final-VOR reconciliation and original primer-to-trimming applicability remain open, rather than an established published error.

Main181 states HiSeq2500 PE50; R1 assigns genes and R2 measures terminal sequence. Three A residues supplied by the U-select primer are not independent observed RNA U-tag evidence. Full random-priming boundaries, index membership, duplex/end geometry and sequencing-primer landings remain unassembled; no finished schematic. Synthetic U-tag standard gel recovery is not imported into biological-library bead cleanup.

## Component ancestry and limits

The main G/I-tail method cites a G-tailing predecessor, explicitly adapted with inosines; the precursor was not followed in this capped pass, and no recipe import occurs. No catalogue parent edge is asserted from common Illumina primers. Target RBP changes and synthetic controls remain branches of one protocol. No second-reader absence audit or accepted missing-information claim.
