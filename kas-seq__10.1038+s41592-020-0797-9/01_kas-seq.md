# KAS-seq: enrich DNA carrying a reversible guanine adduct

Wu et al. (2020), [Kethoxal-assisted single-stranded DNA sequencing captures global transcription dynamics and enhancer activity in situ](https://doi.org/10.1038/s41592-020-0797-9).

## Sources and identity

First capped pass, one reader, 28 October 2026, three originals: the article-associated accepted manuscript PMC7205578 (Results chemical principle; Methods bulk and low-input preparations, sequencing and references17/43); **published publisher Supplementary Information MOESM1**, complete Supplementary Protocol sections1–3; [Author Correction](https://doi.org/10.1038/s41592-020-0881-1), published31May2020, complete notice read on the primary publisher record. The publisher article identity/abstract was checked; its complete version-of-record Methods was not obtained in this pass. The PMC manifest explicitly says author manuscript and its text retains the corrected error. It is not silently treated as the corrected published PDF.

The mandatory oligo scraper was run on the complete bundle. Its only hit, `CHRNA10` in a gene table, is a false positive. That observation establishes no sequence-absence claim. The published Supplementary Protocol has complete extracted prose; it does not present an image-only oligo table requiring reconstruction. Publisher MOESM2 reporting summary and MOESM4 source-data figures remain outside this three-original chemistry pass.

## Correction and version-qualified differences

🟢 The correction replaces **100mM with100µM** for the synthetic DNA oligo in the *in vitro* characterization Methods. The archived manuscript line106 retains100mM. This is a source correction to a validation assay, not a change to the5mM live-cell N3-kethoxal labeling condition. The model uses the latter and never propagates the incorrect synthetic-oligo concentration.

The schematic follows the **published Supplementary Protocol**, whose conditions differ from the accepted-manuscript main Methods in several places. These are version-qualified alternatives, not assertions of additional publication errors:

| Component | Publisher Supplementary Protocol | Accepted-manuscript main Methods |
|---|---|---|
| Bulk click input | 2µg total DNA;10µL10×PBS | 1µg genomic DNA |
| Bulk RNaseA treatment | 15min | 5min |
| Low-input extraction | Quick-DNA Microprep Plus,D4074 | Quick gDNA mini plus,D4068 |
| Low-input Tn5 identifier | 15027865; buffer15027866 | 10527865 |
| Low-input input fraction | 5µL of55µL;50µL captured | No separate volume assignment in main paragraph |

## Molecular meaning

🟢 Article Results: N3-kethoxal reacts with guanine N1/N2 positions exposed in single-stranded nucleic acid. Watson–Crick pairing protects these positions. The installed azide supports DBCO-PEG4-biotin click coupling. RNaseA removes RNA contamination before DNA recovery; biotin-tagged DNA is enriched on streptavidin. Heating95°C removes the reversible kethoxal label and permits DNA recovery/amplification.

🟡 The model shows a **representative accessible guanine**, not a known genomic coordinate, a claim that every genomic G is accessible, or complete conversion of every accessible molecule. Adduct identity and genomic-site assignment are separate fields of the local state. No covalent atom mapping, reaction kinetics, capture efficiency or polymerase-bypass mechanism is claimed checked by drawing the tag.

A captured fragment can contain duplex and single-stranded intervals. Sequencing reads its genomic DNA sequence, not a chemically encoded per-base KAS signal. A peak is fragment-level enrichment and is not proof of an individual labeled base or PolII-specific origin: other sources of exposed DNA can contribute.

## Published bulk route

1. Label1–5million cells with5mM N3-kethoxal for10min37°C. Isolate DNA with PureLink; borate buffer preserves the adduct. The mouse-liver front end homogenizes tissue, removes large pieces and collects cells before the same labeling/preparation; no additional library chemistry is assumed.
2. Click2µg DNA with DBCO-PEG4-biotin at37°C90min; add RNaseA15min37°C. DNA Clean & Concentrator purification follows.
3. Sonicate1µg biotinylated DNA to150–350bp. Reserve5µL of100µL as input; enrich95µL on10µL streptavidinC1 beads15min at room temperature. Wash five times.
4. Heat beads95°C10min to elute. Heat the input in parallel and immediately ice both. Both fractions proceed to **Swift Accel-NGS Methyl-seq30024** library preparation. The kit name does not imply bisulfite conversion of KAS samples. No unreported conversion, Adaptase tail, adapter bases or Swift formulation is supplied from another protocol.
5. Supplement:7–8PCR cycles for input,12–14for enriched DNA; library200–500bp; single-end Illumina sequencing. Associated manuscript specifies NextSeq500 single-end80bp, with30million reads per library. The exact run mask/index assignments and run-primer molecules remain open.

## Published low-input route

1. Label1000/5000/10000cells as above; extract using the supplement-specific D4074 kit. Click biotin in half-volume; RNaseA treatment is5min, then cleanup.
2. Tagment biotinylated DNA with source-named Tagmentation DNA Enzyme15027865 and Buffer15027866 at37°C30min. Enzyme volumes1.5/2/5µL correspond to1000/5000/10000cells. Recover55µL after DNA cleanup.
3. Reserve5µL input; enrich50µL using5µL streptavidinC1 beads and five washes. This is **not a5% split**. The source carries bead-bound enriched DNA directly to PCR; do not insert the bulk heat-elution step here.
4. Prepare20µL input or bead suspension,2.5µL each i5/i7 primer from Illumina20027213,25µL UltraIIQ5 mix.72°C5min precedes95°C10min and15cycles98°C10s/60°C30s/72°C1min. Beads are removed after PCR; MinElute purifies200–600bp libraries; sequencing is single-end. Do not infer exact index/loaded-adapter oligos from catalogue numbers or impose bulk Swift chemistry.

## Open source and model boundaries

Exact bulk kit architecture, low-input loaded adapters and selected i5/i7 members are unassembled in this pass. The source-named kits and PCR operations are retained with explicit limitations rather than arbitrary canonical substitutions. No oligo/reaction detail is asserted absent. Historical vendor originals, the complete corrected article, linked [Protocol Exchange version1](https://doi.org/10.21203/rs.2.24141/v1) and an independent second reader are deferred.

Article reference17 explicitly supplies **N3-kethoxal synthesis** from Weng et al., *Keth-seq for transcriptome-wide RNA structure mapping*, DOI10.1038/s41589-019-0459-3. This supports reagent-synthesis ancestry, not inheritance of Keth-seq RNA library chemistry. The KAS source establishes its own DNA capture route; ChIP-seq, ATAC-seq and GRO-seq comparisons do not make those assay libraries its parents.

## Public preparation and inventory wiring

The public page retains the published-protocol/accepted-manuscript differences and correction above. N3-kethoxal and DBCO-biotin are molecular contributors, with conditions declared separately. Actual selected DNA contributions preserve stock remainders. C1 bead prewashing and bulk residual beads are explicit. Low-input PCR components derive from each actual25µlQ5 master-mix allocation; bulk Swift supply classes remain source-qualified and unassembled. Input and enriched library loading portions leave separate prepared stocks. No historical kit sequence or whole polymer architecture is supplied by default.
