# E/L Repli-seq — research notes

## Citation

Marchal C et al. *Genome-wide analysis of replication timing by next-generation sequencing with E/L Repli-seq.* Nature Protocols 13, 819–839 (2018). [doi:10.1038/nprot.2017.148](https://doi.org/10.1038/nprot.2017.148)

## Sources read

- 🟢 The complete defining protocol in [PMC6044726](https://pmc.ncbi.nlm.nih.gov/articles/PMC6044726/) (PMC Cloud JATS; PROCEDURE Steps 1–81 re-read 2026-10-17).
- 🟢 Its publisher Supplementary Data 1 and 2 (MOESM34) print the NEBNext adaptor, index-primer families, BrdU-IP qPCR primers and re-amplification primers. The oligo scraper was run on the original article and supplements; relevant hits were checked against those source lines.

## 1. What it is

🟢 E/L Repli-seq compares newly synthesized DNA from early- and late-S populations. It uses a BrdU pulse plus sorting and immunoprecipitation; BrdU is a selection mark, not a base conversion read directly by sequencing.

## 2. Construct-changing steps

Article PROCEDURE step numbers (PMC6044726 JATS lines 306–601):

1. 🟢 Add BrdU to 100 μM and incubate 2 h at 37 °C (Steps 1–2).
2. 🟢 Fix in 75% ethanol, stain with PI/RNase A and sort 120,000 early- and 120,000 late-S cells (Steps 3–18).
3. 🟢 Lyse in SDS-PK buffer (0.2 mg/ml proteinase K, 56 °C 2 h) and purify on a Zymo Quick-DNA Microprep column (Steps 19–24).
4. 🟢 Covaris sonication to ~200 bp (175 W, 10% duty, 200 cycles/burst, 120 s; Steps 25–28).
5. 🟢 One NEBNext Ultra End Prep reaction (End Prep Enzyme Mix, 20 °C 30 min then 65 °C 30 min; Steps 30–32). 🟡 The single 3′ dA is the kit's dA-tailing design; the article does not print it.
6. 🟢 Ligate the NEBNext Adaptor (0.45 μM; diluted 1:10 below 100 ng) with Blunt/TA Ligase Master Mix plus Ligation Enhancer, 20 °C 15 min (Steps 33–34). 🟡 The adaptor's 12-bp terminal stem and 3′-T overhang are computed from the printed oligo (`lib/illumina.NEBNEXT_HAIRPIN_MODEL`, checked equal to the Supplementary Data 1 transcription at construction).
7. 🟢 Add 3 μl USER, 37 °C 15 min (Steps 35–36), opening the deoxyU loop; then DNA Clean & Concentrator-5 (Step 37).
8. 🟢 Library construction precedes the IP. Denature at 95 °C 5 min, ice 2 min; IP with mouse anti-BrdU then rabbit anti-mouse IgG, pellet at 16,000 g without beads, wash in IP buffer; the supernatant is removed (Steps 38–51).
9. 🟢 Proteinase K overnight at 37 °C plus 60 min at 56 °C; DNA Clean & Concentrator-5 into 16 μl (Steps 52–55).
10. 🟢 Optional adaptor qPCR (Step 56). Dual-index PCR with NEBNext Q5 Hot Start HiFi Master Mix and one unique i5/i7 pair per library: 98 °C 30 s; 14 cycles of 98 °C 10 s / 65 °C 75 s; 65 °C 5 min (Steps 57–58). Box 1 re-amplification uses TS-Oligo 1/2 when yield is low.
11. 🟢 AMPure XP 90 μl per 100 μl PCR (0.9×), three 80% ethanol washes, elute in 33 μl (Steps 60–72); QC, pool at 10 nM (Steps 73–80).

## 3. Final library and reads

🟢 These are conventional dual-indexed Illumina genomic libraries. Step 81: sequence the pool on a HiSeq; "Generally, 50bp single-end reads is sufficient", and longer reads may be needed for SNP phasing. The analysis pipeline accepts single- or paired-end data. The model records Read 1 = 50 cycles plus both index reads. 🟡 Both index reads are required by the unique-i5/i7 design (Step 57); their acquired cycle counts are not printed, although the Supplementary Data 1 indexes are 8 nt.

🟢 Supplementary Data 1 prints the NEBNext Adaptor for Illumina as `5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC/ideoxyU/ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`, with an internal deoxyU and terminal phosphorothioate bond. It lists i5 (i501–i508) and i7 (i701–i712) index-primer families. Supplementary Methods 1 and Step 56 print BrdU-IP adaptor qPCR primers `ACACTCTTTCCCTACACGACGC` and `GACTGGAGTTCAGACGTGTGC`; Supplementary Data 2 prints the two short re-amplification oligos.

🟡 The page uses a canonical TruSeq-type endpoint for the final library. Hairpin-ligated, USER-opened and IP intermediates are drawn as labelled regions rather than base-resolved duplexes.

## Bounded pass 2026-10-17 (round breadth-first-09)

Corrected the graph from the article's PROCEDURE: combined End Prep reaction, hairpin adaptor end geometry, the missing USER step and column cleanups, both denatured strands entering the IP (the unlabelled strand was previously discarded at denaturation), proteinase K release, the Q5 14-cycle PCR, AMPure cleanup, and 50-bp single-end reads (Read 2 no longer marked as acquired). Open questions: acquired index-read cycles; the NEB E7370 manual (End Prep composition and adaptor dU geometry beyond the printed oligo) was not read; Box 2 / Ryba et al. enrichment primer sets were not followed.

## 4. Boundary

🟢 This page does not merge high-resolution multi-fraction Repli-seq or copy-number-only Sort-seq into E/L Repli-seq.

## Native ownership migration, 2026-10-18

Read the published article JATS Procedure30–81 (library preparation, USER/column cleanup, BrdU IP, indexing, final purification and sequencing) and published supplementary PDF Methods1/Data1/Data2. Data1 primer columns were also inspected with layout-preserving extraction of original PDF pages9–10: physical i501–i508 and i701–i712 index families are8nt; the adjacent Expected Index Read column is not extra ordered sequence. Scraper hits from superseded preprints were excluded. No vendor original or predecessor was added to this capped two-original pass. Earlier pulse/FACS/extraction conditions retain their existing source scope.

Removed all conformation_replication_specs and batch_ngs constructors, including the shared NascentPulse state that concatenated parental and nascent strands. Local states keep BrdU-containing nascent and unlabelled strands separate. End Prep remains one reaction; the single dA and hairpin fold are derived, not established native intermediate geometry. Both published adapter ordering modifications and the computed12bp/3′T fold remain explicit. USER excision products are molecular inputs to Step37 cleanup; only its resulting column flowthrough is a physical discard. Denaturation creates distinct single-strand Scenes, and the BrdU strand propagates through IP/purification into PCR without using the expected completed library as input.

The local endpoint uses the source primer-family scaffolds and8nt physical index windows, with no chosen index members or inferred acquired index cycles. Early and late remain separate terminal sequencing plans selecting source SR50; optional longer/paired reads, qPCR optimization and re-amplification remain source alternatives, not invented molecular steps. Complete folded ligation/USER/first-cycle PCR geometry, individual index assignments and optional QC/re-amplification graphs remain unresolved implementation targets. No publication-absence claim is accepted.
