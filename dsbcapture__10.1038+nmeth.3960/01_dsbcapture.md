# DSBCapture — research notes

Lensing et al., *Nature Methods* (2016), [doi:10.1038/nmeth.3960](https://doi.org/10.1038/nmeth.3960), PMCID: PMC5045719.

## Assay target

The sequenced material is human genomic DNA from fixed nuclei, captured at native and nuclease-induced double-strand breaks. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (modified P5/P7 Illumina duplexes and their PCR primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Full author manuscript and Methods from PMC (`article.xml`, PMC Cloud package).
- 🟢 Supplementary Figures and Supplementary Tables from the same package, read as text (`supplementary_figures.doc`, `supplementary_tables.pdf`).
- 🟢 The original authors’ [Protocol Exchange detailed method](https://doi.org/10.1038/protex.2016.52), PDF pp. 7–9 (`supplement_0.pdf.txt`, lines 568–648), prints the modified P5/P7/control strands, PCR primers, and six-base single-index design as text.

## Oligo provenance

🟢 **The original detailed method resolves figure-only evidence.** The paper’s Fig. 1a and Supplementary Fig. 1c are raster artwork, but Lensing et al.’s Protocol Exchange method (doi:10.1038/protex.2016.52), PDF p. 8–9, prints the adapter duplexes, control strand and PCR primer family. Local extracted text lines 608–648 can be checked base by base. The method also states “single index, 6 bp” at line 574.

🟢 **What the set is.** Every base of it is a canonical Illumina TruSeq part, and the duplexes and primers are internally consistent:

| drawn oligo | identity | length |
|---|---|---|
| modified P5 top | `illumina.TRUSEQ_P5_FULL` | {{= len(illumina.TRUSEQ_P5_FULL) }} nt |
| modified P5 bottom | plain reverse complement of the top strand less its 3′ T | {{= len(illumina.TRUSEQ_P5_FULL) - 1 }} nt |
| modified P7 top | `illumina.INDEX1_PRIMER` | {{= len(illumina.INDEX1_PRIMER) }} nt |
| modified P7 bottom | reverse complement of the top strand plus the ligating 3′ T | {{= len(illumina.INDEX1_PRIMER) + 1 }} nt |
| PCR F | a prefix of the modified P5 top strand | {{= 44 }} nt |
| PCR R (AD006) | `illumina.P7` + `{{= revcomp("GCCAAT") }}` (reverse complement of TruSeq index AD006 `{{= "GCCAAT" }}`) + {{= 21 }} bases complementary to the modified P7 top strand | {{= len(illumina.P7) + 6 + 21 }} nt |

🟢 **Control modified P5 adapter.** The detailed method prints the control strand without the 5′ biotin, otherwise base-identical to the modified P5 top strand, and shares the same phosphorylated complement (lines 612–617). This is a negative-control reagent; no product library is made on that branch.

## Molecular interpretation

1. 🟢 DSBs in fixed nuclei are blunt-ended and A-tailed. A biotinylated, phosphorothioate-protected modified P5 duplex with a 3′ T overhang is ligated at the native break. (The Methods state the blunting kit, the Klenow exo− A-tailing and that "the modified P5 Illumina adapter, previously annealed in 1 × T4 DNA ligase buffer, was then ligated to DNA ends"; the duplex bases are printed in the Protocol Exchange oligonucleotide list; the Methods anneal each duplex at 10 µM in 1 × T4 DNA ligase buffer, 95 °C 10 min, then slow cooling to room temperature.) In situ conditions (Protocol Exchange steps 42–63): Quick Blunting Kit 25 °C 45 min; Klenow exo− + dATP 37 °C 45 min; T4 DNA ligase 16 °C 15–20 h, with a parallel no-ligase control.
2. 🟢 Lambda exonuclease (50 U, 37 °C 30 min) removes excess adapter before DNA extraction. Nuclei are then lysed with proteinase K (55 °C 30 min, 65 °C 30 min) and DNA is isopropanol-precipitated (steps 76–102). The genomic DNA is sonicated to 200–500 bp (Bioruptor HIGH, 30 s on/off, ~25 cycles) and 50 µg is bound to 5 µl MyOne Streptavidin C1 beads at 4 °C for 45 min (steps 103–112).
3. 🟢 The sonication-created distal end is blunt-ended, A-tailed and ligated to the modified P7 duplex. Thus P5 is uniquely attached to the native-break side and P7 to the later random end.
4. 🟢 PCR F and PCR R complete the library in a single on-bead amplification ("amplified by PCR with PCR F and PCR R primers and NEB Next PCR mix for 15 cycles", 65 °C annealing). The beads are pelleted and discarded, and the supernatant is purified with MinElute and size-selected at 250–1200 bp on a BluePippin (Protocol Exchange steps 124–128 give a 2 % gel-excision alternative). Supplementary Figure 1d states that the first sequencing primer "has complementarity to the P5 Illumina adapter and therefore sequencing is initiated from the P5 end", so the first sequenced base is the break site. (Main Figure 1d is the DSBCapture-vs-γH2AX comparison, not the flow-cell statement.)
5. 🟢 The original detailed method explicitly calls the library TruSeq single-index with a six-base index (`supplement_0.pdf.txt`, line 574), prints eight PCR R examples (lines 632–648), and describes barcode-based multiplexing. It does not specify the instrument’s separate Index 1 acquisition settings.
6. 🟢 EcoRV pilot libraries were sequenced 2×50 bp on MiSeq; AID-DIvA and NHEK libraries were sequenced 2×75 bp on NextSeq 500. Only proximal Read 1 was retained for break mapping; paired-end data helped distinguish PCR duplicates and compare with BLESS.

## Sequence orientation

The figures draw duplex bottom strands left-to-right as 3′→5′. The ordering list on the public page reverses those bottom-strand characters into actual 5′→3′ synthesis order. It does not copy the diagram text backwards.

`B` in the source denotes biotin, `P` a terminal phosphate and `*` a phosphorothioate bond.

## Open questions

- 🟢 The cited Protocol Exchange entry is now located: doi:10.1038/protex.2016.52; its original PDF is available via the source manifest and prints the complete oligo set.
- 🔴 The source specifies a six-base single index but does not give the instrument’s separate index-read settings.
- 🟢 Source detail, not a discrepancy: the printed AD19 PCR R carries eight bases (`CGTTTCAC`) between P7 and the Read 2 arm, while AD01–AD12 carry six. The model draws the AD06 example only.

## Review trail

- 2026-10-17 bounded pass (round 2026-10-17-breadth-first-09): re-read the Online Methods and Protocol Exchange steps 42–131 plus its oligo list. The graph now names the source reagents and conditions. It adds the cross-link reversal and DNA extraction step, the digested-adapter waste and the on-bead PCR with bead discard. It also adds MinElute/BluePippin 250–1200 bp size selection. The PCR finish is now local rather than the shared `_finish_pcr` helper.

## Source-owned validation migration, 2026-10-18

Local endpoint now derives the six-base AD06 example directly from the declared modifiedP5/P7 and printed PCRR, without shared break-library or TruSeq defaults. Added paired source adapter scenes and native/distal dA-tail states; distal preparation adds an overhang at the sonication-created end only. Exact PCR geometry now uses the preceding modifiedP5/P7 ligation product with the printed44nt PCRF and AD06 PCRR, rather than the final indexed library. Existing physical waste, bead disposal and size-selection steps are retained. Graph selects the2×75 NextSeq acquisition branch; molecular six-base index content no longer creates an acquired index read.

Source scope: Published article PMC5045719 Online Methods sequencing passage (2×50 MiSeq EcoRV pilot/2×75 NextSeq500); original Protocol Exchange supplement_0.pdf.txt steps124–131 lines568–605 and oligo list lines610–649.

Deferred: Other seven AD primers and2×50 EcoRV pilot remain additional declared variants rather than selected graph runs; Complete native-break repair, sonication and ligation input-output relations remain unresolved. This bounded check establishes neither complete chemistry validation nor an accepted missing-information claim.
