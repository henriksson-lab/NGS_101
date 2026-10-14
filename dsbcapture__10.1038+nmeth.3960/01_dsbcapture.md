# DSBCapture — research notes

Lensing et al., *Nature Methods* (2016), [doi:10.1038/nmeth.3960](https://doi.org/10.1038/nmeth.3960), PMCID: PMC5045719.

## Assay target

The sequenced material is human genomic DNA from fixed nuclei, captured at native and nuclease-induced double-strand breaks. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (modified P5/P7 Illumina duplexes and their PCR primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Full author manuscript and Methods from PMC (`article.xml`, PMC Cloud package).
- 🟢 Supplementary Figures and Supplementary Tables from the same package, read as text (`supplementary_figures.doc`, `supplementary_tables.pdf`).
- 🟡 Main Fig. 1a and Supplementary Fig. 1c, which carry the oligo set, are **raster artwork**. Only their legends are text. The bases in this repository therefore rest on a human reading of an image that cannot be re-derived from anything on disk — see "Oligo provenance" below.

## Oligo provenance

🟡 **The entire DSBCapture oligo set is image-only evidence.** `grep -oiE "[ACGT]{15,}"` over all three extracted `.txt` twins of the fetched sources returns nothing; the only bases in any of them are the two BLESS linker tags `TCGAGGTAGTA` and `TCGAGACGACG`, quoted in the article text when BLESS is described. `strings` over the raw `supplementary_figures.doc` finds no oligo either — the figure is an embedded image, not a drawn table. Supplementary Fig. 1c's legend confirms what the picture contains without repeating any of it: "Sequences of the modified P5, modified P7 and control modified P5 Illumina adapters as well as DSBCapture PCR primers (forward (PCR F) and reverse (PCR R)). AD identifies the Illumina adapter barcode sequence, three example reverse primers are shown; further primers can be created by substituting the barcode sequence."

🟢 **What the set is.** Every base of it is a canonical Illumina TruSeq part, and the duplexes and primers are internally consistent:

| drawn oligo | identity | length |
|---|---|---|
| modified P5 top | `illumina.TRUSEQ_P5_FULL` | {{= len(dsb_capture.DSB_P5) }} nt |
| modified P5 bottom | plain reverse complement of the top strand less its 3′ T | {{= len(dsb_capture.DSB_P5_BOTTOM) }} nt |
| modified P7 top | `illumina.INDEX1_PRIMER` | {{= len(dsb_capture.DSB_P7) }} nt |
| modified P7 bottom | reverse complement of the top strand plus the ligating 3′ T | {{= len(dsb_capture.DSB_P7_BOTTOM) }} nt |
| PCR F | a prefix of the modified P5 top strand | {{= len(dsb_capture.DSB_PCR_F) }} nt |
| PCR R (AD006) | `illumina.P7` + `{{= revcomp(dsb_capture.DSB_AD006) }}` (reverse complement of TruSeq index AD006 `{{= dsb_capture.DSB_AD006 }}`) + {{= len(dsb_capture.DSB_PCR_R_AD06) - len(illumina.P7) - len(dsb_capture.DSB_AD006) }} bases complementary to the modified P7 top strand | {{= len(dsb_capture.DSB_PCR_R_AD06) }} nt |

🟡 **Consequence for this repository.** A transcription error is improbable: the set is self-consistent and matches published Illumina parts exactly. But for precisely that reason, a faithful transcription of the figure and a reconstruction assembled out of canonical constants are *indistinguishable* from the sources available here. `lib/dsb_capture.py` makes that explicit by deriving the set from `lib/illumina.py` and from each duplex's own top strand instead of retyping it, so no base can drift; `SOURCE_OLIGO_EXCEPTIONS` in the same module records, per graph material, that `tools/selftest_source_oligos.py` has no text to compare against. Nothing here is a second, independent witness to the figure. Converting this evidence from 🟡 to 🟢 requires a text source — see "Open questions".

🔴 **Control modified P5 adapter.** Supplementary Fig. 1c also publishes a non-biotinylated *control* modified P5 adapter; Supplementary Fig. 1e uses it for the negative control lane "(C)", where no library is produced. It is not modelled as bases. The artwork is its only source, and no text states whether it differs from the biotinylated modified P5 by nothing but the omitted 5′ biotin, so writing it out would manufacture an identity the source does not assert. The schematic page lists it by name with that reason instead. It is a reagent of a control reaction, not of the library path, so it is not a graph material.

## Molecular interpretation

1. 🟢 DSBs in fixed nuclei are blunt-ended and A-tailed. A biotinylated, phosphorothioate-protected modified P5 duplex with a 3′ T overhang is ligated at the native break. (The Methods state the blunting kit, the Klenow exo− A-tailing and that "the modified P5 Illumina adapter, previously annealed in 1 × T4 DNA ligase buffer, was then ligated to DNA ends"; the duplex's bases are the image evidence above.)
2. 🟢 Lambda exonuclease removes excess adapter before DNA extraction. The genomic DNA is then sonicated to 200–500 bp and biotin-positive fragments are captured on streptavidin beads.
3. 🟢 The sonication-created distal end is blunt-ended, A-tailed and ligated to the modified P7 duplex. Thus P5 is uniquely attached to the native-break side and P7 to the later random end.
4. 🟢 PCR F and PCR R complete the library in a single amplification ("amplified by PCR with PCR F and PCR R primers and NEB Next PCR mix for 15 cycles"). Supplementary Figure 1d states that the first sequencing primer "has complementarity to the P5 Illumina adapter and therefore sequencing is initiated from the P5 end", so the first sequenced base is the break site. (Main Figure 1d is the DSBCapture-vs-γH2AX comparison, not the flow-cell statement.)
5. 🟡 The library is single-index. The Supplementary Fig. 1c legend establishes that PCR R carries an "Illumina adapter barcode sequence" and that further primers come from substituting it, but **no text anywhere states a barcode length, an index read or demultiplexing**. The six-base window in `lib/dsb_capture.py`, the `Index 1 (i7): 6` read length and the identity of the three drawn example primers are all measured from the Fig. 1c artwork. The schematic therefore shades the index window, captions the PCR step `INFERRED — …`, and leaves the index read's handling undetermined.
6. 🟢 EcoRV pilot libraries were sequenced 2×50 bp on MiSeq; AID-DIvA and NHEK libraries were sequenced 2×75 bp on NextSeq 500. Only proximal Read 1 was retained for break mapping; paired-end data helped distinguish PCR duplicates and compare with BLESS.

## Sequence orientation

The figures draw duplex bottom strands left-to-right as 3′→5′. The ordering list on the public page reverses those bottom-strand characters into actual 5′→3′ synthesis order. It does not copy the diagram text backwards.

`B` in the source denotes biotin, `P` a terminal phosphate and `*` a phosphorothioate bond.

## Open questions

- 🔴 **The cited detailed protocol has no URL.** Both the Methods and the BLESS comparison section point to Protocol Exchange with a literal unfilled placeholder: "A detailed protocol can be found on the protocol exchange (Insert URL)" and "the step by step method on the protocol exchange (Insert URL)". The paper as published therefore cites a document it does not identify. No source on disk names the entry, and the DSBCapture URL that the article *does* give, `http://breakome.utmb.edu`, is cited for the BLESS method it was adapted from. A located Protocol Exchange entry is the single change that would convert the whole oligo set from figure artwork to text, so it is recorded as wanted in `ref/MANIFEST.md` rather than guessed at here.
- 🔴 Whether the control modified P5 adapter is base-identical to the modified P5 adapter apart from the biotin.
- 🟡 Whether the index is six bases and whether a separate index read was acquired (see interpretation step 5).
