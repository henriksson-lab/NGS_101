# FLAM-seq — full-length RNA and retained poly(A)

Legnini and colleagues, *FLAM-seq: full-length mRNA sequencing reveals principles of poly(A) tail length control*, Nature Methods (2019), DOI [10.1038/s41592-019-0503-y](https://doi.org/10.1038/s41592-019-0503-y).

Capped source-positive collection, private incomplete graph; no public promotion or complete-source claim. The published supplementary protocol MOESM3 was read through its preparative, SMRTbell and sequencing setup sections. The published supplementary-information MOESM1 pp1–2 was checked for the specific anchoring/readout limitation. The original protocol p2 was rendered for its oligos after mandatory bundle scraping. The defining final main text was not obtained; Nature returned an inaccessible identity-provider redirect, and the standard collector recorded no open full text. Its preprint is superseded and was not used as chemistry authority.

🟢 The published protocol positively establishes two-round TruSeq poly(A) selection, G/I tailing using the named USB tail-length assay mix, bead cleanup, C-rich anchored reverse transcription with SMARTScribe and an iso-TSO, DNA bead recovery, Advantage 2 amplification, DNA cleanup, SMRTbell 1.0 SPv3 damage/end repair, adapter ligation, Exo III/VII treatment and two DNA cleanups. RNA and its RT primer are heated before the enzyme/TSO master mix is added at 42 °C. Sequencing setup separately conditions the v3 sequencing-primer stock at 80 °C, then anneals it to the template at 20 °C, binds Sequel polymerase 3.0, and cleans the bound complex using Magbead buffer with an explicit instruction to omit ethanol washing. Sequel SMRT Cell 1M v3 and sequencing kit 3.0 are named; movies span 360–600 minutes with diffusion loading, and CCS processing uses at least three passes. The recovered preparation is assayed by Qubit, then mixed with an independent diluted internal control and loading reagents before loading. These are source-specific historical components, not current-kit defaults.

🟢 Selected oligos as written in the published protocol:

| Oligo | As written |
|---|---|
| RT primer 1 | `GGTAATACGACTCACTATAGCGAGANNNNNNNNNNCCCCCCCCCTTT` |
| RT primer 2, alternative | `TGAGTCGGCAGAGAACTGGCGAANNNNNNNNNNCCCCCCCCCTTT` |
| Template-switch oligo | `iCiGiCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG` |
| PCR primer 1, paired with RT 1 | `GGTAATACGACTCACTATAGCGAG` |
| PCR primer 2, paired with RT 2 | `TGAGTCGGCAGAGAACTGGCGAA` |

The source describes the `i` bases as stereoisomers of dC/dG; the physical noncanonical bases remain distinct from ordinary CGC. Terminal `rGrGrG` is RNA. The N10 region is retained as a degenerate source pattern; its identifier role is not assigned by this selected implementation. The protocol also names 5′ PCR Primer II A in the reaction; the exact kit oligo/landing is not inferred from similarity to a template-switch body.

🟡 For a selected canonical poly(A)-ending RNA and an explicitly resolved primer, the helper extends an anchored `C9TTT` primer by copying the remaining native tail and upstream RNA. The three anchoring Ts already occupy the final three native As; the helper neither double-counts them nor copies an invented entire artificial G/I tail. It does not model inosine pairing, switch-junction chemistry, kit composition, full molecule recovery or yield.

🟢 The published SI explicitly says terminal G/I addition plus the anchored three Ts excludes native terminal modifications from the selected library, and native guanylation cannot be distinguished from added G/I. Non-anchored alternatives are discussed but not implemented as the default. Thus a full-transcript/poly(A) preparation does not establish complete native 3′ end modification coverage.

RNA and other molecular companions are not destroyed by labeling a product “cDNA”: no RNA hydrolysis step is invented. Enzymatic digestion produces fragments before actual cleanup. Analytical aliquots remain separate from production material. Full GI enzyme/nucleotide identity, stochastic tailing, switch/copying geometry, PCR landing, SMRTbell adapter bases and sequencing-primer geometry remain unassembled; actual named historical kit originals and the linked PacBio loading reference are deferred. No missing component is claimed confirmed absent.
