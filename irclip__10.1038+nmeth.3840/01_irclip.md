# irCLIP: activated infrared/biotin adapter and efficient CLIP recovery

Zarnegar et al., *Nature Methods* (2016), [doi:10.1038/nmeth.3840](https://doi.org/10.1038/nmeth.3840).

The first source pass records the printed oligo inventory. A targeted local follow-up now represents the C18-bearing RT primer, circular cDNA and one conditional first-extension path. The source-described two-PCR finish is now declared at a coarse occupied-product class level. Exact full PCR populations/read geometry and applicability to the final published procedure remain unresolved.

## Sources inspected

- Accepted author-manuscript JATS associated with the published article `nmeth.3840_PMC5477425.xml.txt`: introduction and methodological comparisons 31–45; adapter preparation 53–55; crosslink/IP/nuclease/ligation 61–63; proteinase recovery 65–67; all library and custom sequencing oligos 69–73; RT, streptavidin capture, RNA degradation/elution, circularization and two PCR stages 75–77; computational processing 79–81; citations 137–159.
- Publisher supplement `nmeth.3840_41592_2016_BFnmeth3840_MOESM217_ESM.pdf`: original Fig. S1, PDF p. 1, visually inspected alongside its extracted caption (adapter chemistry and detection); Figures S2–S5 extracted captions inspected for nuclease, recovery and RT comparisons. Remaining original panels and supplementary tables are deferred.
- The NIHMS supplement PDF is separately gated. The manifest calls the publisher file “probably redundant”; that is a retrieval observation, not proof of complete document equivalence. The gated counterpart is not independently inspected or declared irrelevant.

The full bundle scraper returned eight oligo hits; the modifications and wrapped sequence strings were verified against the original Methods text. No source absence is claimed.

## Oligos, modifications and chemistry

🟢 The purchased substrate is `5′-OH-AGATCGGAAGAGCGGTTCAGAAAAAAAAAAAA/iAzideN/AAAAAAAAAAAA/3Bio/-3′`. It is DNA with an internal azide-bearing nucleotide and terminal biotin, not an RNA L3 linker. Original Methods specifies the internal modification at position 32. The local model retains the internal modified-nucleotide token without inventing a canonical base at that position.

🟢 OptiKinase phosphorylates the substrate with ATP for 1 h at 37 °C. A 5′ DNA Adenylation Kit/Mth RNA ligase activates it for 2 h at 65 °C. Strain-promoted azide/DBCO conjugation attaches IR800CW, followed by column cleanup. The resulting preadenylated adapter ligates RNA 3′ ends with T4 RNA ligase 1 overnight at 16 °C, with PEG400. The local `Adapter` state prevents adenylating an unphosphorylated substrate or repeating azide conjugation; 3′ biotin is retained.

**Source-version comparison:** publisher Fig. S1a labels the click step “1hr”; the inspected accepted-manuscript Online Methods 55 specifies 2 h at 37 °C. The private model represents that manuscript procedure. Final published procedural identity has not been reconciled, so ordinary revision differences or figure shorthand remain possible. No within-publication discrepancy or error is established.

🟢 Full printed library inventory (Methods 71–73):

| Component | As written, 5′ to 3′ |
|---|---|
| cDNA-barcode1 primer | `/5phos/WWWNNNATCACGNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGC` |
| P3short elution oligo | `CTGAACCGCTCTTCCGATCT` |
| PCR1 P3tall | `GCATTCCTGCTGAACCGCTCTTCCGATCT` |
| PCR1 P6tall | `TTTCCCCTTGTGTGTGAAGCGAAGGGTA` |
| PCR2 P3solexa | `CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT` |
| PCR2 P6solexa | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA` |
| Custom P6 sequencing primer | `CACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA` |

The source writes the two internal carbon spacers as `/iSp 18/` in the displayed primer and `/iSpl8/` in adjacent extracted prose. The normalized inventory token above is `/iSp18/`; it is not a nucleotide sequence. The phosphorylated primer's degenerate 5′ region contains the inline six-base TruSeq identifier, not necessarily a separate Illumina index read. Barcode1's printed sequence is retained; no other barcode member is selected in the model. `W` is not unrestricted `N`.

## Full original route and the local model boundary

1. 🟢 Crosslink HeLa cells with UV-C; immunoprecipitate, digest on beads and attach the activated adapter. RNase A/1 routes use PNK dephosphorylation; S1 has a separate branch where that operation is not applied. The local starting RNA scopes the RNase-A/1 route and does not conflate these alternatives.
2. 🟢 Resolve RNPs by SDS-PAGE and transfer to nitrocellulose. Proteinase K in SDS at 50 °C releases RNA-peptide complexes; phenol/chloroform and precipitation recover adapter-bearing RNA.
3. 🟢 Anneal cDNA primer; TGIRT-III extends at 60 °C for 2 h. This is source-specific primer-driven RT, not an assumed DM-tRNA-seq template-switch starter.
4. 🟢 Streptavidin beads capture RNA–cDNA hybrids through adapter biotin. P3short/MnCl₂ heating elutes cDNA/degrades RNA. CircLigase II acts for 1.5 h at 60 °C, followed by bead/isopropanol cleanup.
5. 🟢 PCR1 uses Phusion HF with P3tall/P6tall and real-time detection; reactions stop around the exponential-growth knee rather than assuming all 15 programmed cycles. Two bead size selections precede three PCR2 cycles with P3solexa/P6solexa and final cleanup.

The inspected manuscript explicitly says the carbon spacers enable linear PCR products from circularized cDNA. The local `spacer_backbone.py` now retains ordered nucleotide blocks and typed non-nucleotide C18 nodes in the actual primer and cDNA-circle state. Converting that backbone to a plain DNA sequence is rejected. The conditional copying mechanism below remains inferred; no shared implementation was changed.

## Readout, metadata and ancestry

🟢 The source names a custom P6 sequencing oligo and trims 5′ barcodes/3′ adapters during processing (73,81). The local model does not assign run cycles, separate index reads or exact library primer geometry before deriving the spacer-dependent product. No missing run-length disclosure is asserted.

Proposed metadata: RNA; bulk cultured cells; UV-C crosslink/IP; on-bead nuclease; preadenylated DNA adapter ligation; TGIRT RT; cDNA circularization; two-stage PCR; inline barcode/degenerate identifier; custom sequencing primer; published instance. Use explicit source-qualified facets and the normal publication gate; molecularly unresolved products remain visible rather than claiming full topology/readout validation.

Explicit component ancestry candidate `iclip__10.1038+nsmb.1838` → `irclip__10.1038+nmeth.3840`: original 37–43 cites König 2010 as ref. 8 (147), retains crosslink-dependent RT-stop mapping and circularized-cDNA recovery, then describes its own improved adapter, RT and purification chemistry. Scope the conceptual recovery strategy; original iCLIP's oligos, BamHI opening and buffers are not inherited. The 2015 FAST-iCLIP source, DOI `10.1261/rna.047803.114`, is also explicitly named in 43 and cited at 159 as a procedural predecessor but remains uninspected here. Shared adapter motifs alone establish no lineage.

## Source-version limit on the timing comparison

The inspected procedural text is the accepted manuscript, whereas Fig. S1a is a publisher original. Its 2 h versus 1 h timing comparison may reflect a revision difference or figure shorthand. It is not an established discrepancy within the final publication, and no publication error is asserted. The private partial model represents the inspected manuscript procedure; final published procedural identity and applicability require follow-up.


## Targeted C18-backbone modelling follow-up

This follow-up adds no new source retrieval or sections. Its source scope remains the inspected accepted-manuscript Methods 71–77; a second independent source pass and final published procedural comparison are still open.

🟢 **Printed structure/function:** the source primer has two C18 spacers around `GGATCC`, a 5′ phosphate, degenerate identifier blocks, a fixed barcode1 identifier and a terminal RT-landing sequence. The source says the spacers are required for a linear PCR product from the circular cDNA. These statements do not give an exact polymerase stopping coordinate, quantitative readthrough probability or per-cycle product distribution.

🟡 **Derived RT geometry:** the longest exact terminal suffix of the printed RT tail matches the specified adapter arm. The constructor preserves the full primer, including both C18 nodes, adds only the reverse-complemented upstream adapter interval, and retains the captured-RNA copy as a variable-length/sequence DNA block. This is a candidate sequence landing, not evidence of binding affinity or an exclusive RT product. Closing the cDNA joins its existing 3′ OH to 5′ phosphate without removing or flattening any linker.

🟡 **Conditional first extension:** an explicit `StopHypothesis` assumes Phusion extension stops before a C18 node. A P6tall primer landing on the P6-handle block initiates copying toward the original cDNA 5′ end, copies the degenerate/index blocks, crosses the circle's closure, copies the RNA-derived block in the opposite orientation, and continues through copied-adapter/RT-tail blocks until the second C18. The source circle still contains both linkers and `GGATCC`; it is not excised into a plain DNA template. The newly synthesized conditional product contains copied nucleotide blocks, with the non-nucleotide stop node explicitly retained in its traversal record rather than silently skipped.

The output preserves structured sample-index and composite degenerate-identifier features through copying. The `umi` semantic assignment for the degenerate blocks is an inference from degenerate synthesis and duplicate handling; the printed source calls this a degenerate barcode. The W pattern is retained, and captured-RNA sequence/length stay unassigned. The constructor never supplies a made-up transcript or converts C18 into a base.

The graph retains RNA alongside the nascent cDNA in a typed mixture, models streptavidin separation, and carries RNA-degradation products/P3short in the elution/circularization mixture until the source's physical cleanup. Exact residual partition is not invented. Its analytical first-copy terminal remains **unresolved**. The source-described actual PCR1/PCR2 finish is now a separate coarse product-class branch; primer competition, complete cycle-dependent products and acquired run geometry remain unassembled.

`tools/selftest_spacers.py` uses independent synthetic templates to test reverse-complement copying, partial-block intervals, closure crossing, immediate linker stops, rejection of unqualified linker traversal, rejection of linker-free rolling-circle assumptions and invalid primer/closure inputs. It also checks both source spacer tokens against the cached original when available, with normal skip behavior when that external source is absent. These checks validate the stated local relations, not the experimental polymerase-stop hypothesis or final VOR applicability.

The private `irclip_partial.html` inventory renders printed primer and inferred circle/product states; it has no complete-protocol marker. Every inferred state has a caveat, inferred caption and inferred base styling. The open-for-display circle remains a typed closed chemical backbone, not a flat DNA polymer.

## Source-qualified public preparation, 30 October 2026

The public graph preserves the accepted-manuscript authority and publisher Fig. S1a timing comparison. The adapter factory now retains the actual half-reaction aliquot/remainder, column-array recovery/pool and selected ligation contribution. ATP and covalently contributing IR800CW–DBCO are molecular structural supplies. The adapter depiction shows two nucleotide tracts of **one covalent oligo** connected through the internal `/iAzideN/` chemical site; it neither adds a canonical nucleotide there nor assumes direct DNA adjacency across that chemistry. The source calls the modification nucleotide position 32; its exact atom/sugar coordinate convention remains unassembled.

Actual source dot-blot aliquots, phase extraction, linear-polyacrylamide carrier, RNA stock selection, primer annealing, enzyme preincubation before fresh dNTP addition, retained same-tube catalyst/primer components, real-time PCR1, two bead size selections, elution directly into PCR2 primer mix, fresh master mix, three PCR2 cycles, final QC and selected loading stock are explicit. Source-prepared mixtures are decomposed into copresent functional inventory classes, never invented fresh separate doses.

The conditional P6/C18 first-copy reference remains **nonphysical analytical bookkeeping** (`changes_molecule=False`) with an unresolved terminal. It is not a second physical PCR/aliquot and its hypothetical output is **not** the actual PCR1 template. The actual source PCR1/PCR2 branch retains known whole primers with unavailable occupied-product geometry. No exact linker readthrough/stopping distribution, acquired read count/cycles, read mask or final-version equivalence is invented. The full workflow and warnings remain visible under an expandable section. Original build/selftest paths and all independent synthetic spacer tests remain.

This bounded follow-up used only the closed original source set above. It is not a complete chemistry or missing-information audit.
