# Stapleton synthetic long reads — two-end barcode pairing

Stapleton et al. (2016), *Haplotype-Phased Synthetic Long Reads from Short-Read Sequencing*, [10.1371/journal.pone.0147229](https://doi.org/10.1371/journal.pone.0147229), PLOS ONE 11:e0147229.

## Selected originals and scope

Three published originals: PMC4720449 main article Results/Figure 1 overview and Methods **Library preparation for synthetic long read assembly**, **Two-tube barcode pairing**, **Single-tube barcode pairing**, **Complexity determination**, **Library preparation … from mRNA**, and **Multiplexed sample preparation** (text 183–212); **S13 Table** (`s030.docx`) oligo inventory/key; and **S1 Note** (`s017.docx`) alternative barcode-pairing scope. All are published sources with assigned manifest URLs/hashes; the fetched preprint is excluded. The mandatory full-bundle oligo scraper was run and selected oligos verified against S13. Figure details and other supplements remain unreviewed, rather than confirmed absent.

One method includes genomic-DNA and mRNA front ends; these are not separate protocol-count applications. The module renders independent **DNA two-tube**, **RNA two-tube** and **DNA single-tube** source variants. They are not three preparations physically mixed in one experiment. S1 Note assigns the single-tube route to MG1655; multiplexed E. coli, mRNA and env used two-tube pairing, while other genome experiments did not use pairing. No source-specific performance estimate is turned into an invariant.

## Source preparation

🟢 Source genomic DNA is sheared to approximately 10 kb with g-TUBE, end-repaired with NEBNext, column-purified and **dT-tailed** with Klenow exo-, dTTP and PNK. This is not silently changed to the more familiar dA-tailing step. Barcode adapters are prepared separately by annealing Oligo 1 with Oligo 2 and Klenow extension, then ligated to the dT-tailed DNA. Optional gel purification precedes quantification and dilution to approximately 100,000 doubly ligated molecules. Oligo 1 contains a 16-base random molecule barcode; the separate N3 window is left functionally unassigned. Its two molecule-end instances need not have the same barcode value.

🟢 A **single primer species**, source Oligo 2, amplifies this preparation across eight LongAmp Taq reactions for the source 40-cycle program, followed by gel purification. It contains physical deoxyuracils. The model chooses the unambiguous S13 internal-dU variant for its reagent inventory and does not assert which unassigned alternative was used in a particular experiment. Other printed Oligo 2 alternatives, including `/3deoxyU/` within a longer written ordering string, are preserved as deferred source-specific notation/application questions; they are not silently converted to ddC or a nucleotide-free blocker.

🟢 USER treatment precedes fragmentase. The source calls this removal of priming regions; the model records physical digestion/cleavage products and does not reduce it to simple backbone-preserving base deletion. Fragmentase creates roughly one break **on average**, not exactly one break in every molecule. EDTA quench and 0.8× AMPure cleanup follow.

🟢 E. coli DNA ligase treatment precedes T4 DNA polymerase/Klenow treatment with **biotin-dCTP alone** for 10 min. Only then are dGTP/dTTP/dATP added for a further 15 min; the existing enzymes remain present. A 1× AMPure cleanup precedes dilute T4-ligase intramolecular circularization at 1 ng/µL for 16–48 h. T5 exonuclease digestion and column purification enrich circles; Covaris fragmentation to approximately 500 bp precedes streptavidin capture/washing. Neither perfect circle purity nor exact biotin incorporation positions are asserted.

🟢 Main bead-bound material undergoes end repair/wash, dA-tailing/wash, Oligo 3/4 adapter ligation/wash, then Oligo 5/6 lcPCR. The **PCR supernatant is retained and beads discarded**, followed by source alternatives 0.7× AMPure cleanup **or** 500–900 bp gel recovery. The alternatives are not silently combined into one recipe.

🟢 The **two-tube** route splits approximately 15% of captured beads before main end repair. This separate fraction undergoes SapI digestion/wash, Oligo 7/8 adapter ligation/wash and Oligo 6/9 lcPCR. Its supernatant is retained, beads discarded and the pairing library purified with 1.8× AMPure. The pairing library is physically mixed with main library at 1–5% of molecules for sequencing; computational barcode pairing subsequently links read groups, not DNA substrates. Oligo 9 must choose a different sequencer index from Oligo 5; unassigned N6 values do not imply that their concrete values are equal.

🟢 The **single-tube** route uses a mixture of two adapter sequencing-region variants and adds Oligo 10 during lcPCR. It has no invented 15% bead split or separate SapI branch. S1 Note explains selective amplification versus same-ended hairpin species qualitatively; full hairpin/copy geometry is not calculated here.

## mRNA front end and physical oligos

🟢 The published mRNA Methods explicitly replaces oligo-dT and TSO in its cited full-length RT procedure with Oligos 12 and 13, then enters the common source workflow at quantification. Exact cited predecessor reagents/version remain unreviewed; no current Smart-seq2 recipe is imported. Oligo 12's source T30VN is an anchored poly-dT primer with a separate molecule-barcode window. Oligo 13 retains its two rG and terminal locked-nucleic-acid G residues. These sugar and linkage annotations apply to the physical reagent; polymerase-made copies are not assigned inherited RNA/LNA or phosphorothioate chemistry by assumption. Full reverse-transcribed molecular products remain unassembled.

| Selected source input | Verbatim physical ordering notation |
|---|---|
| Oligo 2 selected variant | `{{= chem.PRIMER2.source_notation }}` |
| Oligo 3 | `{{= chem.OL3.source_notation }}` |
| Oligo 4 | `{{= chem.OL4.source_notation }}` |
| Oligo 6 | `{{= chem.OL6.source_notation }}` |
| Oligo 7 | `{{= chem.OL7.source_notation }}` |
| Oligo 8 | `{{= chem.OL8.source_notation }}` |
| Oligo 10 | `{{= chem.OL10.source_notation }}` |
| Oligo 12 | `{{= chem.RT12.source_notation }}` |
| Oligo 13 | `{{= chem.TSO13.source_notation }}` |

S13 independently defines **Index = 6-base Illumina TruSeq Small RNA multiplexing index sequence**. Its Oligo 5/9 notation uses `(Index)`. The molecule barcode, optional inline sample multiplexing index, and six-base sequencer index are distinct identifiers. Optional multiplexed adapters from S14 remain unassembled; they are not represented by relabeling the molecule barcode.

## Output, lineage and limits

🟢 Source Methods uses MiSeq, NextSeq or HiSeq with standard primer mixtures. Results reports 150×150 for the single-tube MG1655 example; that length is not assigned to every source variant. Full final molecular primer binding/read geometry is unassembled and the preview remains private despite the operational sequencing endpoint.

The publication explicitly describes modified full-length RT/Smart-seq2 ancestry for the mRNA front end only. LoopSeq later explicitly says it builds upon earlier academic work including this original, supporting broad conceptual ancestry without establishing equal adapters, enzymes or proprietary barcode-distribution mechanism. No BAsE-Seq, LoopSeq or Smart-seq constructor supplies hidden reactions.

Deferred: all full strand/end/USER-cut/biotin/circle/ligation/read geometry; applicability of remaining oligo variants; S13's other reagents and S14 multiplexing oligos; source graphics S1–S3/S12 and single-tube diagrams, cited exact full-length RT original and relevant historical kit sources; remaining linked material and independent source review. No unavailable sequence/reagent is claimed confirmed absent.
