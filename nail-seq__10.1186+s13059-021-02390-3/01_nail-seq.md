# NAIL-seq

[Published defining article](https://doi.org/10.1186/s13059-021-02390-3), *Transcription shapes DNA replication initiation to preserve genome integrity* (2021).

🟢 Three originals inspected: published PMC8188667 main Methods169–183; published Additionalfile1 original FigureS1c/d (renderedp1, captionp2); published Additionalfile5 TableS4 resource workbook, whole source adapter/primer rows verified in the original ZIP/XML. Mandatory whole-bundle scraper ran; relevant pointers were verified against those original documents. Other fetched workbooks and the later source references remain unreviewed.

## Selected sequential-label branch

🟢 K562 cells:5µMpalbociclib36h, wash with prewarmedRPMI1640 and release2.5h,10µMEdU15min, wash, then50µMBrdU15min. GM12878 uses1µMpalbociclib and3h release; this is a variant of the same method. EdU/HU-only and mES branches are not merged into the selected two-label route.

🟢 Harvest/PBSwash, genomic DNA extraction with ProteinaseK is delegated to reference70. PurifiedgDNA undergoes100µMbiotin-PEG-azide/2mMCuSO4/4mMTHPTA/10mMsodiumascorbate25°C1h, followed by DNA precipitation. The source explicitly says copper hydrolyzes gDNA into~200–400bp fragments: the graph does not add a sonication stage. Denature95°C5min, T4PNK37°C30min. Exact donor/stock recipe and product-end populations remain unassembled; a known intended5P class is distinguished from quantitative phosphorylation of every molecule.

🟢 Anti-BrdU(BU1/75) plus ProteinG beads first isolate the BrdU fraction. Actual supernatants are denatured again and incubated with streptavidinC1 to recover the EdU fraction. All enrichment steps areRT4h. FigureS1d separately tests antibody crossreaction: EdU crossreactivity becomes undetectable after click in that assay; this is not a proof of perfect specificity or zero mixed-labelled fragments in every sequencing sample. Computational EdU-minus-BrdU subtraction is not a physical subtraction/pool.

## Source adapter inventory

🟢 TableS4 BridgeAdapter1up: `5'P-CCACGCGTGCTCTACANNNTNNNNANNNNTNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGAGT-NH2`; down: `5'-TGTAGAGCACGCGTGGNNNNNN-NH2`.

🟢 BridgeAdapter2up: `5'-ACTGGAGTTCAGACGTGTGCTCTTCCGATCT`; down: `5'-NNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGT-NH2`. The source5P and3′amines remain on the actual whole orderings. Exact local duplex relations can be computed from these strands, while full occupied ssDNA/adaptor geometry and sealed products remain unassembled. The N-rich bridge region is called barcoded in FigureS1; its study identifier purpose is left unassigned rather than automatically declared a UMI. Random6nt splint windows are not sample indices.

🟢 Source ligates both bridge types simultaneously on each bead-retained ssDNA fraction, at1µM each bridge/1600UT4DNAligase/15%PEG8000,RT≥4h. Wash ligated DNA, then tag with sourceP5-I5/P7-I7 as delegated to reference70. WholePCR/library steps and source stock quantities/cycles remain unassembled, not confirmed absent.

🟢 P5-I5: `5'-AATGATACGGCGACCACCGAGATCTACACACACTCTTTCCCTACACGACGC`; P7-I7-Index: `5'-CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGC`. Despite the P5-I5 name, its printed sequence has no variable i5 window; the model retains its fixed ordering. P7 has an explicitN6 sample-index window, without invented actual index values. SourceHiSeq2×150bp; exact primer landing/index runmask and whole library populations are not assembled.

## Deferred

Reference70 extraction/tagging protocol and corresponding historical reagent instructions are directly linked ordinary unread sources, outside the3-original cap. TableS4 contains separate ChIP-seq oligos; they are not imported into NAIL merely because they share handles. FigureS1gating/antibody assays are separate validation experiments; they are not production QC aliquots invented inside the source sequencing preparation. Known physical two-fraction partition, exact blocked adapter ends and positive source readout are retained, with no complete source/reaction or accepted absence audit.

🟡 Local end geometry: Bridge1up5′P and Bridge2up unmodifiedDNA3′OH interpretation are typed splinted interfaces adjacent to the sourceN6 windows. This retains known adapter chemistry without claiming the occupied target strand, sequence-specific six-base pairing, source-target3′ends or whole sealed reporter. Source3′amine blocks remain on both Bridge1strands and Bridge2down.
