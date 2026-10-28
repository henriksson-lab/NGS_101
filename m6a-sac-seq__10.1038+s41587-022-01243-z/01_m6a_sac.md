# m6A-SAC-seq: selective allylation, cyclization and mutational readout

Hu et al., *Nature Biotechnology* (2022), [doi:10.1038/s41587-022-01243-z](https://doi.org/10.1038/s41587-022-01243-z).

This bounded pass separates the accepted-manuscript procedure from publisher-supplement inventory. The graph models the chemical-control branches on prepared bead-bound RNA and stops before complete RT/library/read reconstruction. Final VOR procedural applicability is unverified; a source-version comparison is not a publication-error claim.

## Sources inspected

- Accepted-manuscript `s41587-022-01243-z_PMC9378555.xml.txt`: strategy 119–125; allyl assay 191; three controls and library preparation 209–227; predecessor reference R53. The manifest identifies this as an author manuscript.
- Publisher `s41587-022-01243-z_41587_2022_1243_MOESM1_ESM.pdf`: Supplementary Fig. 1 caption opening (biochemistry); original Table 1, PDF p. 29/printed p. 28, visually inspected for wrapped oligos, internal biotin-dT and terminal C3. Extracted table lines 446–551 were compared with the original. Remaining panels and quantitative workbooks were not systematically reviewed.
- The full bundle scraper was run after retrieval as well as during initial reading. It found the real SI oligo family and low-scoring result-table false positives. The NIHMS attachments have separate gated URLs; publisher files were fetched, but complete counterpart equivalence was not established. Short text extractions from image-heavy additional PDFs are not treated as inspected original figures.

## Chemical states and controls

🟢 MjDim1 with allylic-SAM transfers an allyl group preferentially to m6A, producing N6-allyl,N6-methyladenosine (`a6m6A`); background A can form N6-allyl-adenosine (`a6A`). Iodine then cyclizes these adducts. HIV RT produces different mutational readouts from the two chemically modified substrates. The RNA is not physically converted into canonical U/C/G at this step.

The local `AdenosinePath` retains the original A/m6A identity, requires allylation before cyclization, and refuses to flatten an adduct into a canonical RNA sequencing base. These are selected reacted pathways, not assumptions of complete conversion, fixed mutation probabilities or zero background. Residual native material and FTO erasure fractions remain unassigned.

| Control | Published role | Local scope |
|---|---|---|
| Input | Library without MjDim1 or iodine; mutation reference | No chemical labeling introduced |
| FTO− | Allylation/cyclization of target RNA | A and m6A reacted pathways retained separately |
| FTO+ | FTO pretreatment, then same allylation/cyclization | Demethylated background control; complete erasure not assumed |

The source's biochemical assay at 191 uses different labeling schedules from the biological-library procedure. The graph follows **Methods 219**: 50 °C 1 h, followed by six fresh 20-min labeling rounds, then iodine in the dark at room temperature for 1 h, quench and washing. It does not import the assay's three-round description as the library recipe. Soluble reaction products are retained until the actual bead washes; rejected washes are physical material, not automatically generated waste.

## Library route read in the manuscript

🟢 Poly(A)+ or rRNA-depleted RNA undergoes oligo(dT)/RNase H treatment, DNase I cleanup, sonication to roughly 150-nt fragments and PNK end treatment. Calibration RNA is added, followed by biotin-modified RNA 3′ adapter ligation, excess-adapter deadenylase/RecJf digestion, RT-primer annealing and C1 bead capture. FTO treatment, where applicable, precedes labeling. HIV RT acts for 3 h on labeled samples versus the source's shorter input RT. RNase H and heat elution/cleanup precede cDNA 3′ adapter ligation, Q5 amplification and final purification. Reported sequencing is HiSeq X Ten, paired-end 2×150.

The current graph begins at prepared bead-bound RNA. It does not pretend the omitted upstream/downstream operations are molecule-resolved. Exact copied products, amplification primers/index choices and acquired geometry remain deferred.

## Publisher-SI oligo inventory, not silently assigned to manuscript states

| Component | Verbatim normalized inventory |
|---|---|
| RNA 3′ adapter | `rAppNNNNNATCACGAGATCGGAAGAGCACACGTC/iBiodT//3SpC3/` |
| cDNA 3′ adapter | `Phos-NNNNNGATCGTCGGACTGTAGAACTCTGAAC/3SpC3/` |
| RT primer | `AGACGTGTGCTCTTCCGATCT` |

`/iBiodT/` is retained as biotin-modified dT and `/3SpC3/` as a non-nucleotide terminal block. The local publisher-inventory scene preserves their identities, the printed six-base index and random window. Its UMI semantic role remains an explicit inference pending complete read-processing review. The stock adapter inventory is not inserted into the accepted-manuscript branch as though complete revision equivalence had been verified, and no sample assignment is invented.

## Lineage, metadata and deferrals

The manuscript explicitly says its library strategy modifies **m1A-MAP** (Methods 215; R53), [doi:10.1016/j.molcel.2017.10.019](https://doi.org/10.1016/j.molcel.2017.10.019). Scope: adapter/RT/cDNA-library strategy with local changes; it does not inherit m1A detection chemistry. That predecessor was not followed in this cap. eTAM-seq later explicitly cites this method for a library subroute; its detection enzyme and RNA chemistry differ.

Proposed metadata: RNA; bulk poly(A)+/rRNA-depleted input; fragmentation/end treatment; affinity-tethered chemical labeling; covalent m6A-sensitive conversion; RNA/cDNA adapter ligation; RT/PCR; calibration and three controls; Illumina PE150; published instance. Keep `notes`/`pending_source_review`. Deferred: final published procedure, complete version/oligo applicability, remaining SI, quantitative calibration, complete copying/PCR/read layout. The later detailed protocol and its correction are not used to revise this original study without a scoped follow-up. No confirmed missing-information claim.
