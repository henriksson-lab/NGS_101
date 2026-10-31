# m6A-SAC-seq: selective allylation, cyclization and mutational readout

Hu et al., *Nature Biotechnology* (2022), [doi:10.1038/s41587-022-01243-z](https://doi.org/10.1038/s41587-022-01243-z).

This bounded pass separates the accepted-manuscript procedure from publisher-supplement inventory. The graph now represents three independently prepared source control libraries from RNA entry through adapter preparation, treatment, RT, library amplification and sequencing allocation. Whole occupied molecules and read landing remain unassembled. Final VOR procedural applicability is unverified; a source-version comparison is not a publication-error claim.

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

The current graph includes the actual upstream and downstream material stages. Exact copied products, amplification primer bodies/index choices and acquired geometry remain unassembled; source PE150 is retained.

## Publisher-SI oligo inventory, not silently assigned to manuscript states

| Component | Verbatim normalized inventory |
|---|---|
| RNA 3′ adapter | `rAppNNNNNATCACGAGATCGGAAGAGCACACGTC/iBiodT//3SpC3/` |
| cDNA 3′ adapter | `Phos-NNNNNGATCGTCGGACTGTAGAACTCTGAAC/3SpC3/` |
| RT primer | `AGACGTGTGCTCTTCCGATCT` |

`/iBiodT/` is retained as biotin-modified dT and `/3SpC3/` as a non-nucleotide terminal block. The local publisher-inventory scene preserves their identities, the printed six-base index and random window. Its UMI semantic role remains an explicit inference pending complete read-processing review. The accepted manuscript explicitly points to Supplementary Table 1 for adapter components. The graph uses that component-qualified stock inventory while preserving its publisher identity; this does not establish whole-document revision equivalence or an assigned study sample.

## Lineage, metadata and deferrals

The manuscript explicitly says its library strategy modifies **m1A-MAP** (Methods 215; R53), [doi:10.1016/j.molcel.2017.10.019](https://doi.org/10.1016/j.molcel.2017.10.019). Scope: adapter/RT/cDNA-library strategy with local changes; it does not inherit m1A detection chemistry. That predecessor was not followed in this cap. eTAM-seq later explicitly cites this method for a library subroute; its detection enzyme and RNA chemistry differ.

Proposed metadata: RNA; bulk poly(A)+/rRNA-depleted input; fragmentation/end treatment; affinity-tethered chemical labeling; covalent m6A-sensitive conversion; RNA/cDNA adapter ligation; RT/PCR; calibration and three controls; Illumina PE150; published instance. The public page is a source-qualified partial published instance with explicit unknown metadata. Deferred: final published procedure, complete version/oligo applicability, remaining SI, quantitative calibration, complete copying/PCR/read layout. The later detailed protocol and its correction are not used to revise this original study without a scoped follow-up. No confirmed missing-information claim.

## Bounded source-positive implementation after the initial model hold

The initial three/four-step component-model assessment correctly rejected adding
a sequencing label to modified RNA. A subsequent bounded implementation used
the same two closed originals and added the actual physical preparation stages;
no new source, reference, vendor recipe or whole-version equivalence was adopted.

The three independent source control libraries now include RNA entry, oligo(dT)/
RNase H/DNase I preparation, fragmentation, PNK, calibration mixture, RNA-adapter
ligation and sequential adapter digestion, primer annealing and C1 capture. FTO
remains a separate control. The model uses an initial MjDim1/allylic-SAM mixture followed by six fresh
twenty-minute mixture allocations and removed supernatants before molecular iodine, the source
printed Na2S2SO3 quench and washing. Selected reacted pathways retain residual
native/background populations; chemical adducts are never replaced by basecalls.

The actual carried RT primer is not added again. HIV RT supplies and the shorter
input RT are explicit, followed by washing, RNase H, heat eluate versus residual
beads, DNA cleanup, cDNA-adapter ligation and another cleanup. The 21 µl eluate
has distinct 1 µl qPCR, 15 µl production and nominal 5 µl untransferred-stock
fates. Actual Q5 master-mix components and unassembled historical NEBNext
functional primer classes feed production PCR; 0.8× AMPure and loading versus
stock complete the source-qualified PE150 route. No final library QC aliquot
or exact PCR cycles were invented.

The three graphs pass their normal sequencing publication gates. Full workflows
remain visible in expandable sections. Whole transcript/adduct occupancy,
mutation fractions, exact copied molecules, historical PCR primer bodies,
acquired index cycles and read landing remain unassembled. The intended RNA/
cDNA acceptor classes are explicitly qualified inferred views; adapter chemistry
is limited to free source stocks, not a fabricated mature library.

A narrow source-transcription exception retains the actual nucleotide represented
by `/iBiodT/`: the generic source parser drops that modified dT, while the model
retains its physical T. Whole modified source-order and physical-backbone checks
verify the mapping. `/3SpC3/` remains a non-nucleotide blocker, not an added C,
and is not propagated into PCR copies. Random windows retain source-pattern
semantics without declaring an unsupported UMI interpretation.

### Exact repeat-count wording and local interpretation

Accepted-manuscript `s41587-022-01243-z_PMC9378555.xml.txt` line 219,
“m6A site labeling and reverse transcription”, first describes the 50 °C
one-hour reaction, then removal of supernatant and a fresh 4+1+1+2+2 µl
mixture incubated for twenty minutes. Its parenthesis says “this labeling step
was repeated six times”. The local graph interprets this as **six total fresh
twenty-minute rounds after the initial one-hour reaction**, giving seven
mixture allocations. The prose can also be read as six repeats after the first
separately described twenty-minute round. That alternative total is unresolved
in this capped pass. The model and public page disclose this interpretation;
neither calls it a publication error nor asserts that the wording proves a
unique total. Each control is an independent source preparation; input receives
no MjDim1/I2, and no control pool or shared labeling stock is fabricated.
